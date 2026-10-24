"""Simulation experiment of the study Huet1983."""

from typing import ClassVar, override

from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan, Simulation

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class Huet1983(ICGSimulationExperiment):
    """Simulation experiment of Huet1983.

    Blood clearance of indocyanine green after an intravenous bolus of
    0.5 mg/kg against the hepatic blood flow (Q, from the clearance and the
    extraction ratio of indocyanine green) of 53 patients with cirrhosis
    (69.5 kg) and 12 patients with chronic hepatitis (66.3 kg) (Fig 3B), with
    scans of the hepatic blood flow of the model (`Qh` per body weight) for
    every stage of cirrhosis; the control is the curve of chronic hepatitis.
    The normal volunteers of Huet1980 are shown at the apparent hepatic blood
    flow of lidocaine, a different measure of the flow.

    The study fits one exponential to the arterial samples up to 15 min and
    gives the systemic blood clearance, kel * D/C0 / (1 - Hct) (`cl_kel` /
    (1 - `HCT`), the window 2 - 15 min).
    """

    #: the body weight of the groups [kg], Table 1
    bodyweights: ClassVar[dict[str, float]] = {"hepatitis": 66.3, "cirrhosis": 69.5}

    def _bodyweight(self, key: str) -> float:
        """Get the body weight of a key of `cirrhosis_map` [kg]."""
        group = "hepatitis" if key == "Control" else "cirrhosis"
        return self.bodyweights[group]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig3B")
        for label in df.label.unique():
            # Fig3B_<group>_<measurement>
            dsets[self.to_sid(label)] = DataSet.from_df(
                df[df.label == label], self.ureg
            )

        # the normal volunteers of Huet1980, at the apparent hepatic blood flow
        # of lidocaine
        frames = self.load_dataframes_by_substance("Tab2", study="Huet1980")
        for substance, measurement, key in [
            ("lidocaine", "hbf", "Q"),
            ("indocyanine green", "clearance", "ICGclearance"),
        ]:
            df = frames[substance]
            dset = DataSet.from_df(df[df.measurement_type == measurement], self.ureg)
            dsets[f"Huet1980_Tab2_{key}_controls"] = dset[dset.group == "controls"]
        return dsets

    @override
    def simulations(self) -> dict[str, Scan]:
        scans = {}
        for key, value in self.cirrhosis_map.items():
            bodyweight = self._bodyweight(key)
            scans[self.to_sid(f"scan_bloodflow_{key}")] = Scan(
                simulation=Simulation(
                    time_unit="min",
                    end=30,
                    preinit_changes={
                        **self.default_changes(),
                        "BW": Q(bodyweight, "kg"),
                        "IVDOSE_icg": Q(0.5 * bodyweight, "mg"),
                        "f_tissue_loss": Q(value, "dimensionless"),
                        "f_shunts": Q(value, "dimensionless"),
                    },
                ),
                dimensions=[
                    Dimension(
                        "dim_bloodflow", values=self.bloodflow_scan_changes(stop=3.5)
                    ),
                ],
            )
        return scans

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk()
        return {"fig_cl": self._figure_clearance_vs_bloodflow()}

    def _figure_clearance_vs_bloodflow(self) -> FigureMPL:
        """Figure of the blood clearance against the hepatic blood flow per body weight.

        The dashed line is the hepatic blood flow of the model per body weight.
        """
        f, ax = base_figure()
        unit = "ml/min/kg"
        ax.set_xlabel(f"Hepatic blood flow [{unit}]", fontdict=self.font)
        ax.set_ylabel(f"{self.pk_labels['cl_blood']} [{unit}]", fontdict=self.font)
        ax.axvline(
            x=(self.hepatic_bloodflow() / Q(75, "kg")).to(unit).magnitude,
            color="grey",
            linestyle="--",
        )

        for key in self.cirrhosis_map:
            sim_key = self.to_sid(f"scan_bloodflow_{key}")
            df = self.pk_dfs[sim_key]
            xres = self.results[f"task_{sim_key}"]
            bodyweight = Q(self._bodyweight(key), "kg")
            hct = float(self.scan_values(xres, "HCT")[0])
            # the blood clearance per body weight
            factor = 1 / ((1 - hct) * bodyweight.magnitude)
            x = self.scan_values(xres, "Qh", "ml/min") / bodyweight.magnitude
            values, lower, upper = self.pk_band(df, "cl_kel", "ml/min")
            self.plot_band(
                ax,
                x,
                values * factor,
                None if lower is None else lower * factor,
                None if upper is None else upper * factor,
                color=self.cirrhosis_colors[key],
                label="Control (chronic hepatitis)" if key == "Control" else key,
            )

        # data
        for group, color, fill in [
            ("cirrhosis", "purple", "purple"),
            ("hepatitis", "black", "white"),
        ]:
            ax.plot(
                self._datasets[f"Fig3B_{group}_hbf"].get_quantity("value").magnitude,
                self._datasets[f"Fig3B_{group}_cl"].get_quantity("value").magnitude,
                marker="s",
                linestyle="",
                color=color,
                markerfacecolor=fill,
                label=f"Huet 1983 {group}",
            )
        ax.plot(
            self._datasets["Huet1980_Tab2_Q_controls"].get_quantity("value").magnitude,
            self._datasets["Huet1980_Tab2_ICGclearance_controls"]
            .get_quantity("value")
            .magnitude,
            marker="^",
            linestyle="",
            color="grey",
            label="Huet 1980 normals\n(Q of lidocaine)",
        )

        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0, top=35)
        ax.legend(
            fontsize=self.legend_font_size - 2,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
        )
        return f


if __name__ == "__main__":
    MODEL.run_experiments(Huet1983, output_dir=Huet1983.__name__)
