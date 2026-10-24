"""Simulation experiment of the study Gadano1997."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    DataType,
    Dosing,
    Fasting,
    Health,
    PKPDData,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Change, Dimension, Scan, Simulation
from sbmlsim.units import Quantity

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class Gadano1997(ICGSimulationExperiment):
    """Simulation experiment of Gadano1997.

    Extraction ratio and clearance of indocyanine green at the steady state of
    an infusion (transjugular catheterization, peripheral and hepatic venous
    samples 20 - 24 min after the start) in 59 controls, 36 patients with
    hepatic fibrosis and 356 patients with cirrhosis, and the dependency of the
    clearance on the hepatic blood flow, simulated as a scan of the blood flow.

    The controls are patients with chronic hepatitis B or C and a normal liver
    architecture (no fibrosis in the biopsy), not healthy subjects. Every
    patient was sedated with meperidine (pethidine) 50 mg i.m. The infusion
    was 4 mg/min for 6 min, then 0.8 mg/min in the controls and in fibrosis,
    and 4 mg/min for 3 min, then 0.4 mg/min in cirrhosis. The plasma clearance
    is the infusion rate over the peripheral plasma concentration
    (`CLinfusion_icg`), the hepatic blood flow HBF = CL / E / (1 - Hct).
    Fibrosis has no model counterpart, its data is shown without a curve.
    """

    groups: ClassVar[list[str]] = ["Controls", "Hepatic Fibrosis", "Cirrhosis"]

    group_colors: ClassVar[dict[str, str]] = {
        "Controls": "black",
        "Hepatic Fibrosis": "tab:orange",
        "Cirrhosis": "purple",
    }

    #: the infusion of the controls and of fibrosis: rate 1 [mg/min], its end
    #: [min], rate 2 [mg/min]
    infusion_control: ClassVar[tuple[float, float, float]] = (4.0, 6.0, 0.8)
    #: the infusion of cirrhosis: rate 1 [mg/min], its end [min], rate 2 [mg/min]
    infusion_cirrhosis: ClassVar[tuple[float, float, float]] = (4.0, 3.0, 0.4)

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Tab3")
        for group in self.groups:
            # a data set holds the rows of one measurement, each in its unit
            for measurement, key in [("extraction_ratio", "er"), ("clearance", "cl")]:
                df_subset = df[(df.measurement == measurement) & (df.group == group)]
                dset = DataSet.from_df(df_subset, self.ureg)
                if measurement == "extraction_ratio":
                    dset.unit_conversion("mean", factor=1 / Q(100, "percent"))
                dsets[self.to_sid(f"Tab3_{key}_{group}")] = dset

        df = self.load_dataframe("Fig1")
        for group in self.groups:
            dset = DataSet.from_df(df, self.ureg)
            dsets[self.to_sid(f"Fig1_hbfvsCl_{group}")] = dset[dset.group == group]
        return dsets

    @property
    def cirrhosis_keys(self) -> list[str]:
        """The stages of cirrhosis of `cirrhosis_map`, without the control."""
        return [key for key in self.cirrhosis_map if key != "Control"]

    def _infusion(self, key: str, end: float) -> Simulation:
        """Get the infusion of indocyanine green of the control or of a stage.

        Args:
            key: the key of `cirrhosis_map`, the control has the infusion of the
                controls, a stage of cirrhosis the infusion of cirrhosis.
            end: end of the simulation [min].
        """
        changes: dict[str, Quantity | float] = dict(self.default_changes())
        if key == "Control":
            rate1, t_change, rate2 = self.infusion_control
        else:
            rate1, t_change, rate2 = self.infusion_cirrhosis
            changes["f_tissue_loss"] = Q(self.cirrhosis_map[key], "dimensionless")
            changes["f_shunts"] = Q(self.cirrhosis_map[key], "dimensionless")
        changes["Ri_icg"] = Q(rate1, "mg/min")
        return Simulation(
            time_unit="min",
            end=end,
            preinit_changes=changes,
            changes=[Change(t_change, {"Ri_icg": Q(rate2, "mg/min")})],
        )

    def _sid(self, key: str) -> str:
        """Get the part of the simulation keys of a key of `cirrhosis_map`."""
        return "control" if key == "Control" else self.to_sid(key)

    @override
    def simulations(self) -> dict[str, Simulation | Scan]:
        sims: dict[str, Simulation | Scan] = {}
        for key in self.cirrhosis_map:
            sims[f"icg_inf_{self._sid(key)}"] = self._infusion(key, end=76)

        # the scans of the hepatic blood flow, the steady state at the end
        dimension_bloodflow = Dimension(
            "dim_bloodflow", values=self.bloodflow_scan_changes(stop=3.5)
        )
        for key in self.cirrhosis_map:
            sims[f"scan_bloodflow_{self._sid(key)}"] = Scan(
                simulation=self._infusion(key, end=106),
                dimensions=[dimension_bloodflow],
            )
        return sims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        # chronic hepatitis B or C without fibrosis, sedated with pethidine
        metadata = PKPDMappingMetaData(
            tissue=Tissue.PLASMA,
            route=Route.IV,
            application_form=ApplicationForm.SOLUTION,
            dosing=Dosing.CONSTANT_INFUSION,
            health=Health.HEPATITIS,
            fasting=Fasting.FASTED,
            coadministration=Coadministration.PETHIDINE,
            data=PKPDData.PK,
            data_type=DataType.MEAN,
        )
        return {
            f"fm_{key}_hepatitis": FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"Tab3_{key}_Controls",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task="task_icg_inf_control", xid="time", yid=yid
                ),
                metadata=metadata,
            )
            for key, yid in [("er", "ER_icg"), ("cl", "CLinfusion_icg")]
        }

    @override
    def figures(self) -> dict[str, Figure]:
        return {"fig_icg": self._figure_icg()}

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        return {"fig_1": self._figure_bloodflow_vs_clearance()}

    def _figure_bloodflow_vs_clearance(self) -> FigureMPL:
        """Figure of the clearance against the hepatic blood flow (Fig 1).

        The model curves are the steady state of the infusion of the study at
        the end of the scans; the dashed line is the hepatic blood flow of the
        model at 75 kg.
        """
        f, ax = base_figure()
        ax.axvline(
            x=self.hepatic_bloodflow().to(self.unit_bloodflow).magnitude,
            color="grey",
            linestyle="--",
        )

        for key in self.cirrhosis_map:
            xres = self.results[f"task_scan_bloodflow_{self._sid(key)}"]
            # the blood flow is constant, the value at the last time point
            x = self.scan_values(xres, "Qh", self.unit_bloodflow)
            self.plot_band(
                ax,
                x,
                *self.scan_at_end(xres, "CLinfusion_icg", self.unit_clearance),
                color=self.cirrhosis_colors[key],
                label=key,
            )
        ax.set_xlabel(f"Hepatic blood flow [{self.unit_bloodflow}]", fontdict=self.font)
        ax.set_ylabel(
            f"{self.pk_labels['cl']} [{self.unit_clearance}]", fontdict=self.font
        )

        # data
        for group in self.groups:
            dset = self._datasets[self.to_sid(f"Fig1_hbfvsCl_{group}")]
            ax.plot(
                dset.get_quantity("x").magnitude,
                dset.get_quantity("y").magnitude,
                marker="s",
                linestyle="",
                color=self.group_colors[group],
                label=f"Gadano 1997 {group.lower()}",
                ms=3,
                alpha=0.7,
            )

        ax.set_xlim(left=0, right=4)
        ax.set_ylim(bottom=0, top=1050)
        ax.legend(fontsize=self.legend_font_size - 3, loc="upper right")
        return f

    def _figure_icg(self) -> Figure:
        """Figure of the extraction ratio of indocyanine green (Table 3)."""
        fig = Figure(experiment=self, sid="Tab3", width=6, height=6)
        plots = fig.create_plots(xaxis=Axis("Time", unit=self.unit_time), legend=True)
        plots[0].set_yaxis("ICG extraction ratio", self.unit_er, min=0, max=1.1)

        # simulation
        for key in self.cirrhosis_map:
            plots[0].add_data(
                task=f"task_icg_inf_{self._sid(key)}",
                xid="time",
                yid="ER_icg",
                label=key,
                color=self.cirrhosis_colors[key],
            )

        # data
        for group in self.groups:
            plots[0].add_data(
                dataset=self.to_sid(f"Tab3_er_{group}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"Gadano 1997 {group.lower()}",
                marker="s",
                color=self.group_colors[group],
            )
        return fig


if __name__ == "__main__":
    MODEL.run_experiments(Gadano1997, output_dir=Gadano1997.__name__)
