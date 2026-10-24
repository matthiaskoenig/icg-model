"""Simulation experiment of the study Caesar1961."""

from typing import ClassVar, Literal, override

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from pkpdsim.metadata import (
    ApplicationForm,
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

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure, scan_color
from pkdb_models.models.icg.icg_pk import summarize


class Caesar1961(ICGSimulationExperiment):
    """Simulation experiment of Caesar1961.

    Indocyanine green in 10 fasting controls without liver disease (hospital
    patients) and 17 patients with cirrhosis (mostly cryptogenic, several with
    a portacaval anastomosis or ascites).

    Plasma disappearance (Table 1): an intravenous bolus of 0.5 mg/kg, venous
    samples every 3 min for 21 min; the disappearance rate K is the log-linear
    regression, the retention R20 the concentration at 20 min relative to the
    back extrapolation P0, both from the regression of 3 - 21 min of the model
    (`window`, `kel` and `icgr20`). Simulated as a scan of the degree of
    cirrhosis (tissue loss and shunts) and of the bilirubin for every stage of
    cirrhosis.

    Hepatic catheterization (Table 3): after the bolus of 0.5 mg/kg and 21 min
    of samples a second loading dose of about 6 mg and a constant infusion of
    about 0.5 mg/min (controls 0.56, cirrhosis 0.49 mg/min, the means of
    Table 3), 15 min of equilibration and paired femoral arterial (P) and
    hepatic venous (H) samples every 10 min for up to 60 min; the extraction
    ratio is E = (P - H) / P. The means of Table 3 are compared at 60 min. The
    estimated hepatic blood flow of the controls, 1460 ml/min, is a literature
    value of the model's hepatic blood flow. The body weight is not reported,
    the model's 75 kg are used.
    """

    colors: ClassVar[dict[str, str]] = {
        "Controls": "black",
        "Cirrhosis": "purple",
    }
    #: the window of the log-linear regression [min]
    window: ClassVar[tuple[float, float]] = (3.0, 21.0)
    #: the infusion rates of Table 3, the means of the groups [mg/min]
    infusion_rates: ClassVar[dict[str, float]] = {
        "Controls": 0.56,
        "Cirrhosis": 0.49,
    }
    #: the start of the infusion and the second loading dose [min]
    infusion_start: ClassVar[float] = 21.0
    #: the second loading dose [mg]
    loading_dose: ClassVar[float] = 6.0

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets: dict[str, DataSet] = {}

        # Table 1: the disappearance rate, the retention and the bilirubin of
        # every subject in one row
        df = self.load_dataframe("Tab1")
        columns = {
            "kel": "kel",
            "r20": "retention ratio",
            "bilirubin": "concentration_ss",
        }
        frames = [
            df[df.measurement == measurement][["subject", "group", "value", "unit"]]
            .rename(columns={"value": key, "unit": f"{key}_unit"})
            .set_index(["subject", "group"])
            for key, measurement in columns.items()
        ]
        df_subjects = pd.concat(frames, axis=1, join="inner").reset_index()
        for group in self.colors:
            dsets[f"Tab1_{group}"] = DataSet.from_df(
                df_subjects[df_subjects.group == group], self.ureg
            )

        # Table 3: the means of the infusion
        df = self.load_dataframe("Tab3Mean")
        for key, measurement, tissue in [
            ("ER", "extraction_ratio", "plasma"),
            ("Chv", "concentration_ss", "hepatic vein plasma"),
            ("Car", "concentration_ss", "arterial plasma"),
        ]:
            for group in self.colors:
                dset = DataSet.from_df(
                    df[
                        (df.measurement == measurement)
                        & (df.tissue == tissue)
                        & (df.group == group)
                    ],
                    self.ureg,
                )
                if measurement == "concentration_ss":
                    dset.unit_conversion("mean", 1 / self.Mr.icg)
                dsets[f"Tab3_{key}_{group}"] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation | Scan]:
        sims: dict[str, Simulation | Scan] = {
            "scan_cirrhosis_bolus": Scan(
                simulation=self.normal_timecourse(),
                dimensions=[
                    Dimension(
                        "dim_cirrhosis",
                        values=self.cirrhosis_scan_changes(stop=0.95),
                    ),
                ],
            ),
            "four_cirrhosis_scan": Scan(
                simulation=self.normal_timecourse(),
                dimensions=[
                    Dimension(
                        "dim_cirrhosis",
                        values={
                            "f_tissue_loss": Q(
                                list(self.cirrhosis_map.values()), "dimensionless"
                            ),
                            "f_shunts": Q(
                                list(self.cirrhosis_map.values()), "dimensionless"
                            ),
                        },
                    )
                ],
            ),
        }

        tcs = self.cirrhosis_timecourses()
        for key, value in self.cirrhosis_map.items():
            sims[self.to_sid(f"scan_bilirubin_{key}")] = Scan(
                simulation=tcs[key],
                dimensions=[
                    Dimension("dim_bil", values=self.bilirubin_scan_changes()),
                ],
            )
            rate = self.infusion_rates["Controls" if key == "Control" else "Cirrhosis"]
            sims[self.to_sid(f"icg_inf_{key}")] = Simulation(
                time_unit="min",
                end=100,
                preinit_changes={
                    **self.default_changes(),
                    "IVDOSE_icg": Q(37.5, "mg"),
                    "f_tissue_loss": Q(value, "dimensionless"),
                    "f_shunts": Q(value, "dimensionless"),
                },
                changes=[
                    self.add_doses(
                        self.infusion_start,
                        {"IVDOSE_icg": Q(self.loading_dose, "mg")},
                    ),
                    Change(self.infusion_start, {"Ri_icg": Q(rate, "mg/min")}),
                ],
            )
        return sims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for key, yid, dataset in [
            ("fm_er_icg", "ER_icg", "Tab3_ER_Controls"),
            ("fm_Chv_icg", "[Chv_icg]", "Tab3_Chv_Controls"),
            ("fm_Car_icg", "[Car_icg]", "Tab3_Car_Controls"),
        ]:
            mappings[key] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=dataset,
                    xid="time",
                    yid="mean",
                    count="count",
                    yid_sd="mean_sd",
                ),
                observable=FitData(
                    self, task="task_icg_inf_Control", xid="time", yid=yid
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.IV,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.CONSTANT_INFUSION,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    data=PKPDData.PK,
                    data_type=DataType.MEAN,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {"fig_infusion": self._figure_infusion()}

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk(
            scans=[
                "scan_cirrhosis_bolus",
                "four_cirrhosis_scan",
                *(self.to_sid(f"scan_bilirubin_{key}") for key in self.cirrhosis_map),
            ],
            window=self.window,
        )
        return {
            "fig_kel_vs_r20": self._figure_kel_vs_r20(),
            **self._figures_pk_vs_bilirubin(),
        }

    def _median(self, sim_key: str, key: str) -> np.ndarray:
        """Get a parameter of every point of a scan, the median over the draws."""
        df = self.pk_dfs[sim_key]
        if self.uncertainty is None:
            return df[key].to_numpy(dtype=float)
        return summarize(df, key, self.uncertainty.quantiles)["median"].to_numpy()

    def _figure_kel_vs_r20(self) -> FigureMPL:
        """Figure of the retention after 20 min against the elimination rate."""
        x = self._median("scan_cirrhosis_bolus", "kel")
        y = self._median("scan_cirrhosis_bolus", "icgr20")
        x_stars = self._median("four_cirrhosis_scan", "kel")
        y_stars = self._median("four_cirrhosis_scan", "icgr20")

        f, ax = base_figure()
        ax.set_xlabel(
            f"{self.pk_labels['kel']} [{self.pk_units['kel']}]", fontdict=self.font
        )
        ax.set_ylabel(
            f"{self.pk_labels['icgr20']} [{self.pk_units['icgr20']}]",
            fontdict=self.font,
        )
        ax.scatter(
            x,
            y,
            marker="o",
            c=np.linspace(0, 1, len(x)),
            label="Scan: cirrhosis",
            cmap=scan_color(9),
            s=16,
        )
        ax.plot(x, y, linestyle="-", color="grey", zorder=0)
        for k, (key, color) in enumerate(self.cirrhosis_colors.items()):
            ax.plot(
                x_stars[k],
                y_stars[k],
                marker="*",
                linestyle="",
                color=color,
                label=key,
                ms=20,
            )

        # data
        for group, color in self.colors.items():
            dset = self._datasets[f"Tab1_{group}"]
            ax.plot(
                dset.get_quantity("kel").to(self.pk_units["kel"]).magnitude,
                dset.get_quantity("r20").to(self.pk_units["icgr20"]).magnitude,
                marker="s",
                linestyle="",
                color=color,
                label=f"Caesar 1961 {group.lower()}",
            )

        # the legend in the order: scan, stages of cirrhosis, data
        handles, labels = ax.get_legend_handles_labels()
        order = [
            "Scan: cirrhosis",
            *self.cirrhosis_colors,
            *(f"Caesar 1961 {group.lower()}" for group in self.colors),
        ]
        by_label = dict(zip(labels, handles, strict=True))
        ax.legend(
            [by_label[label] for label in order],
            order,
            fontsize=self.legend_font_size,
        )
        ax.set_yscale("log")
        return f

    def _plot_scan(
        self, ax: Axes, sim_key: str, key: str, color: str, label: str
    ) -> None:
        """Plot a parameter of a scan of the bilirubin, with the band of the draws."""
        x = self.bilirubin_scan_changes()["[LI__bil_ext]"].to(self.unit_bilirubin)
        df = self.pk_dfs[sim_key]
        unit_pk = df[f"{key}_unit"].iloc[0]

        def convert(values: np.ndarray) -> np.ndarray:
            """Get the values in the unit of the figure."""
            return Q(values, unit_pk).to(self.pk_units[key]).magnitude

        if self.uncertainty is None:
            y = convert(df[key].to_numpy(dtype=float))
        else:
            summary = summarize(df, key, self.uncertainty.quantiles)
            y, lower, upper = (
                convert(summary[column].to_numpy())
                for column in ["median", "lower", "upper"]
            )
            ax.fill_between(
                x.magnitude, lower, upper, color=color, alpha=0.2, linewidth=0
            )
        ax.plot(x.magnitude, y, marker="o", linestyle="-", color=color, label=label)

    def _figures_pk_vs_bilirubin(self) -> dict[str, FigureMPL]:
        """Figures of the elimination rate and the retention against the bilirubin."""
        figures: dict[str, FigureMPL] = {}
        # the upper limit of the axis leaves room for the legend above the data
        settings: list[tuple[str, str, float, Literal["upper right", "upper left"]]] = [
            ("kel", "kel", 0.5, "upper right"),
            ("icgr20", "r20", 130, "upper left"),
        ]
        for pk_key, data_key, top, loc in settings:
            f, ax = base_figure()
            ax.set_xlabel("Plasma bilirubin [mmole/l]", fontdict=self.font)
            ax.set_ylabel(
                f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                fontdict=self.font,
            )
            # the upper limit of normal bilirubin, 1.2 mg/dl
            ax.axvline(
                x=self.bilirubin_mM(1.2, "mg/dl").magnitude,
                color="grey",
                linestyle="--",
            )
            for key in self.cirrhosis_map:
                self._plot_scan(
                    ax,
                    sim_key=self.to_sid(f"scan_bilirubin_{key}"),
                    key=pk_key,
                    color=self.cirrhosis_colors[key],
                    label=key,
                )

            # data
            for group, color in self.colors.items():
                dset = self._datasets[f"Tab1_{group}"]
                x = dset.get_quantity("bilirubin")
                ax.plot(
                    self.bilirubin_mM(x.magnitude, str(x.units)).magnitude,
                    dset.get_quantity(data_key).to(self.pk_units[pk_key]).magnitude,
                    marker="s",
                    linestyle="",
                    color=color,
                    label=f"Caesar 1961 {group.lower()}",
                )

            ax.set_xscale("log")
            ax.set_ylim(bottom=0, top=top)
            ax.legend(fontsize=self.legend_font_size - 2, loc=loc)
            figures[f"fig_{pk_key}_vs_bilirubin"] = f
        return figures

    def _figure_infusion(self) -> Figure:
        """Figure of the arterial and hepatic venous concentrations and the ER."""
        fig = Figure(
            experiment=self,
            sid="Fig_infusion",
            num_cols=3,
            width=21,
            height=6,
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(xaxis=Axis("Time", unit=self.unit_time), legend=True)
        plots[0].set_yaxis("ICG arterial plasma", self.unit_icg)
        plots[1].set_yaxis("ICG hepatic venous plasma", self.unit_icg)
        plots[2].set_yaxis("ICG extraction ratio", self.unit_er, min=0, max=1)

        for plot, yid, key in zip(
            plots,
            ["[Car_icg]", "[Chv_icg]", "ER_icg"],
            ["Car", "Chv", "ER"],
            strict=True,
        ):
            # simulation
            for stage in self.cirrhosis_map:
                plot.add_data(
                    task=self.to_sid(f"task_icg_inf_{stage}"),
                    xid="time",
                    yid=yid,
                    color=self.cirrhosis_colors[stage],
                    label=stage,
                    yunit=self.unit_er if key == "ER" else self.unit_icg,
                )
            # data
            for group, color in self.colors.items():
                plot.add_data(
                    dataset=f"Tab3_{key}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    color=color,
                    label=f"Caesar 1961 {group.lower()}",
                    yunit=self.unit_er if key == "ER" else self.unit_icg,
                )
        return fig


if __name__ == "__main__":
    MODEL.run_experiments(Caesar1961, output_dir=Caesar1961.__name__)
