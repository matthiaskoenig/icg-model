"""Simulation experiment of the study Leevy1962."""

from typing import ClassVar, override

import numpy as np
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
from pkpdsim.scan import end_values
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan, Simulation
from sbmlsim.units import Quantity

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class Leevy1962(ICGSimulationExperiment):
    """Simulation experiment of Leevy1962.

    Extraction ratio and concentrations of indocyanine green during a constant
    infusion in healthy subjects and in cirrhosis (hepatic vein catheter, Fick
    principle), the dependency of the extraction ratio on the dose, on the
    hepatic blood flow and on the plasma bilirubin (with the data of
    Caesar1961).

    - Table I: 14 healthy men (BSA 1.46 - 1.95 m^2, mean 1.784), 0.3
      mg/m^2/min after a loading dose of 10 mg, arterial and hepatic venous
      samples every 5 min at steady state (45 - 60 min); the "peripheral
      conc." is arterial plasma in mg%. The removal rate of Table I (0.99
      mg/min) is 1.8 times the stated infusion (0.54 mg/min), so that the
      concentrations cannot be described with the stated infusion: the
      extraction ratio (independent of the rate) is the robust result.
    - Table IV: 0.3 mg/m^2/min is replaced by 1.5 mg/m^2/min after 10 mg in 4
      controls (BSA 1.91 m^2) and 3 cirrhosis patients (1.75 m^2), with
      bromsulphalein (BSP) given simultaneously, a competitor of the hepatic
      uptake which is not modelled.
    - Table IX: 2 healthy subjects, boluses of 0.15, 0.25 and 0.5 mg/kg,
      samples every 2 min for 14 min; the time of the extraction ratio is
      not reported (8 min, the middle of the samples). The body weight is
      not reported, the model's 75 kg are used.
    - Fig 1: one healthy subject, 0.3 mg/m^2/min; the concentrations of the
      figure are 10 times below the axis of the curation and 45 times below
      Table I: the scale is unknown.
    """

    #: mean body surface area of the subjects of Table I [m^2]
    bsa: ClassVar[float] = 1.784
    #: mean body surface area of the subjects of Table IV [m^2], the controls
    #: and the cirrhosis patients
    bsa_tab4: ClassVar[dict[str, float]] = {"Control": 1.91, "cirrhosis": 1.75}
    #: the boluses of Table IX [mg/kg]
    doses: ClassVar[list[float]] = [0.15, 0.25, 0.5]
    #: the body weight of the boluses [kg]
    bodyweight: ClassVar[float] = 75.0

    @override
    def datasets(self) -> dict[str, DataSet]:
        # a dataset has a single unit: the tables are split by measurement
        dsets = {}

        df = self.load_dataframe("Tab1Mean")
        for measurement in df.measurement.unique():
            dset = DataSet.from_df(df[(df.measurement == measurement)], self.ureg)
            if measurement == "concentration_ss":
                # [mg/dl] = [mg%]
                dset.unit_conversion("mean", 1 / self.Mr.icg)
            elif measurement == "extraction_ratio":
                dset.unit_conversion("mean", 1 / Q(100, "percent"))
            dsets[self.to_sid(f"Tab1Mean_{measurement}")] = dset

        # Table I and II by measurement and substance
        data_map = {
            ("concentration_ss", "albumin"): "alb",
            ("concentration_ss", "bilirubin"): "bil",
            ("concentration_ss", "indocyanine green"): "icg",
            ("extraction_ratio", "indocyanine green"): "ER",
            ("hbf", "indocyanine green"): "HBF",
        }
        for fig_id in ["Tab1", "Tab2"]:
            df = self.load_dataframe(fig_id)
            for (measurement, substance), key in data_map.items():
                df_subset = df[
                    (df.measurement == measurement) & (df.substance == substance)
                ]
                if df_subset.empty:
                    continue
                dset = DataSet.from_df(df_subset, self.ureg)
                if key == "icg":
                    dset.unit_conversion("value", 1 / self.Mr.icg)
                dsets[self.to_sid(f"{fig_id}_{key}")] = dset

        # the arterial and the hepatic venous concentrations of Fig 1
        df = self.load_dataframe("Fig1")
        for label in ["ar", "hv"]:
            dset = DataSet.from_df(df[df.label == f"Fig1_{label}"], self.ureg)
            dset.unit_conversion("value", 1 / self.Mr.icg)
            dsets[self.to_sid(f"Fig1_icg_{label}")] = dset

        df = self.load_dataframe("Tab9")
        for intervention in df.intervention.unique():
            dset = DataSet.from_df(df[df.intervention == intervention], self.ureg)
            dsets[self.to_sid(f"Tab9_ER_{intervention}")] = dset

        df = self.load_dataframe("Tab4")
        for measurement in df.measurement.unique():
            for group in df.group.unique():
                dset = DataSet.from_df(
                    df[(df.measurement == measurement) & (df.group == group)],
                    self.ureg,
                )
                if measurement == "concentration_ss":
                    # [mg/dl] = [mg%]
                    dset.unit_conversion("value", 1 / self.Mr.icg)
                elif measurement == "extraction_ratio":
                    dset.unit_conversion("value", 1 / Q(100, "percent"))
                dsets[self.to_sid(f"Tab4_{measurement}_{group}")] = dset

        # extraction ratio and bilirubin of Caesar1961
        df_caesar = self.load_dataframe("Tab3", study="Caesar1961")
        caesar_groups = ["Cirrhosis", "Controls"]
        for measurement in df_caesar.measurement.unique():
            for group in caesar_groups:
                df_subset = df_caesar[
                    (df_caesar.measurement == measurement)
                    & (df_caesar.group == group)
                    & (df_caesar.substance == "indocyanine green")
                ]
                dset = DataSet.from_df(df_subset, self.ureg)
                dsets[self.to_sid(f"Caesar1961_Tab3_ICG_{measurement}_{group}")] = dset
        for group in caesar_groups:
            df_subset = df_caesar[
                (df_caesar.group == group) & (df_caesar.substance == "bilirubin")
            ]
            dset = DataSet.from_df(df_subset, self.ureg)
            dsets[self.to_sid(f"Caesar1961_Tab3_bilirubin_{group}")] = dset

        return dsets

    def _infusion(self, rate: float, bsa: float, **changes: Quantity) -> Simulation:
        """Get the simulation of an infusion after a loading dose of 10 mg.

        Args:
            rate: the infusion rate [mg/min/m^2].
            bsa: the body surface area [m^2].
            changes: further changes of the simulation.
        """
        return Simulation(
            time_unit="min",
            end=150,
            preinit_changes={
                **self.default_changes(),
                "IVDOSE_icg": Q(10, "mg"),
                "Ri_icg": Q(rate, "mg/min/m^2") * Q(bsa, "m^2"),
                **changes,
            },
        )

    @override
    def simulations(self) -> dict[str, Simulation | Scan]:
        tcsims: dict[str, Simulation | Scan] = {}
        # 0.3 mg/min/m^2 * 1.784 m^2 = 0.54 mg/min after 10 mg
        tcsims["icg_inf03"] = self._infusion(rate=0.3, bsa=self.bsa)

        for dose in self.doses:
            tcsims[self.to_sid(f"icg_iv{dose}")] = Simulation(
                time_unit="min",
                end=30,
                preinit_changes={
                    **self.default_changes(),
                    "IVDOSE_icg": Q(dose * self.bodyweight, "mg"),
                },
            )

        for key, value in self.cirrhosis_map.items():
            cirrhosis_changes = {
                "f_tissue_loss": Q(value, "dimensionless"),
                "f_shunts": Q(value, "dimensionless"),
            }
            bsa = self.bsa_tab4["Control" if key == "Control" else "cirrhosis"]
            tcsims[self.to_sid(f"icg_inf15_{key}")] = self._infusion(
                rate=1.5, bsa=bsa, **cirrhosis_changes
            )
            tcsims[self.to_sid(f"scan_bilirubin_{key}")] = Scan(
                simulation=self._infusion(rate=0.3, bsa=self.bsa, **cirrhosis_changes),
                dimensions=[
                    Dimension("dim_bil", values=self.bilirubin_scan_changes()),
                ],
            )
            tcsims[self.to_sid(f"scan_bloodflow_{key}")] = Scan(
                simulation=self._infusion(rate=0.3, bsa=self.bsa, **cirrhosis_changes),
                dimensions=[
                    Dimension("dim_bloodflow", values=self.bloodflow_scan_changes())
                ],
            )
        return tcsims

    @staticmethod
    def _metadata(dosing: Dosing, data_type: DataType) -> PKPDMappingMetaData:
        """Get the metadata of a fit mapping of the healthy subjects."""
        return PKPDMappingMetaData(
            tissue=Tissue.PLASMA,
            route=Route.IV,
            application_form=ApplicationForm.SOLUTION,
            dosing=dosing,
            health=Health.HEALTHY,
            fasting=Fasting.NR,
            data=PKPDData.PK,
            data_type=data_type,
        )

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        infusion_mean = self._metadata(Dosing.CONSTANT_INFUSION, DataType.MEAN)
        infusion_individual = self._metadata(
            Dosing.CONSTANT_INFUSION, DataType.INDIVIDUAL
        )
        mappings = {
            # Table I, the extraction ratio of the 14 subjects
            "fm_er_icg": FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset="Tab1Mean_extraction_ratio",
                    xid="time",
                    yid="mean",
                    count="count",
                    yid_sd="mean_sd",
                ),
                observable=FitData(
                    self, task="task_icg_inf03", xid="time", yid="ER_icg"
                ),
                metadata=infusion_mean,
            ),
            # Table I, the arterial concentration, inconsistent with the
            # stated infusion rate
            "fm_Car_icg_tab1": FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset="Tab1Mean_concentration_ss",
                    xid="time",
                    yid="mean",
                    count="count",
                    yid_sd="mean_sd",
                ),
                observable=FitData(
                    self, task="task_icg_inf03", xid="time", yid="[Car_icg]"
                ),
                metadata=infusion_mean,
            ),
            # Fig 1, the timecourse of one subject with an unknown scale
            "fm_Chv_icg_fig1": FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset="Fig1_icg_hv",
                    xid="time",
                    yid="value",
                    count="count",
                ),
                observable=FitData(
                    self, task="task_icg_inf03", xid="time", yid="[Chv_icg]"
                ),
                metadata=infusion_individual,
            ),
            "fm_Car_icg_fig1": FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset="Fig1_icg_ar",
                    xid="time",
                    yid="value",
                    count="count",
                ),
                observable=FitData(
                    self, task="task_icg_inf03", xid="time", yid="[Car_icg]"
                ),
                metadata=infusion_individual,
            ),
        }
        # Table IX, the dose dependency of the extraction ratio (2 subjects)
        for dose in self.doses:
            mappings[self.to_sid(f"fm_er_{dose}_icg")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"Tab9_ER_ICG{dose}"),
                    xid="time",
                    yid="value",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=self.to_sid(f"task_icg_iv{dose}"),
                    xid="time",
                    yid="ER_icg",
                ),
                metadata=self._metadata(Dosing.SINGLE, DataType.INDIVIDUAL),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self._figure_dose(),
            **self._figures_icg03(),
            **self._figures_icg15(),
        }

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        return {
            **self._figure_er_vs_hbf(),
            **self._figure_er_vs_bil(),
        }

    def _figure_dose(self) -> dict[str, Figure]:
        """Extraction ratio after the intravenous doses."""
        fig = Figure(experiment=self, sid="Tab9", width=6, height=6)
        plots = fig.create_plots(xaxis=Axis("Time", unit=self.unit_time), legend=True)
        plots[0].set_yaxis("ICG extraction ratio", self.unit_er)

        colors = {0.15: "orange", 0.25: "green", 0.5: "brown"}
        for dose, color in colors.items():
            # simulation
            plots[0].add_data(
                task=self.to_sid(f"task_icg_iv{dose}"),
                xid="time",
                yid="ER_icg",
                label=f"Sim: {dose} mg/kg",
                color=color,
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(f"Tab9_ER_ICG{dose}"),
                xid="time",
                yid="value",
                color=color,
                label=f"Leevy 1962 {dose} mg/kg",
                linestyle="",
            )
        return {"fig_dose": fig}

    def _figures_icg03(self) -> dict[str, Figure]:
        """Concentrations and extraction ratio during the infusion of 0.3 mg/min/m^2."""
        figures = {}
        for plot in ["car", "er"]:
            fig = Figure(experiment=self, sid="Tab1", width=6, height=6)
            plots = fig.create_plots(
                xaxis=Axis("Time", unit=self.unit_time), legend=True
            )
            if plot == "car":
                plots[0].set_yaxis("ICG", self.unit_icg)
                # simulation
                plots[0].add_data(
                    task="task_icg_inf03",
                    xid="time",
                    yid="[Car_icg]",
                    label="Sim: ICG (arterial plasma)",
                    color="red",
                    linestyle="-",
                )
                plots[0].add_data(
                    task="task_icg_inf03",
                    xid="time",
                    yid="[Chv_icg]",
                    label="Sim: ICG (hepatic venous plasma)",
                    color="blue",
                )

                # data
                plots[0].add_data(
                    dataset="Tab1Mean_concentration_ss",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    color="red",
                    label="Leevy 1962 ar (Table I)",
                )
                plots[0].add_data(
                    dataset="Fig1_icg_ar",
                    xid="time",
                    yid="value",
                    count="count",
                    color="red",
                    marker="o",
                    label="Leevy 1962 ar (Fig 1)",
                )
                plots[0].add_data(
                    dataset="Fig1_icg_hv",
                    xid="time",
                    yid="value",
                    count="count",
                    color="blue",
                    marker="o",
                    label="Leevy 1962 hv (Fig 1)",
                )

            else:
                plots[0].set_yaxis("ICG extraction ratio", self.unit_er)

                # simulation
                plots[0].add_data(
                    task="task_icg_inf03",
                    xid="time",
                    yid="ER_icg",
                    label="Sim: ICG extraction ratio",
                    color="black",
                )

                # data
                plots[0].add_data(
                    dataset="Tab1Mean_extraction_ratio",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    color="black",
                    label="Leevy 1962",
                )
                plots[0].add_data(
                    dataset="Tab1_ER",
                    xid="time",
                    yid="value",
                    color="black",
                    label="Leevy 1962",
                    linestyle="",
                )
            figures[self.to_sid(f"fig_icg03{plot}")] = fig
        return figures

    def _figures_icg15(self) -> dict[str, Figure]:
        """Concentrations and extraction ratio during the infusion of 1.5 mg/min/m^2."""
        figures = {}
        groups = {"Controls": "black", "Cirrhosis": "purple"}
        for plot in ["ar", "er"]:
            fig = Figure(experiment=self, sid="Tab4", width=6, height=6)
            plots = fig.create_plots(
                xaxis=Axis("Time", unit=self.unit_time), legend=True
            )
            if plot == "ar":
                plots[0].set_yaxis("ICG (arterial plasma)", self.unit_icg)
                # simulation
                for key in self.cirrhosis_map:
                    plots[0].add_data(
                        task=self.to_sid(f"task_icg_inf15_{key}"),
                        xid="time",
                        yid="[Car_icg]",
                        label=key,
                        color=self.cirrhosis_colors[key],
                    )
                # data
                for group, color in groups.items():
                    plots[0].add_data(
                        dataset=self.to_sid(f"Tab4_concentration_ss_{group}"),
                        xid="time",
                        yid="value",
                        label=f"Leevy 1962 ar {group.lower()}",
                        marker="s",
                        color=color,
                        linestyle="",
                    )
            else:
                plots[0].set_yaxis("ICG extraction ratio", self.unit_er)
                # simulation
                for key in self.cirrhosis_map:
                    plots[0].add_data(
                        task=self.to_sid(f"task_icg_inf15_{key}"),
                        xid="time",
                        yid="ER_icg",
                        label=key,
                        color=self.cirrhosis_colors[key],
                    )
                # data
                for group, color in groups.items():
                    plots[0].add_data(
                        dataset=self.to_sid(f"Tab4_extraction_ratio_{group}"),
                        xid="time",
                        yid="value",
                        label=f"Leevy 1962 {group.lower()}",
                        marker="s",
                        color=color,
                        linestyle="",
                    )
            figures[self.to_sid(f"fig_icg15{plot}")] = fig
        return figures

    def _er_end(self, sim_key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get the extraction ratio at the end of the simulations of a scan.

        Returns:
            The values at the points of the scan (the median of the draws in a
            run over them), the lower and the upper quantile of the draws.
        """
        er = end_values(self.results[self.to_sid(f"task_{sim_key}")], "ER_icg")
        if self.uncertainty is None:
            values = np.asarray(er.values, dtype=float)
            return values, values, values
        across = self.uncertainty.id
        low, high = self.uncertainty.quantiles
        return (
            np.asarray(er.median(dim=across).values, dtype=float),
            np.asarray(er.quantile(low, dim=across).values, dtype=float),
            np.asarray(er.quantile(high, dim=across).values, dtype=float),
        )

    def _plot_er(self, ax: Axes, x: np.ndarray, sim_key: str, key: str) -> None:
        """Plot the extraction ratio of a scan of a stage of cirrhosis."""
        y, lower, upper = self._er_end(sim_key)
        color = self.cirrhosis_colors[key]
        if self.uncertainty is not None:
            ax.fill_between(x, lower, upper, color=color, alpha=0.2, linewidth=0)
        ax.plot(x, y, marker="o", linestyle="-", color=color, label=key)

    def _figure_er_vs_hbf(self) -> dict[str, FigureMPL]:
        """Extraction ratio of indocyanine green vs hepatic blood flow."""
        f, ax = base_figure()
        ax.set_xlabel("Hepatic blood flow/BSA\n[l/min/m$^2$]", fontdict=self.font)
        ax.set_ylabel("ICG extraction ratio [-]", fontdict=self.font)
        # the hepatic blood flow of the model
        bsa = Q(self.bsa, "m^2")
        ax.axvline(
            x=(self.hepatic_bloodflow() / bsa).to("l/min/m^2").magnitude,
            color="grey",
            linestyle="--",
        )
        # the hepatic blood flow of the points of the scan
        f_bloodflow = self.bloodflow_scan_changes()["f_bloodflow"].magnitude
        x = np.array(
            [
                (self.hepatic_bloodflow(f_bloodflow=value) / bsa)
                .to("l/min/m^2")
                .magnitude
                for value in f_bloodflow
            ]
        )
        for key in self.cirrhosis_map:
            self._plot_er(ax, x, f"scan_bloodflow_{key}", key)

        # data: percent -> dimensionless, ml -> l
        for fig_id, color, label in [
            ("Tab2", "purple", "Leevy 1962 cirrhosis"),
            ("Tab1", "black", "Leevy 1962 healthy"),
        ]:
            ydset = self._datasets[self.to_sid(f"{fig_id}_ER")]
            xdset = self._datasets[self.to_sid(f"{fig_id}_HBF")]
            ax.plot(
                xdset.get_quantity("value").to("l/min/m^2").magnitude,
                ydset.get_quantity("value").to("dimensionless").magnitude,
                marker="s",
                linestyle="",
                color=color,
                label=label,
            )
        # the legend above the data
        ax.set_ylim(bottom=0, top=1.4)
        ax.legend(fontsize=self.legend_font_size - 2, loc="upper right", ncols=2)

        return {"fig_hbf": f}

    def _figure_er_vs_bil(self) -> dict[str, FigureMPL]:
        """Extraction ratio of indocyanine green vs plasma bilirubin."""
        f, ax = base_figure()
        ax.set_xlabel("Plasma bilirubin [mmole/l]", fontdict=self.font)
        ax.set_ylabel("ICG extraction ratio [-]", fontdict=self.font)
        # the bilirubin of the model
        ax.axvline(x=0.01, color="grey", linestyle="--")

        x = (
            self.bilirubin_scan_changes()["[LI__bil_ext]"]
            .to(self.unit_bilirubin)
            .magnitude
        )
        for key in self.cirrhosis_map:
            self._plot_er(ax, x, f"scan_bilirubin_{key}", key)

        # data: datasets of bilirubin (x) and extraction ratio (y), marker, color
        plot_map = {
            "Leevy 1962 cirrhosis": ("Tab2_bil", "Tab2_ER", "s", "purple"),
            "Caesar 1961 healthy": (
                "Caesar1961_Tab3_bilirubin_Controls",
                "Caesar1961_Tab3_ICG_extraction_ratio_Controls",
                "D",
                "black",
            ),
            "Caesar 1961 cirrhosis": (
                "Caesar1961_Tab3_bilirubin_Cirrhosis",
                "Caesar1961_Tab3_ICG_extraction_ratio_Cirrhosis",
                "D",
                "purple",
            ),
        }
        for label, (xdset_key, ydset_key, marker, color) in plot_map.items():
            xdset = self._datasets[self.to_sid(xdset_key)]
            ydset = self._datasets[self.to_sid(ydset_key)]
            bilirubin = xdset.get_quantity("value")
            x_data = self.bilirubin_mM(bilirubin.magnitude, str(bilirubin.units)).to(
                self.unit_bilirubin
            )
            ax.plot(
                x_data.magnitude,
                ydset.get_quantity("value").to("dimensionless").magnitude,
                marker=marker,
                linestyle="",
                color=color,
                label=label,
            )

        ax.set_xscale("log")
        # the legend above the data
        ax.set_ylim(bottom=0, top=1.4)
        ax.legend(fontsize=self.legend_font_size - 2, loc="upper right", ncols=2)

        return {"fig_bil": f}
