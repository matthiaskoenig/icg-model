"""Simulation experiment of the study Grundmann1992."""

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
from sbmlsim.simulation import Dimension, Scan, Simulation

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class Grundmann1992(ICGSimulationExperiment):
    """Simulation experiment of Grundmann1992.

    Timecourses of indocyanine green after two intravenous boli of 0.3 mg/kg,
    45 min apart, in healthy volunteers (Fig 2, the pilot experiment of the
    study, the number of volunteers is not reported), and the clearance and
    the half-life of 40 male patients (ASA I, healthy liver, 18 - 45 y) before
    elective surgery (Table 1), against a scan of the hepatic blood flow.

    The study fits a monoexponential decline to the samples of 2 - 14 min
    (`kel`, `thalf_kel` of the window 2 - 14 min), the clearance is the D/AUC
    of this decline, kel * D/C0 (`cl_kel`), the volume D/C0 (`v0`). The
    hepatic blood flow of Table 1 is CL / (1 - Hct), i.e. an extraction ratio
    of 1 is assumed: it is the blood clearance, plotted at the axis of the
    hepatic blood flow of the model. The body weight is not reported (75 kg).
    The tables of the anesthetics (Tables 2 - 5) are not simulated.
    """

    #: the window of the monoexponential fit of the study [min]
    window: ClassVar[tuple[float, float]] = (2.0, 14.0)

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        fig_id = "Fig2A"
        for substance, df in self.load_dataframes_by_substance(fig_id).items():
            dset = DataSet.from_df(df, self.ureg)
            if substance == "indocyanine green":
                dset.unit_conversion("mean", 1 / self.Mr.icg)
            for intervention in df.intervention.unique():
                intervention_type = "ICG1" if "ICG1" in intervention else "ICG2"
                dsets[self.to_sid(f"{fig_id}_{substance} {intervention_type}")] = dset[
                    dset.intervention == intervention
                ]

        tab_id = "Tab1"
        df = self.load_dataframe(tab_id)
        for measurement_type in ["thalf", "clearance", "hbf"]:
            # a data set holds the rows of one measurement type, each in its
            # unit (sbmlsim reads a column in one unit)
            dset = DataSet.from_df(
                df[df.measurement_type == measurement_type], self.ureg
            )
            dsets[self.to_sid(f"{tab_id}_{measurement_type}")] = dset[
                dset.intervention == "ICG1"
            ]
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation | Scan]:
        return {
            "icg_iv": self.normal_timecourse(dose=0.3 * 75),
            "scan_bloodflow": Scan(
                simulation=self.normal_timecourse(dose=0.3 * 75),
                dimensions=[
                    Dimension("dim_bloodflow", values=self.bloodflow_scan_changes()),
                ],
            ),
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        # the number of volunteers of the pilot experiment is not reported:
        # the mean curve is weighted as one subject
        return {
            self.to_sid(f"fm_Cve_icg{k}"): FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"Fig2A_indocyanine green ICG{k}"),
                    xid="time",
                    yid="mean",
                ),
                observable=FitData(
                    self, task="task_icg_iv", xid="time", yid="[Cve_icg]"
                ),
                weight=1.0,
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.IV,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.NONE,
                    data=PKPDData.PK,
                    data_type=DataType.MEAN,
                ),
            )
            for k in [1, 2]
        }

    @override
    def figures(self) -> dict[str, Figure]:
        return {"fig1": self._figure_icg()}

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk(
            scans=["scan_bloodflow"], window=self.window
        )
        return self._figures_pk_vs_bloodflow()

    def _figure_icg(self) -> Figure:
        """Figure of indocyanine green in the venous plasma (Fig 2)."""
        fig = Figure(experiment=self, sid="Fig2", width=6, height=6)
        plots = fig.create_plots(xaxis=Axis("Time", unit=self.unit_time), legend=True)
        plots[0].set_yaxis("ICG", self.unit_icg)

        # simulation
        plots[0].add_data(
            task="task_icg_iv",
            xid="time",
            yid="[Cve_icg]",
            label="Sim: ICG (venous plasma)",
            color="black",
        )

        # data, the first and the second bolus of the same volunteers
        # (open squares, filled squares as in the figure of the study)
        for k, fill in [(1, "white"), (2, "black")]:
            plots[0].add_data(
                dataset=self.to_sid(f"Fig2A_indocyanine green ICG{k}"),
                xid="time",
                yid="mean",
                color="black",
                marker="s",
                markerfacecolor=fill,
                label=f"Grundmann 1992 bolus {k}",
            )
        return fig

    def _figures_pk_vs_bloodflow(self) -> dict[str, FigureMPL]:
        """Figures of the clearance and the half-life against the blood flow.

        The model is the scan of the hepatic blood flow `Qh`; the data is at the
        hepatic blood flow of the study, CL / (1 - Hct), which assumes an
        extraction ratio of 1 and is the blood clearance (Table 1).
        """
        xres = self.results["task_scan_bloodflow"]
        df = self.pk_dfs["scan_bloodflow"]
        x = self.scan_values(xres, "Qh", self.unit_bloodflow)

        figures: dict[str, FigureMPL] = {}
        for fig_key, pk_key, data_key in [
            ("fig_cl", "cl_kel", "clearance"),
            ("fig_thalf", "thalf_kel", "thalf"),
        ]:
            f, ax = base_figure()
            unit = self.pk_units[pk_key]
            self.plot_pk(
                ax, x, df, pk_key, unit, color="black", label="Blood flow scan"
            )
            ax.set_xlabel(
                f"Hepatic blood flow [{self.unit_bloodflow}]", fontdict=self.font
            )
            ax.set_ylabel(f"{self.pk_labels[pk_key]} [{unit}]", fontdict=self.font)

            # data
            xdset = self._datasets["Tab1_hbf"]
            ydset = self._datasets[self.to_sid(f"Tab1_{data_key}")]
            ax.errorbar(
                xdset.get_quantity("mean").magnitude,
                ydset.get_quantity("mean").to(unit).magnitude,
                yerr=Q(ydset["mean_sd"].values, ydset.uinfo["mean"]).to(unit).magnitude,
                xerr=Q(xdset["mean_sd"].values, xdset.uinfo["mean"]).magnitude,
                marker="s",
                linestyle="",
                color="black",
                label="Grundmann 1992 (n=40)\nat CL/(1 - Hct)",
            )
            ax.set_ylim(bottom=0, top=1200 if pk_key == "cl_kel" else None)
            ax.legend(fontsize=self.legend_font_size)
            figures[fig_key] = f
        return figures


if __name__ == "__main__":
    MODEL.run_experiments(Grundmann1992, output_dir=Grundmann1992.__name__)
