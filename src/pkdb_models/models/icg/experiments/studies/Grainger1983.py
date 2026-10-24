"""Simulation experiment of the study Grainger1983."""

from typing import override

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
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment


class Grainger1983(ICGSimulationExperiment):
    """Simulation experiment of Grainger1983.

    Extraction ratio of indocyanine green of individual subjects after an
    intravenous bolus of 0.25 mg/kg (over 10 s) in 5 patients (4 men, 1 woman,
    30 - 70 y) with a normal wedge pressure and a normal liver architecture,
    who underwent a diagnostic transjugular liver biopsy. The direct extraction
    ratio (peripheral - hepatic venous)/peripheral is given at the time of the
    maximal liver content (6 - 10 min, Table 1). The body weight is not
    reported (75 kg).
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        return {
            f"{fig_id}_er": DataSet.from_df(self.load_dataframe(fig_id), self.ureg)
            for fig_id in ["Tab1", "Tab1Mean"]
        }

    @override
    def simulations(self) -> dict[str, Simulation]:
        return {"icg_iv": self.normal_timecourse(dose=0.25 * 75)}

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        return {
            "fm_er": FitMapping(
                self,
                reference=FitData(
                    self, dataset="Tab1_er", xid="time", yid="value", count="count"
                ),
                observable=FitData(self, task="task_icg_iv", xid="time", yid="ER_icg"),
                # the patients were sedated with diazepam i.v., which has no
                # member in the vocabulary and no known effect on the hepatic
                # uptake of indocyanine green
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.IV,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.NONE,
                    data=PKPDData.PK,
                    data_type=DataType.INDIVIDUAL,
                ),
            )
        }

    @override
    def figures(self) -> dict[str, Figure]:
        return {"fig1": self._figure_icg()}

    def _figure_icg(self) -> Figure:
        """Figure of the extraction ratio of indocyanine green."""
        fig = Figure(experiment=self, sid="Fig1", width=6, height=6)
        plots = fig.create_plots(xaxis=Axis("Time", unit=self.unit_time), legend=True)
        plots[0].set_yaxis("ICG extraction ratio", self.unit_er, min=0, max=1.05)

        # simulation
        plots[0].add_data(
            task="task_icg_iv",
            xid="time",
            yid="ER_icg",
            label="Sim: ICG extraction ratio",
            color="black",
        )

        # data, the individuals and their mean
        for dset_key, yid, yid_sd, count, label, marker in [
            ("Tab1_er", "value", None, None, "Grainger 1983 individuals", "o"),
            ("Tab1Mean_er", "mean", "mean_sd", "count", "Grainger 1983", "s"),
        ]:
            plots[0].add_data(
                dataset=dset_key,
                xid="time",
                yid=yid,
                count=count,
                color="black",
                label=label,
                marker=marker,
                linestyle="",
                yid_sd=yid_sd,
            )
        return fig


if __name__ == "__main__":
    MODEL.run_experiments(Grainger1983, output_dir=Grainger1983.__name__)
