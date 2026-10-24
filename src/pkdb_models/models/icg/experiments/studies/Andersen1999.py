"""Simulation experiment of the study Andersen1999."""

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
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Scan, Simulation

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment


class Andersen1999(ICGSimulationExperiment):
    """Simulation experiment of Andersen1999.

    Timecourse of indocyanine green in the venous plasma of a healthy control
    after an intravenous bolus of 0.5 mg/kg, given after an overnight fast
    together with antipyrine (1000 mg) and oxazepam (15 mg) by mouth, a
    cocktail of probe drugs which do not interact with indocyanine green.
    Fig 1 shows one individual of each group; the individual of the patients
    with chronic pancreatitis (58.2 kg, opioids) is no population of the model
    and is not used, nor are the timecourses of antipyrine and oxazepam. The
    body weight is the mean of the controls (83.4 kg), the weight of the
    individual is not reported. The clearance of the controls in Table I,
    771 ml/min (95 % CI 677 - 865), is a literature value of the model.
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        fig_id = "Fig1"
        df = self.load_dataframe(fig_id)
        df = df[df.individual == "ICGControl"]
        dset = DataSet.from_df(df, self.ureg)
        dset.unit_conversion("value", 1 / self.Mr.icg)
        return {f"{fig_id}_ICGControl": dset}

    @override
    def simulations(self) -> dict[str, Simulation | Scan]:
        return self._simulation_icg()

    def _simulation_icg(self) -> dict[str, Simulation | Scan]:
        """Intravenous bolus of indocyanine green, 0.5 mg/kg."""
        bodyweight = Q(83.4, "kg")
        icg_iv = Q(0.5, "mg/kg") * bodyweight
        return {
            "icg_iv": Simulation(
                time_unit="min",
                end=20,
                preinit_changes={
                    **self.default_changes(),
                    "BW": bodyweight,
                    "IVDOSE_icg": icg_iv,
                },
            )
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        return {
            "fm_Cve_icg_individual": FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset="Fig1_ICGControl",
                    xid="time",
                    yid="value",
                    count="count",
                ),
                observable=FitData(
                    self, task="task_icg_iv", xid="time", yid="[Cve_icg]"
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.IV,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.COCKTAIL,
                    data=PKPDData.PK,
                    data_type=DataType.INDIVIDUAL,
                ),
            )
        }

    @override
    def figures(self) -> dict[str, Figure]:
        return {"fig1": self._figure_icg()}

    def _figure_icg(self) -> Figure:
        """Figure of indocyanine green in the venous plasma."""
        fig = Figure(experiment=self, sid="Fig1", width=6, height=6)
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

        # data
        plots[0].add_data(
            dataset="Fig1_ICGControl",
            xid="time",
            yid="value",
            count="count",
            color="black",
            label="Andersen 1999",
        )
        return fig


if __name__ == "__main__":
    MODEL.run_experiments(Andersen1999, output_dir=Andersen1999.__name__)
