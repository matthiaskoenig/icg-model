"""Simulation experiment of the study Kamimori2000."""

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
from sbmlsim.simulation import Simulation

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment


class Kamimori2000(ICGSimulationExperiment):
    """Simulation experiment of Kamimori2000.

    Timecourses of indocyanine green in healthy women (58.6 kg) in the
    follicular, the ovulatory and the luteal phase of the menstrual cycle after
    an intravenous bolus of 0.5 mg/kg, a randomized crossover of the same 13
    women (11 completed), fasted, with triazolam 0.25 mg orally in every phase.
    The means of Fig 2 are the means of the 7 individual profiles shown per
    phase (6 at 5 min in the follicular phase); the means at 20 and 30 min are
    means of the 1 - 4 subjects with a detectable concentration and are not
    part of the data. The study found no effect of the phase.
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        fig_id = "Fig2A"
        for substance, df in self.load_dataframes_by_substance(fig_id).items():
            dset = DataSet.from_df(df, self.ureg)
            if substance == "indocyanine green":
                dset.unit_conversion("mean", 1 / self.Mr.icg)
            for group in df.group.unique():
                if "Follicular" in group:
                    group_type = "Follicular"
                elif "Ovulatory" in group:
                    group_type = "Ovulatory"
                else:
                    group_type = "Luteal"
                dsets[self.to_sid(f"{fig_id}_{substance} {group_type}")] = dset[
                    dset.group == group
                ]
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        return self._simulation_icg()

    def _simulation_icg(self) -> dict[str, Simulation]:
        """Intravenous bolus of indocyanine green, 0.5 mg/kg."""
        # dose of 0.5 mg/kg * average bodyweight of the group
        icg_iv = 0.5 * 58.6
        return {
            "icg_iv": Simulation(
                time_unit="min",
                end=30,
                preinit_changes={
                    **self.default_changes(),
                    "IVDOSE_icg": Q(icg_iv, "mg"),
                    "BW": Q(58.6, "kg"),
                },
            )
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        return {
            self.to_sid(f"fm_Cve_icg_{key}"): FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(f"Fig2A_indocyanine green {phase}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    # the smallest number of subjects of the time points
                    count=int(
                        self._datasets[self.to_sid(f"Fig2A_indocyanine green {phase}")][
                            "count"
                        ].min()
                    ),
                ),
                observable=FitData(
                    self, task="task_icg_iv", xid="time", yid="[Cve_icg]"
                ),
                # triazolam 0.25 mg orally has no member in the vocabulary and
                # no known effect on the hepatic uptake of indocyanine green
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.IV,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE,
                    data=PKPDData.PK,
                    data_type=DataType.MEAN,
                ),
            )
            for key, phase in [
                ("Fol", "Follicular"),
                ("Ovu", "Ovulatory"),
                ("Lut", "Luteal"),
            ]
        }

    @override
    def figures(self) -> dict[str, Figure]:
        return {"fig1": self._figure_icg()}

    def _figure_icg(self) -> Figure:
        """Figure of indocyanine green in the venous plasma."""
        fig = Figure(experiment=self, sid="Fig2A", width=6, height=6)
        plots = fig.create_plots(xaxis=Axis("Time", unit=self.unit_time), legend=True)
        plots[0].set_yaxis("ICG", self.unit_icg)

        # simulation
        plots[0].add_data(
            task="task_icg_iv",
            xid="time",
            yid="[Cve_icg]",
            xunit=self.unit_time,
            yunit=self.unit_icg,
            label="Sim: ICG (venous plasma)",
            color="black",
        )

        # data
        for phase, color in [
            ("Follicular", "blue"),
            ("Ovulatory", "red"),
            ("Luteal", "green"),
        ]:
            plots[0].add_data(
                dataset=self.to_sid(f"Fig2A_indocyanine green {phase}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                color=color,
                label=f"Kamimori 2000 {phase.lower()}",
            )
        return fig


if __name__ == "__main__":
    MODEL.run_experiments(Kamimori2000, output_dir=Kamimori2000.__name__)
