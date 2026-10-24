"""Simulation experiment of the study Niemann2000."""

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
from sbmlsim.simulation import Simulation
from sbmlsim.units import Quantity

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment


class Niemann2000(ICGSimulationExperiment):
    """Simulation experiment of Niemann2000.

    Arterial timecourse of one of four healthy, fasted men (80 kg) after an
    intravenous bolus of 10 mg indocyanine green at -0.05 min (with antipyrine
    and fentanyl as markers), without and during a target controlled infusion
    of propranolol (Fig 2B). The data points before 2 min (the first pass and
    the recirculation peaks, which the model has no transit times for) are not
    used.

    Propranolol lowered the cardiac output from 10.6 to 4.1 l/min (0.39) and
    the elimination (blood) clearance of indocyanine green from 1.20 to 0.83
    l/min (0.69, Table II). The beta blockade reduces the splanchnic and the
    hepatic blood flow, which limits the clearance of indocyanine green. In the
    model the hepatic blood flow is a fixed fraction of the cardiac output, so
    the cardiac output of the data (0.39) would reduce the hepatic blood flow
    by 61 % (model: kel 0.092 instead of 0.195 /min, data 0.125 and 0.158 /min,
    a ratio of 0.79); the propranolol arm is simulated with the hepatic blood
    flow reduced by the measured ratio of the clearance (`f_bloodflow` 0.69,
    model kel ratio 0.75) and the cardiac output unchanged, which acts only on
    the mixing in the first minutes.
    """

    #: hematocrit of the conversion of the blood concentrations of the paper to
    #: plasma: the authors multiplied the plasma concentrations by 1 - Hct of
    #: every subject (not reported); the reference hematocrit of young men
    hematocrit: ClassVar[float] = 0.44
    #: factor of the hepatic blood flow during propranolol, the ratio of the
    #: elimination clearances 0.83/1.20 (Table II)
    f_bloodflow_propranolol: ClassVar[float] = 0.69

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        fig_id = "Fig2B"
        for substance, df in self.load_dataframes_by_substance(fig_id).items():
            if substance != "indocyanine green":
                continue
            # the data are blood concentrations, the plasma concentrations of
            # the paper times (1 - Hct); back to plasma, the observable
            df["value"] = df["value"] / (1 - self.hematocrit)
            dset = DataSet.from_df(df, self.ureg)
            dset.unit_conversion("value", 1 / self.Mr.icg)
            for intervention, suffix in [
                ("ICG, APY, FEN", ""),
                ("PROP, ICG, APY, FEN", "_propranolol"),
            ]:
                # only the data points after the mixing (> 2 min)
                dsets[self.to_sid(f"{fig_id}_{substance}{suffix}")] = dset[
                    (dset.intervention == intervention) & (dset.time >= 2.0)
                ]
        return dsets

    def _simulation(self, changes: dict[str, Quantity]) -> Simulation:
        """Get the simulation of the bolus of 10 mg at -0.05 min."""
        return Simulation(
            time_unit="min",
            start=-0.05,
            end=10,
            preinit_changes={
                **self.default_changes(),
                "IVDOSE_icg": Q(10, "mg"),
                "BW": Q(80, "kg"),
                **changes,
            },
        )

    @override
    def simulations(self) -> dict[str, Simulation]:
        return {
            "icg_iv": self._simulation({}),
            "icg_iv_propranolol": self._simulation(
                {"f_bloodflow": Q(self.f_bloodflow_propranolol, "dimensionless")}
            ),
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for suffix, coadministration in [
            ("", Coadministration.NONE),
            # pkpdsim has no vocabulary of propranolol
            ("_propranolol", Coadministration.COMEDICATION),
        ]:
            mappings[f"fm_Car_icg_individual{suffix}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"Fig2B_indocyanine_green{suffix}",
                    xid="time",
                    yid="value",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_icg_iv{suffix}", xid="time", yid="[Car_icg]"
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.IV,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=coadministration,
                    data=PKPDData.PK,
                    data_type=DataType.INDIVIDUAL,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(experiment=self, sid="Fig2B", width=6, height=6)
        plots = fig.create_plots(xaxis=Axis("Time", unit=self.unit_time), legend=True)
        plots[0].set_yaxis("ICG", self.unit_icg, scale="log", min=0.1)

        for suffix, label, color in [
            ("", "control", "black"),
            ("_propranolol", "propranolol", "tab:red"),
        ]:
            plots[0].add_data(
                task=f"task_icg_iv{suffix}",
                xid="time",
                yid="[Car_icg]",
                xunit=self.unit_time,
                yunit=self.unit_icg,
                label=f"Sim: {label} (arterial plasma)",
                color=color,
            )
            plots[0].add_data(
                dataset=f"Fig2B_indocyanine_green{suffix}",
                xid="time",
                yid="value",
                count="count",
                xunit=self.unit_time,
                yunit=self.unit_icg,
                color=color,
                label=f"Niemann 2000 {label}",
            )

        return {"fig2b": fig}
