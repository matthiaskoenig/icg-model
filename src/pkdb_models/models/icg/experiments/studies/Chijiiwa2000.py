"""Simulation experiment of the study Chijiiwa2000."""

from typing import ClassVar, override

import pandas as pd
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
from sbmlsim.simulation import Scan, Simulation
from sbmlsim.units import Quantity

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment


class Chijiiwa2000(ICGSimulationExperiment):
    """Simulation experiment of Chijiiwa2000.

    Biliary excretion of indocyanine green in 14 patients after the relief of a
    complete extrahepatic biliary obstruction by percutaneous transhepatic
    biliary drainage (PTBD), after an intravenous bolus of 0.5 mg/kg in the
    early morning after an overnight fast. The whole bile was collected hourly
    for 5 h through the drainage; Fig 1 gives the excretion rate in % of the
    dose per hour (mean +- SD, n=14), the text the excretion of 5 h, 40 +- 18 %
    of the dose. The serum bilirubin at the test was 1.2 +- 0.5 mg/dl and the
    ICG-R15 11.5 +- 4.7 % (normal < 10 %): the uptake and the canalicular
    transport of these livers recover from the cholestasis, which the model
    does not describe (no impairment is simulated besides the bilirubin), so
    the data are a validation of a limitation of the model. The body weight is
    not reported, the simulation is of 75 kg; the rates are converted with the
    dose of the simulation, so that they are fractions of the dose.
    """

    bodyweight: ClassVar[Quantity] = Q(75, "kg")
    #: the serum bilirubin at the test [mg/dl]
    bilirubin: ClassVar[float] = 1.2
    #: the excretion of 5 h, mean and SD of the text [% of the dose]
    excretion_5h: ClassVar[tuple[float, float]] = (40.0, 18.0)

    @classmethod
    def dose(cls) -> Quantity:
        """Get the dose of the simulation, 0.5 mg/kg."""
        return (Q(0.5, "mg/kg") * cls.bodyweight).to("mg")

    @override
    def datasets(self) -> dict[str, DataSet]:
        fig_id = "Fig1"
        df = self.load_dataframe(fig_id)
        df = df[
            (df.measurement == "excretion rate") & (df.substance == "indocyanine green")
        ]
        # the rate in mg/kg/h is the fraction of the dose of 0.5 mg/kg per hour
        dset = DataSet.from_df(df, self.ureg)
        dset.unit_conversion("mean", self.bodyweight / self.Mr.icg)

        # the cumulative excretion at the end of every hour, the sum of the
        # hourly fractions (no SD), and the excretion of 5 h of the text
        dose = self.dose().magnitude
        fractions = df.meanraw.to_numpy() / 100
        cumulative = pd.DataFrame(
            {
                "time": df.time.to_numpy() + 0.5,
                "time_unit": "hr",
                "mean": fractions.cumsum() * dose,
                "unit": "mg",
                "count": df["count"].to_numpy(),
            }
        )
        dset_cumulative = DataSet.from_df(cumulative, self.ureg)
        dset_cumulative.unit_conversion("mean", 1 / self.Mr.icg)
        mean_5h, sd_5h = self.excretion_5h
        total = pd.DataFrame(
            {
                "time": [5.0],
                "time_unit": "hr",
                "mean": [mean_5h / 100 * dose],
                "sd": [sd_5h / 100 * dose],
                "unit": "mg",
                "count": [14],
            }
        )
        dset_total = DataSet.from_df(total, self.ureg)
        dset_total.unit_conversion("mean", 1 / self.Mr.icg)
        return {
            f"{fig_id}_excretion_rate": dset,
            f"{fig_id}_excretion_cumulative": dset_cumulative,
            "Text_excretion_5h": dset_total,
        }

    @override
    def simulations(self) -> dict[str, Simulation | Scan]:
        return {
            "icg_iv": Simulation(
                time_unit="min",
                end=360,
                preinit_changes={
                    **self.default_changes(),
                    "BW": self.bodyweight,
                    "IVDOSE_icg": self.dose(),
                    "[LI__bil_ext]": self.bilirubin_mM(self.bilirubin, "mg/dl"),
                },
            )
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        return {
            "fm_LI__ICGLI2BI": FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset="Fig1_excretion_rate",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task="task_icg_iv", xid="time", yid="LI__ICGLI2BI"
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.BILE,
                    route=Route.IV,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.BILIARY_OBSTRUCTION,
                    fasting=Fasting.FASTED,
                    data=PKPDData.PK,
                    data_type=DataType.MEAN,
                ),
            )
        }

    @override
    def figures(self) -> dict[str, Figure]:
        return {"fig1": self._figure_icg()}

    def _figure_icg(self) -> Figure:
        """Figure of the biliary excretion rate and the cumulative excretion."""
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=2,
            width=14,
            height=6,
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(xaxis=Axis("Time", unit="hr"), legend=True)
        plots[0].set_yaxis("ICG biliary excretion rate", self.unit_bile_rate)
        plots[1].set_yaxis("ICG biliary excretion", "µmole")

        # simulation
        for plot, yid in zip(plots, ["LI__ICGLI2BI", "Afeces_icg"], strict=True):
            plot.add_data(
                task="task_icg_iv",
                xid="time",
                yid=yid,
                label="Sim: bilirubin 1.2 mg/dl",
                color="black",
            )

        # data
        plots[0].add_data(
            dataset="Fig1_excretion_rate",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            color="black",
            label="Chijiiwa 2000 (hourly)",
        )
        plots[1].add_data(
            dataset="Fig1_excretion_cumulative",
            xid="time",
            yid="mean",
            count="count",
            color="black",
            label="Chijiiwa 2000 (sum of Fig 1)",
        )
        plots[1].add_data(
            dataset="Text_excretion_5h",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            color="grey",
            label="Chijiiwa 2000 (5 h, text)",
        )
        return fig


if __name__ == "__main__":
    MODEL.run_experiments(Chijiiwa2000, output_dir=Chijiiwa2000.__name__)
