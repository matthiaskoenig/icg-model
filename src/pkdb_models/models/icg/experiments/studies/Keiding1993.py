"""Simulation experiment of the study Keiding1993."""

from typing import ClassVar, override

from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment


class Keiding1993(ICGSimulationExperiment):
    """Simulation experiment of Keiding1993.

    Steady state of a constant infusion of indocyanine green in subject 1
    (Fig 1), a woman of 45 years and 49 kg with autoimmune cirrhosis
    (bilirubin 8 µM, hepatic blood flow 1.78 l/min). The infusion started
    90 min before the first measurement period (the time 0 of Fig 1), the
    concentrations in the artery and the hepatic vein were measured in two
    periods, before and after an infusion of dextran (Macrodex) which diluted
    the plasma proteins; the dilution of the binding proteins is not modelled
    and the second period is not used.

    The subject is encoded from the paper, nothing is fitted:

    - the infusion rate of the subject is not reported (0.07 - 0.16 µmol/min,
      mean 0.10). With Eqn 2 (`K = -Qp ln(1 - E)`), the intrinsic clearance
      K = 0.437 l/min of Table 2 and E = 1 - 0.2391/0.3354 = 0.287 of Fig 1 the
      hepatic plasma flow is Qp = 1.29 l/min, and the removal rate at steady
      state `Qp (ci - co)` = 1.29 * 0.0963 = 0.124 µmol/min is the infusion;
    - the hepatic blood flow 1.78 l/min of Table 2 is the model's blood flow at
      49 kg times `f_bloodflow`; the hematocrit 1 - Qp/Q = 0.27 is not
      reported, it is implied by the plasma flow and the blood flow;
    - the degree of cirrhosis (`f_tissue_loss` = `f_shunts`) is the one at
      which the model has the extraction ratio 0.29 of period 1 at this flow
      (0.3 with the start values of the model, to be derived again after a
      refit).
    """

    #: body weight of subject 1 [kg]
    bodyweight: ClassVar[float] = 49.0
    #: hepatic blood flow of subject 1 [l/min] (Table 2)
    bloodflow: ClassVar[float] = 1.78
    #: hematocrit implied by the plasma flow and the blood flow of subject 1
    hematocrit: ClassVar[float] = 0.27
    #: infusion rate of subject 1 [µmol/min], the removal rate at steady state
    infusion: ClassVar[float] = 0.124
    #: degree of cirrhosis, the extraction ratio of period 1
    cirrhosis: ClassVar[float] = 0.3

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig1")
        # the first period, before the infusion of dextran
        df = df[df.comment != "after macrodex infusion"]
        for label in ["er", "ar", "hv"]:
            dsets[self.to_sid(f"Fig1_ICG_{label}")] = DataSet.from_df(
                df[df.label == f"Fig1_{label}"], self.ureg
            )
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        f_bloodflow = (
            Q(self.bloodflow, "l/min") / self.hepatic_bloodflow(self.bodyweight)
        ).to("dimensionless")
        return {
            "icg_inf": Simulation(
                time_unit="min",
                end=160,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "HCT": Q(self.hematocrit, "dimensionless"),
                    "f_bloodflow": f_bloodflow,
                    "f_shunts": Q(self.cirrhosis, "dimensionless"),
                    "f_tissue_loss": Q(self.cirrhosis, "dimensionless"),
                    "Ri_icg": (Q(self.infusion, "µmole/min") * self.Mr.icg).to(
                        "mg/min"
                    ),
                },
            )
        }

    @override
    def figures(self) -> dict[str, Figure]:
        figures = {}
        for plot in ["hv_ar", "er"]:
            fig = Figure(experiment=self, sid="Fig1", width=6, height=6)
            plots = fig.create_plots(
                xaxis=Axis("Time", unit=self.unit_time, min=0, max=160), legend=True
            )

            if plot == "hv_ar":
                plots[0].set_yaxis("ICG", self.unit_icg)
                # simulation
                plots[0].add_data(
                    task="task_icg_inf",
                    xid="time",
                    yid="[Car_icg]",
                    label="Sim: ICG (arterial plasma)",
                    color="red",
                    linestyle="-",
                )
                plots[0].add_data(
                    task="task_icg_inf",
                    xid="time",
                    yid="[Chv_icg]",
                    label="Sim: ICG (hepatic venous plasma)",
                    color="blue",
                )

                # data
                plots[0].add_data(
                    dataset="Fig1_ICG_ar",
                    xid="time",
                    yid="value",
                    xunit=self.unit_time,
                    yunit=self.unit_icg,
                    count="count",
                    label="Keiding 1993 ar cirrhosis",
                    color="red",
                )
                plots[0].add_data(
                    dataset="Fig1_ICG_hv",
                    xid="time",
                    yid="value",
                    xunit=self.unit_time,
                    yunit=self.unit_icg,
                    count="count",
                    label="Keiding 1993 hv cirrhosis",
                    color="blue",
                )
            else:
                plots[0].set_yaxis("ICG extraction ratio", self.unit_er)
                # simulation
                plots[0].add_data(
                    task="task_icg_inf",
                    xid="time",
                    yid="ER_icg",
                    label="Sim: ICG extraction ratio",
                    color="black",
                )

                # data
                plots[0].add_data(
                    dataset="Fig1_ICG_er",
                    xid="time",
                    yid="value",
                    count="count",
                    label="Keiding 1993 cirrhosis",
                    color="black",
                )

            figures[f"fig_{plot}"] = fig
        return figures
