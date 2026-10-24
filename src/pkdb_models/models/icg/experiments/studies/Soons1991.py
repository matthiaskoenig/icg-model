"""Simulation experiment of the study Soons1991."""

from typing import ClassVar, override

import numpy as np
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
from pkpdsim.scan import point_timecourse
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class Soons1991(ICGSimulationExperiment):
    """Simulation experiment of Soons1991.

    Plasma concentrations of indocyanine green in six healthy men (72 kg,
    overnight fast, HPLC) during constant infusions of 0.5, 1.0 and 2.0 mg/min
    of 40 min each (part B). The steady state concentration is the mean of the
    last 10 min of a period (Fig 5, mean and SD); the blood clearance per body
    weight is the rate over the steady state plasma concentration and 1 - Ht
    (Table 3). The order of the rates was randomized (Latin square); the
    simulation and the times of the data (35, 75, 115 min) follow subject 4
    (2.0, 0.5, 1.0 mg/min). With a half-life of 3.5 min the plateau of the
    last 10 min of a period does not depend on the rate before it.
    """

    bodyweight: ClassVar[float] = 72.0  # [kg]
    #: the rates of the periods in the order of the simulation [mg/min]
    rates: ClassVar[list[float]] = [2.0, 0.5, 1.0]
    #: the duration of a period [min]
    period: ClassVar[float] = 40.0

    @override
    def datasets(self) -> dict[str, DataSet]:
        df = self.load_dataframe("Fig5mean")
        dset = DataSet.from_df(df, self.ureg)
        dset.unit_conversion("mean", 1 / self.Mr.icg)
        return {
            "Fig5mean_icg_ss": dset,
            "Tab3mean_icg_clb": DataSet.from_df(
                self.load_dataframe("Tab3mean"), self.ureg
            ),
        }

    @override
    def simulations(self) -> dict[str, Simulation]:
        # three infusion rates of 40 min each
        return {
            "icg_inf": Simulation(
                time_unit="min",
                end=3 * self.period,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "Ri_icg": Q(self.rates[0], "mg/min"),
                },
                changes=[
                    Change(k * self.period, {"Ri_icg": Q(rate, "mg/min")})
                    for k, rate in enumerate(self.rates)
                    if k > 0
                ],
            )
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        return {
            "fm_Cve_icg_mean": FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset="Fig5mean_icg_ss",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task="task_icg_inf", xid="time", yid="[Cve_icg]"
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
        }

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(experiment=self, sid="Fig5mean", width=6, height=6)
        plots = fig.create_plots(xaxis=Axis("Time", unit=self.unit_time), legend=True)
        plots[0].set_yaxis("ICG", self.unit_icg)

        # simulation
        plots[0].add_data(
            task="task_icg_inf",
            xid="time",
            yid="[Cve_icg]",
            label="Sim: ICG (venous plasma)",
            color="black",
        )

        # data
        plots[0].add_data(
            dataset="Fig5mean_icg_ss",
            xid="time",
            yid="mean",
            yid_sd="mean_sd",
            count="count",
            color="black",
            linestyle="",
            label="Soons 1991",
        )

        return {"fig1": fig}

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        return {**self._figure_blood_clearance()}

    def _blood_clearances(self, **point: int) -> np.ndarray:
        """Get the blood clearance per body weight of every rate [ml/min/kg].

        The rate over the mean plasma concentration of the last 10 min of its
        period and 1 - HCT, as in the paper.
        """
        res = self.results["task_icg_inf"]
        time, values = point_timecourse(res, "[Cve_icg]", **point)
        t_vec = Q(np.asarray(time), res.units["time"]).to("min").magnitude
        c_vec = np.asarray(values)
        c_unit = res.units["[Cve_icg]"]
        _, hct = point_timecourse(res, "HCT", **point)
        hematocrit = float(np.asarray(hct)[0])
        clearances = []
        for k, rate in enumerate(self.rates):
            end = (k + 1) * self.period
            mask = (t_vec >= end - 10) & (t_vec < end)
            css = (Q(float(np.mean(c_vec[mask])), c_unit) * self.Mr.icg).to("mg/l")
            clb = Q(rate, "mg/min") / css / (1 - hematocrit) / Q(self.bodyweight, "kg")
            clearances.append(clb.to("ml/min/kg").magnitude)
        return np.asarray(clearances)

    def _figure_blood_clearance(self) -> dict[str, FigureMPL]:
        """Blood clearance per body weight vs infusion rate, Table 3."""
        f, ax = base_figure()
        ax.set_xlabel("Infusion rate [mg/min]", fontdict=self.font)
        ax.set_ylabel(f"{self.pk_labels['cl_blood']} [ml/min/kg]", fontdict=self.font)
        order = np.argsort(self.rates)
        rates = np.asarray(self.rates)[order]

        run = self.uncertainty
        if run is None:
            clb = self._blood_clearances()[order]
        else:
            res = self.results["task_icg_inf"]
            draws = np.array(
                [
                    self._blood_clearances(**{run.id: k})[order]
                    for k in range(res.ds.sizes[run.id])
                ]
            )
            clb = np.median(draws, axis=0)
            lower, upper = np.quantile(draws, run.quantiles, axis=0)
            ax.fill_between(rates, lower, upper, color="black", alpha=0.2)
        ax.plot(
            rates,
            clb,
            marker="o",
            linestyle="-",
            color="black",
            label="Sim: ICG blood clearance",
        )

        dset = self._datasets["Tab3mean_icg_clb"]
        data_rates = dset["intervention"].str.removeprefix("ICG_").astype(float)
        ax.errorbar(
            data_rates.to_numpy(),
            dset.get_quantity("mean").to("ml/min/kg").magnitude,
            yerr=dset["mean_sd"].to_numpy(),
            marker="s",
            linestyle="",
            color="black",
            label="Soons 1991 (n=6)",
        )
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=self.legend_font_size, loc="lower right")
        return {"fig_tab3_blood_clearance": f}
