"""Simulation experiment of the study Martin1976."""

from typing import ClassVar, override

import numpy as np
import pandas as pd
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan
from sbmlsim.units import Quantity

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure
from pkdb_models.models.icg.icg_pk import summarize


class Martin1976(ICGSimulationExperiment):
    """Simulation experiment of Martin1976.

    Dose dependency of the elimination rate of indocyanine green in the 19
    healthy controls (rapid intravenous bolus of 0.5, 2.0, 3.5 and 5.0 mg/kg).
    The paper's rate is the least squares slope of the log concentrations from
    5 to 12 min (Table 2, mean and SD), the window of the regression of the
    model (`kel`). No body weight is reported, the scan of the dose is at the
    75 kg of the model. The groups with Gilbert's syndrome are not simulated:
    they need an uptake of their own per subgroup.
    """

    #: the window of the regression of the elimination rate [min]
    window: ClassVar[tuple[float, float]] = (5.0, 12.0)
    bodyweight: ClassVar[float] = 75.0  # [kg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for df in self.load_dataframes_by_substance("Tab2").values():
            dset = DataSet.from_df(df, self.ureg)
            for group in df.group.unique():
                dsets[f"Tab2_ICGkel_{group}"] = dset[dset.group == group]
        return dsets

    @override
    def simulations(self) -> dict[str, Scan]:
        return {
            "scan_dose": Scan(
                # the scan sets the dose
                simulation=self.normal_timecourse(dose=0),
                dimensions=[
                    Dimension(
                        "dim_dose",
                        values={"IVDOSE_icg": Q(np.linspace(10, 500, num=20), "mg")},
                    ),
                ],
            )
        }

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk(window=self.window)
        return {**self._figure_icgkel_vs_dose()}

    def _pk_band(
        self, df: pd.DataFrame, key: str, unit: str
    ) -> tuple[Quantity, Quantity | None, Quantity | None]:
        """Get a parameter at the points of a scan, the median and the band of the draws.

        Returns:
            The values (the median of the draws in a run over them), and the
            lower and the upper quantile, `None` without draws.
        """
        source = df[f"{key}_unit"].iloc[0]
        run = self.uncertainty
        if run is None:
            return Q(df[key].to_numpy(dtype=float), source).to(unit), None, None
        summary = summarize(df, key, run.quantiles)
        median, lower, upper = (
            Q(summary[column].to_numpy(), source).to(unit)
            for column in ["median", "lower", "upper"]
        )
        return median, lower, upper

    def _figure_icgkel_vs_dose(self) -> dict[str, FigureMPL]:
        """Elimination rate of indocyanine green vs dose."""
        df = self.pk_dfs["scan_dose"]

        f, ax = base_figure()
        ax.set_xlabel("Dose [mg/kg]", fontdict=self.font)
        ax.set_ylabel(f"{self.pk_labels['kel']} (5-12 min) [1/min]", fontdict=self.font)

        # simulation: [mmole] -> [mg/kg], the dose of every point of the scan
        doses = df.groupby("point")["dose"].first().to_numpy()
        x = Q(doses, df["dose_unit"].iloc[0]) * self.Mr.icg / Q(self.bodyweight, "kg")
        x = x.to("mg/kg")
        y, lower, upper = self._pk_band(df, "kel", "1/min")
        if lower is not None and upper is not None:
            ax.fill_between(
                x.magnitude, lower.magnitude, upper.magnitude, color="black", alpha=0.2
            )
        ax.plot(
            x.magnitude,
            y.magnitude,
            marker="o",
            linestyle="-",
            color="black",
            label="Sim: dose scan (75 kg)",
            ms=4,
        )

        ax.set_ylim(bottom=0, top=0.3)

        # data
        dset = self._datasets["Tab2_ICGkel_ICGControls"]
        ax.errorbar(
            dset.get_quantity("dose").magnitude,
            dset.get_quantity("mean").magnitude,
            yerr=dset["mean_sd"].values,
            marker="s",
            linestyle="",
            color="black",
            label="Martin 1976 controls (n=19)",
        )
        ax.legend(fontsize=self.legend_font_size, loc="lower left")

        return {"fig_kel": f}
