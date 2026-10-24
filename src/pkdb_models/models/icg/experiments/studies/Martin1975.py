"""Simulation experiment of the study Martin1975."""

from typing import ClassVar, override

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure
from pkdb_models.models.icg.icg_pk import summarize


class Martin1975(ICGSimulationExperiment):
    """Simulation experiment of Martin1975.

    Dose dependency of the elimination rate of indocyanine green in healthy
    young women and men (18 - 29 years): 0.5 mg/kg in 51 subjects, 2.0, 3.5
    and 5.0 mg/kg in 19 of them. The elimination rate `k` of the paper is the
    log-linear regression of the plasma concentrations from 5 to 12 min after
    the end of the injection (`window`), the clearance (Fig 2) is `k` times
    the volume of distribution of the back extrapolation (`cl_kel`). The body
    weight is not reported, the model's 75 kg are used. The data are mean and
    standard error.
    """

    #: the window of the log-linear regression of the paper [min]
    window: ClassVar[tuple[float, float]] = (5.0, 12.0)
    #: the doses of the scan [mg/kg]
    doses_perkg: ClassVar[np.ndarray] = np.linspace(0.1, 6.0, num=20)
    #: the body weight of the simulations [kg]
    bodyweight: ClassVar[float] = 75.0

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for df in self.load_dataframes_by_substance("Tab1").values():
            dset = DataSet.from_df(df, self.ureg)
            for group in df.group.unique():
                dsets[self.to_sid(f"Tab1_ICGkel_{group}")] = dset[dset.group == group]
        df = self.load_dataframe("Fig2")
        dsets["Fig2_ICGcl"] = DataSet.from_df(df, self.ureg)
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
                        values={
                            "IVDOSE_icg": Q(self.doses_perkg * self.bodyweight, "mg")
                        },
                    ),
                ],
            )
        }

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk(window=self.window)
        return {**self._figure_icgkel_vs_dose(), **self._figure_cl_vs_dose()}

    def _plot_scan(
        self, ax: Axes, key: str, unit: str, per: str = "dimensionless"
    ) -> None:
        """Plot a parameter of the dose scan, with the band of the draws.

        Args:
            ax: the axes.
            key: the parameter of `calculate_icg_pk`.
            unit: the unit of the figure.
            per: the unit of the body weight for a parameter per body weight
                (`kg`), the parameter is divided by the body weight.
        """
        df: pd.DataFrame = self.pk_dfs["scan_dose"]
        unit_pk = df[f"{key}_unit"].iloc[0]
        divisor = Q(self.bodyweight, per) if per != "dimensionless" else Q(1, per)

        def convert(values: np.ndarray) -> np.ndarray:
            """Get the values in the unit of the figure."""
            return (Q(values, unit_pk) / divisor).to(unit).magnitude

        x = self.doses_perkg
        if self.uncertainty is None:
            y = convert(df[key].to_numpy(dtype=float))
        else:
            summary = summarize(df, key, self.uncertainty.quantiles)
            y, lower, upper = (
                convert(summary[column].to_numpy())
                for column in ["median", "lower", "upper"]
            )
            ax.fill_between(x, lower, upper, color="black", alpha=0.2, linewidth=0)
        ax.plot(
            x, y, marker="o", linestyle="-", color="black", label="Scan: dose", ms=4
        )

    def _figure_icgkel_vs_dose(self) -> dict[str, FigureMPL]:
        """Elimination rate of indocyanine green (5 - 12 min) vs dose."""
        f, ax = base_figure()
        ax.set_xlabel("Dose [mg/kg]", fontdict=self.font)
        ax.set_ylabel(f"{self.pk_labels['kel']} [1/min]", fontdict=self.font)
        self._plot_scan(ax, key="kel", unit="1/min")

        # data: marker and label of the groups
        groups = {
            "allWomen": ("s", "Martin 1975 women (n=33)"),
            "allMen": ("*", "Martin 1975 men (n=18)"),
            "additionalWomen": ("D", "Martin 1975 women (n=10)"),
            "additionalMen": ("^", "Martin 1975 men (n=9)"),
        }
        for group, (marker, label) in groups.items():
            dset = self._datasets[self.to_sid(f"Tab1_ICGkel_{group}")]
            ax.errorbar(
                dset.get_quantity("dose").magnitude,
                dset.get_quantity("mean").magnitude,
                yerr=dset["mean_se"].values,
                marker=marker,
                linestyle="",
                color="black",
                label=label,
            )

        ax.set_ylim(bottom=0, top=0.4)
        ax.legend(fontsize=self.legend_font_size)
        return {"fig_kel": f}

    def _figure_cl_vs_dose(self) -> dict[str, FigureMPL]:
        """Clearance of indocyanine green per body weight (k * V) vs dose."""
        f, ax = base_figure()
        ax.set_xlabel("Dose [mg/kg]", fontdict=self.font)
        ax.set_ylabel(f"{self.pk_labels['clperbw']} [ml/min/kg]", fontdict=self.font)
        self._plot_scan(ax, key="cl_kel", unit="ml/min/kg", per="kg")

        # data: the clearance after 0.5 mg/kg
        dset = self._datasets["Fig2_ICGcl"]
        ax.errorbar(
            [0.5],
            dset.get_quantity("mean").to("ml/min/kg").magnitude,
            yerr=dset["mean_se"].values,
            marker="s",
            linestyle="",
            color="black",
            label="Martin 1975 (n=51)",
        )
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=self.legend_font_size)
        return {"fig_cl": f}
