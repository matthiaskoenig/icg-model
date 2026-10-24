"""Correlation between the pharmacokinetic parameters of indocyanine green."""

from typing import override

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.simulation import Dimension, Scan, Simulation

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment


class ICGPKcorrelation(ICGSimulationExperiment):
    """Correlation between the pharmacokinetic parameters of indocyanine green.

    The parameters of a scan of the dose and of a scan of the body weight are
    plotted against each other.
    """

    @override
    def simulations(self) -> dict[str, Scan]:
        return {
            "scan_dose": Scan(
                simulation=Simulation(
                    time_unit="min",
                    end=30,
                    preinit_changes={**self.default_changes()},
                ),
                dimensions=[
                    Dimension(
                        "dim_dose",
                        values={"IVDOSE_icg": Q(np.logspace(-2, 3, num=20), "mg")},
                    )
                ],
            ),
            "scan_bodyweight": Scan(
                simulation=Simulation(
                    time_unit="min",
                    end=30,
                    preinit_changes={
                        **self.default_changes(),
                        "IVDOSE_icg": Q(37.5, "mg"),
                    },
                ),
                dimensions=[
                    Dimension(
                        "dim_bodyweight",
                        values={"BW": Q(np.linspace(67.5, 82.5, num=3), "kg")},
                    )
                ],
            ),
        }

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()
        return self.fig_icgpk_vs_pk()

    def fig_icgpk_vs_pk(self) -> dict[str, Figure]:
        """Plot the parameters of the liver function test against each other.

        In a run over the draws of the fitted parameters the medians of the
        draws of every point.

        Returns:
            The figures by id, `fig_<parameter>_vs_<parameter>`.
        """
        figures = {}

        plotmap = {
            "icgpdr": ["thalf_kel", "cl_kel", "icgr15"],
            "thalf_kel": ["cl_kel", "icgr15"],
            "cl_kel": ["icgr15"],
        }

        for pk_key, corr_keys in plotmap.items():
            for corr_key in corr_keys:
                f, ax = plt.subplots()
                f.subplots_adjust(wspace=0.3, hspace=0.3)

                for sim_key, df in self.pk_dfs.items():
                    x, _, _ = self.pk_band(df, corr_key, self.pk_units[corr_key])
                    y, _, _ = self.pk_band(df, pk_key, self.pk_units[pk_key])
                    label_token = "Dose" if sim_key == "scan_dose" else "Bodyweight"
                    ax.plot(
                        x,
                        y,
                        label=f"Scan: {label_token}",
                        marker="s",
                        linestyle="-",
                    )
                ax.set_xlabel(
                    f"{self.pk_labels[corr_key]} [{self.pk_units[corr_key]}]",
                    fontdict=self.font,
                )
                ax.set_ylabel(
                    f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                    fontdict=self.font,
                )
                ax.legend()

                figures[f"fig_{pk_key}_vs_{corr_key}"] = f
        return figures
