"""Dependency of the pharmacokinetic parameters of indocyanine green on hepatectomy."""

from typing import override

import numpy as np
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class ICGPKHepatectomyScan(ICGSimulationExperiment):
    """Dependency of the pharmacokinetic parameters on the resection rate.

    The resection rate of the liver is scanned in every stage of cirrhosis.
    """

    @override
    def simulations(self) -> dict[str, Scan]:
        tcs = self.cirrhosis_timecourses()
        return {
            self.to_sid(f"scan_rr_{key}"): Scan(
                simulation=tcs[key],
                dimensions=[
                    Dimension(
                        "dim_rr",
                        values={
                            "resection_rate": Q(
                                np.linspace(0, 0.9, num=20), "dimensionless"
                            )
                        },
                    ),
                ],
            )
            for key in self.cirrhosis_map
        }

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()
        return self.fig_icgpk_vs_resection()

    def fig_icgpk_vs_resection(self) -> dict[str, Figure]:
        """Plot the pharmacokinetic parameters against the resection rate.

        Returns:
            The figures by id, `fig_<parameter>`.
        """
        figures = {}

        for pk_key in self.pk_scan_parameters:
            f, ax = base_figure()

            for key in self.cirrhosis_map:
                sim_key = f"scan_rr_{key}"
                xres = self.results[self.to_sid(f"task_{sim_key}")]
                df = self.pk_dfs[self.to_sid(sim_key)]

                x = self.scan_values(xres, "resection_rate")
                self.plot_pk(
                    ax,
                    x,
                    df,
                    pk_key,
                    self.pk_units[pk_key],
                    color=self.cirrhosis_colors[key],
                    label=key,
                )
                ax.set_xlabel("Resection rate [-]", fontdict=self.font)
                ax.set_ylabel(
                    f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                    fontdict=self.font,
                )
                ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_{pk_key}"] = f
        return figures
