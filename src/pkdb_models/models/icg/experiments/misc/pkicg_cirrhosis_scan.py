"""Dependency of the pharmacokinetic parameters of indocyanine green on cirrhosis."""

from typing import ClassVar, override

import numpy as np
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure, scan_color


class ICGPKCirrhosisScan(ICGSimulationExperiment):
    """Dependency of the pharmacokinetic parameters on the degree of cirrhosis.

    Cirrhosis is the loss of functional tissue together with the shunts around
    the liver; both are scanned together and alone.
    """

    colors: ClassVar[dict[str, str]] = {
        "scan_shunts": "green",
        "scan_tissue_loss": "orange",
    }
    sampling: ClassVar[np.ndarray] = np.linspace(0, 0.9, num=20)

    @override
    def simulations(self) -> dict[str, Scan]:
        return {
            "scan_cirrhosis": Scan(
                simulation=self.normal_timecourse(),
                dimensions=[
                    Dimension(
                        "dim_cirrhosis",
                        values={
                            **self.cirrhosis_scan_changes(),
                        },
                    ),
                ],
            ),
            "scan_shunts": Scan(
                simulation=self.normal_timecourse(),
                dimensions=[
                    Dimension(
                        "dim_shunts",
                        values={"f_shunts": Q(self.sampling, "dimensionless")},
                    ),
                ],
            ),
            "scan_tissue_loss": Scan(
                simulation=self.normal_timecourse(),
                dimensions=[
                    Dimension(
                        "dim_tissue_loss",
                        values={"f_tissue_loss": Q(self.sampling, "dimensionless")},
                    ),
                ],
            ),
        }

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()
        return self.fig_icgpk_vs_cirrhosis()

    def fig_icgpk_vs_cirrhosis(self) -> dict[str, Figure]:
        """Plot the parameters of the liver function test against the degree of cirrhosis.

        Returns:
            The figures by id, `fig_<parameter>`.
        """
        cmap = scan_color(9)
        figures = {}
        x = self.sampling

        for pk_key in self.pk_scan_parameters:
            f, ax = base_figure()
            unit = self.pk_units[pk_key]

            for sim_key, df in self.pk_dfs.items():
                label_token = " ".join(sim_key.split("_")[1:])
                if sim_key == "scan_cirrhosis":
                    y = self.plot_pk(
                        ax,
                        x,
                        df,
                        pk_key,
                        unit,
                        color="black",
                        label=f"Scan: {label_token}",
                        markersize=0,
                    )
                    ax.scatter(x, y, c=x, cmap=cmap, zorder=50)
                    x_stars = list(self.cirrhosis_map.values())
                    y_stars = np.interp(x_stars, x, y)
                    for k, (key, color) in enumerate(self.cirrhosis_colors.items()):
                        ax.plot(
                            x_stars[k],
                            y_stars[k],
                            marker="*",
                            linestyle="",
                            color=color,
                            label=key,
                            zorder=100,
                            ms=20,
                        )
                else:
                    self.plot_pk(
                        ax,
                        x,
                        df,
                        pk_key,
                        unit,
                        color=self.colors[sim_key],
                        label=f"Scan: {label_token}",
                        linestyle="--",
                        linewidth=1,
                        markersize=4,
                    )
            ax.set_xlabel("Degree of cirrhosis [-]", fontdict=self.font)
            ax.set_ylabel(
                f"{self.pk_labels[pk_key]} [{unit}]",
                fontdict=self.font,
            )
            ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_{pk_key}"] = f
        return figures
