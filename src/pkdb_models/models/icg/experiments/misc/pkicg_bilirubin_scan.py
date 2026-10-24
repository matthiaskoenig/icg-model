"""Dependency of the pharmacokinetic parameters of indocyanine green on bilirubin."""

from typing import override

from matplotlib.figure import Figure
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class ICGPKBilirubinScan(ICGSimulationExperiment):
    """Dependency of the pharmacokinetic parameters on the plasma bilirubin.

    The bilirubin concentration is scanned in every stage of cirrhosis.
    """

    @override
    def simulations(self) -> dict[str, Scan]:
        tcs = self.cirrhosis_timecourses()
        return {
            self.to_sid(f"scan_bilirubin_{key}"): Scan(
                simulation=tcs[key],
                dimensions=[
                    Dimension(
                        "dim_bil",
                        values=self.bilirubin_scan_changes(start=-3, stop=-0.5),
                    )
                ],
            )
            for key in self.cirrhosis_map
        }

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()
        return self.fig_icgpk_vs_bil()

    def fig_icgpk_vs_bil(self) -> dict[str, Figure]:
        """Plot the pharmacokinetic parameters against the plasma bilirubin.

        Returns:
            The figures by id, `fig_bil_cir<parameter>`.
        """
        figures = {}

        for pk_key in self.pk_scan_parameters:
            f, ax = base_figure()

            ax.axvline(x=0.02, color="grey", linestyle="--")
            for key in self.cirrhosis_map:
                sim_key = f"scan_bilirubin_{key}"
                xres = self.results[self.to_sid(f"task_{sim_key}")]
                df = self.pk_dfs[self.to_sid(sim_key)]

                x = self.scan_values(xres, "[LI__bil_ext]", self.unit_bilirubin)
                self.plot_pk(
                    ax,
                    x,
                    df,
                    pk_key,
                    self.pk_units[pk_key],
                    color=self.cirrhosis_colors[key],
                    label=key,
                )

            ax.set_xlabel(
                f"Plasma bilirubin [{self.unit_bilirubin}]", fontdict=self.font
            )
            ax.set_ylabel(
                f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                fontdict=self.font,
            )

            ax.set_xscale("log")
            ax.set_ylim(top=self.pk_limits[pk_key])
            ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_bil_cir{pk_key}"] = f
        return figures
