"""Dependency of the pharmacokinetic parameters of indocyanine green on the dose."""

from typing import override

import numpy as np
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class ICGPKDoseScan(ICGSimulationExperiment):
    """Dependency of the pharmacokinetic parameters on the intravenous dose.

    The dose is scanned in every stage of cirrhosis.
    """

    @override
    def simulations(self) -> dict[str, Scan]:
        # the scan sets the dose
        tcs = self.cirrhosis_timecourses(dose=None)
        return {
            self.to_sid(f"scan_dose_{key}"): Scan(
                simulation=tcs[key],
                dimensions=[
                    Dimension(
                        "dim_dose",
                        values={"IVDOSE_icg": Q(np.logspace(-2, 3, num=20), "mg")},
                    )
                ],
            )
            for key in self.cirrhosis_map
        }

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()
        return self.fig_icgpk_vs_dose()

    def fig_icgpk_vs_dose(self) -> dict[str, Figure]:
        """Plot the pharmacokinetic parameters against the dose.

        Returns:
            The figures by id, `fig_dose_cir<parameter>`.
        """
        figures = {}

        for pk_key in self.pk_scan_parameters:
            f, ax = base_figure()
            if pk_key == "cl_kel":
                ax.set_ylim(top=830)
            if pk_key == "icgpdr":
                ax.set_ylim(top=30)

            for key in self.cirrhosis_map:
                # get simulation results
                sim_key = f"scan_dose_{key}"
                xres = self.results[self.to_sid(f"task_{sim_key}")]
                df = self.pk_dfs[self.to_sid(sim_key)]

                x = self.scan_values(xres, "IVDOSE_icg", self.unit_dose)
                self.plot_pk(
                    ax,
                    x,
                    df,
                    pk_key,
                    self.pk_units[pk_key],
                    color=self.cirrhosis_colors[key],
                    label=key,
                )
                ax.set_xscale("log")

                ax.set_xlabel(f"Dose [{self.unit_dose}]", fontdict=self.font)
                ax.set_ylabel(
                    f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                    fontdict=self.font,
                )

            ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_dose_cir{pk_key}"] = f
        return figures
