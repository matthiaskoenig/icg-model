"""Dependency of the pharmacokinetic parameters of indocyanine green on OATP1B3."""

from typing import override

import numpy as np
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class ICGPKOatp1b3Scan(ICGSimulationExperiment):
    """Dependency of the pharmacokinetic parameters on the OATP1B3 activity.

    The activity of the hepatic uptake transporter OATP1B3 is scanned in every
    stage of cirrhosis and compared to the data of the genotypes of Anzai2020,
    Kagawa2017 and the excretory defect of Namihisa1981.
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig1", study="Anzai2020")
        for measurement in df.measurement.unique():
            for group in df.group.unique():
                df_subset = df.copy()[
                    (df.measurement == measurement) & (df.group == group)
                ]
                dset = DataSet.from_df(df_subset, ureg=self.ureg)
                dsets[self.to_sid(f"Anzai_{measurement}_{group}")] = dset

        df = self.load_dataframe("Tab1", study="Kagawa2017")
        dsets["Kagawa_r15"] = DataSet.from_df(df, ureg=self.ureg)

        return dsets

    @override
    def simulations(self) -> dict[str, Scan]:
        tcs = self.cirrhosis_timecourses()
        return {
            self.to_sid(f"scan_f_oatp1b3_{key}"): Scan(
                simulation=tcs[key],
                dimensions=[
                    Dimension(
                        "dim_oatp1b3",
                        values={
                            "LI__f_oatp1b3": Q(
                                np.linspace(0.1, 1.9, num=20), "dimensionless"
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
        return self.fig_icgpk_vs_oatp1b3()

    def fig_icgpk_vs_oatp1b3(self) -> dict[str, Figure]:
        """Plot the pharmacokinetic parameters against the OATP1B3 activity.

        Returns:
            The figures by id, `fig_<parameter>`.
        """
        figures = {}

        for pk_key in self.pk_scan_parameters:
            f, ax = base_figure()
            ax.axvline(x=1.0, color="grey", linestyle="--")

            for key in self.cirrhosis_map:
                sim_key = f"scan_f_oatp1b3_{key}"
                xres = self.results[self.to_sid(f"task_{sim_key}")]
                df = self.pk_dfs[self.to_sid(sim_key)]

                x = self.scan_values(xres, "LI__f_oatp1b3")
                self.plot_pk(
                    ax,
                    x,
                    df,
                    pk_key,
                    self.pk_units[pk_key],
                    color=self.cirrhosis_colors[key],
                    label=key,
                )
                ax.set_xlabel("OATP1B3 [-]", fontdict=self.font)
                ax.set_ylabel(
                    f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                    fontdict=self.font,
                )
            ax.set_ylim(top=self.pk_limits[pk_key])
            if pk_key == "icgr15":
                # Anzai 2020
                wy = self._datasets["Anzai_retention_ratio_wildtype"].get_quantity(
                    "value"
                )
                hety = self._datasets[
                    "Anzai_retention_ratio_heterozygote"
                ].get_quantity("value")
                homoy = self._datasets["Anzai_retention_ratio_homozygote"].get_quantity(
                    "value"
                )
                ax.errorbar(
                    1,
                    np.mean(wy).magnitude,
                    yerr=np.std(wy).magnitude,
                    linestyle="",
                    marker="s",
                    color="black",
                    ms=5,
                    label="Anzai 2020 wildtype",
                    zorder=100,
                )
                ax.errorbar(
                    0.5,
                    np.mean(hety).magnitude,
                    yerr=np.std(hety).magnitude,
                    linestyle="",
                    marker="D",
                    color="black",
                    ms=5,
                    label="Anzai 2020 heterozygote",
                )
                ax.plot(
                    0.0,
                    homoy.magnitude,
                    linestyle="",
                    marker="*",
                    color="black",
                    ms=10,
                    label="Anzai 2020 homozygote",
                )
                # Kagawa 2017
                kagawa_y = self._datasets["Kagawa_r15"].get_quantity("value")
                ax.errorbar(
                    0.0,
                    np.mean(kagawa_y).magnitude,
                    yerr=np.std(kagawa_y).magnitude,
                    linestyle="",
                    marker="^",
                    color="black",
                    ms=7,
                    label="Kagawa 2017 homozygote",
                )

            if pk_key == "icgpdr":
                ax.errorbar(
                    x=0.0,
                    y=2.1,
                    yerr=0.4,
                    linestyle="",
                    marker="",
                    color="black",
                    capsize=4,
                    label="Namihisa 1981 excretory defect",
                )
            ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_{pk_key}"] = f
        return figures
