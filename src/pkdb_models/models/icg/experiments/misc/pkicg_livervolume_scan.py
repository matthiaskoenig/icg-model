"""Dependency of the pharmacokinetic parameters of indocyanine green on liver volume."""

from typing import override

import numpy as np
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class ICGPKLiverVolumeScan(ICGSimulationExperiment):
    """Dependency of the pharmacokinetic parameters on the liver volume.

    The liver volume is scanned in every stage of cirrhosis and compared to the
    data of Roberts1976, Haimerl2016 and Hashimoto2000.
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig1B", study="Roberts1976")
        for measurement in df.measurement.unique():
            df_subset = df.copy()[df.measurement == measurement]
            dset = DataSet.from_df(df_subset, self.ureg)
            dsets[self.to_sid(f"Roberts1976_{measurement}")] = dset

        for sid in ["Hashimoto2000_Tab1", "Haimerl2016_Tab1"]:
            study, fig_id = sid.split("_", 1)
            df = self.load_dataframe(fig_id, study=study)
            for measurement in df.measurement.unique():
                df_subset = df.copy()[
                    (df.measurement == measurement) & (df.group != "all")
                ]
                dset = DataSet.from_df(df_subset, self.ureg)
                dsets[self.to_sid(f"{sid}_{measurement}")] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Scan]:
        tcs = self.cirrhosis_timecourses()
        return {
            self.to_sid(f"scan_FVli_{key}"): Scan(
                simulation=tcs[key],
                dimensions=[
                    Dimension(
                        "dim_vli",
                        values={
                            "FVli": Q(np.linspace(0.1, 1.9, num=20) * 0.021, "l/kg")
                        },
                    ),
                ],
            )
            for key in self.cirrhosis_map
        }

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()
        return self.fig_icgpk_vs_livervolume()

    def fig_icgpk_vs_livervolume(self) -> dict[str, Figure]:
        """Plot the pharmacokinetic parameters against the liver volume.

        Returns:
            The figures by id, `fig_<parameter>`.
        """
        figures = {}

        for pk_key in self.pk_scan_parameters:
            f, ax = base_figure()
            ax.axvline(x=1.575, color="grey", linestyle="--")

            for key in self.cirrhosis_map:
                sim_key = f"scan_FVli_{key}"
                xres = self.results[self.to_sid(f"task_{sim_key}")]
                df = self.pk_dfs[self.to_sid(sim_key)]

                x = self.scan_values(xres, "Vli", "l")
                self.plot_pk(
                    ax,
                    x,
                    df,
                    pk_key,
                    self.pk_units[pk_key],
                    color=self.cirrhosis_colors[key],
                    label=key,
                )
                ax.set_xlabel("Liver volume [l]", fontdict=self.font)
                ax.set_ylabel(
                    f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                    fontdict=self.font,
                )

            if pk_key == "cl_kel":
                # the clearances of Roberts1976 (700 - 1700 ml/min) are taken
                # as blood clearances, the plasma clearance with the hematocrit
                # of the model
                hct = float(np.asarray(xres["HCT"].values).flat[0])
                x_data = self._datasets["Roberts1976_organ_volume"]["value"] / 1000
                y_data = self._datasets["Roberts1976_clearance"]["value"] * (1 - hct)
                ax.plot(
                    x_data,
                    y_data,
                    color="black",
                    linestyle="",
                    marker="s",
                    label="Roberts 1976 controls",
                )
                ax.errorbar(
                    np.mean(x_data),
                    np.mean(y_data),
                    xerr=np.std(x_data),
                    yerr=np.std(y_data),
                    linestyle="",
                    color="black",
                    marker="s",
                )

            if pk_key == "icgpdr":
                lv_dset = self._datasets["Haimerl2016_Tab1_liver_volume"]
                pdr_dset = self._datasets["Haimerl2016_Tab1_pdr"]
                x_data = lv_dset["mean"].to_numpy() / 1000  # ml->l
                x_sd = lv_dset["mean_sd"].to_numpy() / 1000  # ml->l
                y_data = pdr_dset["mean"].to_numpy()
                y_sd = pdr_dset["mean_sd"].to_numpy()
                ax.errorbar(
                    x_data[0],
                    y_data[0],
                    xerr=x_sd[0],
                    yerr=y_sd[0],
                    marker="s",
                    color="black",
                    label="Haimerl 2016 controls",
                    linestyle="",
                )
                for k, color in enumerate(list(self.cirrhosis_colors.values())[1:]):
                    ax.errorbar(
                        x_data[k + 1],
                        y_data[k + 1],
                        xerr=x_sd[k + 1],
                        yerr=y_sd[k + 1],
                        marker="s",
                        color=color,
                        label="Haimerl 2016 CTP A-C" if k == 0 else None,
                        linestyle="",
                    )

            if pk_key in ["icgpdr", "icgr15"]:
                lv_dset = self._datasets["Hashimoto2000_Tab1_liver_volume"]
                x_data = lv_dset["mean"].to_numpy() / 1000  # ml->l
                x_sd = lv_dset["mean_sd"].to_numpy() / 1000  # ml->l
                if pk_key == "icgpdr":
                    pdr_dset = self._datasets["Hashimoto2000_Tab1_kel"]
                    y_data = pdr_dset["mean"].to_numpy() * 100  # kel -> pdr
                    y_sd = pdr_dset["mean_sd"].to_numpy() * 100  # kel -> pdr
                else:
                    r15_dset = self._datasets["Hashimoto2000_Tab1_retention_ratio"]
                    y_data = r15_dset["mean"].to_numpy()
                    y_sd = r15_dset["mean_sd"].to_numpy()
                ax.errorbar(
                    x_data[0],
                    y_data[0],
                    xerr=x_sd[0],
                    yerr=y_sd[0],
                    marker="D",
                    color="black",
                    label="Hashimoto 2000 controls",
                    linestyle="",
                )
                ax.errorbar(
                    x_data[1],
                    y_data[1],
                    xerr=x_sd[1],
                    yerr=y_sd[1],
                    marker="D",
                    color="purple",
                    label="Hashimoto 2000 cirrhosis",
                    linestyle="",
                )

            ax.set_ylim(top=self.pk_limits[pk_key])
            if pk_key == "thalf_kel":
                ax.set_ylim(bottom=-5)

            ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_{pk_key}"] = f
        return figures
