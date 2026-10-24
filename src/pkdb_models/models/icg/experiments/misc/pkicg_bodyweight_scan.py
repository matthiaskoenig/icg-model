"""Dependency of the pharmacokinetic parameters of indocyanine green on body weight."""

from typing import override

import numpy as np
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class ICGPKBodyweightScan(ICGSimulationExperiment):
    """Dependency of the pharmacokinetic parameters on the body weight.

    The body weight is scanned in every stage of cirrhosis with a dose of
    0.5 mg/kg and compared to the data of Huet1983 and Haimerl2016.
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}

        df = self.load_dataframe("Tab1", study="Huet1983")
        for group in ["hepatitis", "cirrhosis"]:
            df_subset = df.copy()[(df.measurement == "clearance") & (df.group == group)]
            dsets[self.to_sid(f"Huet_clearance_{group}")] = DataSet.from_df(
                df_subset, ureg=self.ureg
            )
        df = self.load_dataframe("TabGroups", study="Huet1983")
        for group in ["hepatitis", "cirrhosis"]:
            df_subset = df.copy()[(df.measurement == "weight") & (df.name == group)]
            dsets[self.to_sid(f"Huet_weight_{group}")] = DataSet.from_df(
                df_subset, ureg=self.ureg
            )

        df = self.load_dataframe("Tab1", study="Haimerl2016")
        for measurement in df.measurement.unique():
            df_subset = df.copy()[(df.measurement == measurement) & (df.group != "all")]
            dset = DataSet.from_df(df_subset, self.ureg)
            dsets[self.to_sid(f"Haimerl_{measurement}")] = dset

        df = self.load_dataframe("TabGroups", study="Haimerl2016")
        df_subset = df.copy()[(df.measurement == "weight") & (df.name != "all")]
        dsets["Haimerl_weight"] = DataSet.from_df(df_subset, self.ureg)

        return dsets

    @override
    def simulations(self) -> dict[str, Scan]:
        # the scan sets the dose
        tcs = self.cirrhosis_timecourses(dose=None)
        bodyweights = Q(np.linspace(0.5, 1.5, num=20) * 75, "kg")
        return {
            self.to_sid(f"scan_BW_{key}"): Scan(
                simulation=tcs[key],
                dimensions=[
                    Dimension(
                        "dim_bw",
                        values={
                            "BW": bodyweights,
                            "IVDOSE_icg": bodyweights * Q(0.5, "mg/kg"),
                        },
                    ),
                ],
            )
            for key in self.cirrhosis_map
        }

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()
        return self.fig_icgpk_vs_bodyweight()

    def fig_icgpk_vs_bodyweight(self) -> dict[str, Figure]:
        """Plot the pharmacokinetic parameters against the body weight.

        Returns:
            The figures by id, `fig_<parameter>`.
        """
        figures = {}

        for pk_key in self.pk_scan_parameters:
            f, ax = base_figure()
            ax.axvline(x=75, color="grey", linestyle="--")

            for key in self.cirrhosis_map:
                sim_key = f"scan_BW_{key}"
                xres = self.results[self.to_sid(f"task_{sim_key}")]
                df = self.pk_dfs[self.to_sid(sim_key)]

                x = self.scan_values(xres, "BW", "kg")
                self.plot_pk(
                    ax,
                    x,
                    df,
                    pk_key,
                    self.pk_units[pk_key],
                    color=self.cirrhosis_colors[key],
                    label=key,
                )

            if pk_key == "cl_kel":
                group = "cirrhosis"
                bw_dset = self._datasets[self.to_sid(f"Huet_weight_{group}")]
                cl_dset = self._datasets[self.to_sid(f"Huet_clearance_{group}")]
                # the systemic blood clearance per kg of the study as the
                # plasma clearance with the hematocrit of the model
                hct = float(xres["HCT"].values.flat[0])
                bw = bw_dset["mean"].to_numpy()
                bw_se = bw_dset["mean_se"].to_numpy()
                cl = cl_dset["mean"].to_numpy() * bw * (1 - hct)
                cl_se = cl_dset["mean_se"].to_numpy() * bw * (1 - hct)
                ax.errorbar(
                    bw,
                    cl,
                    xerr=bw_se,
                    yerr=cl_se,
                    linestyle="",
                    marker="s",
                    color="purple",
                    label=rf"Huet 1983 {group} Mean$\pm$SE",
                )

            if pk_key == "icgpdr":
                bw_dset = self._datasets["Haimerl_weight"]
                pdr_dset = self._datasets["Haimerl_pdr"]
                x_data = bw_dset["mean"].to_numpy()
                x_sd = bw_dset["mean_sd"].to_numpy()
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

            ax.set_xlabel("Body weight [kg]", fontdict=self.font)
            ax.set_ylabel(
                f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                fontdict=self.font,
            )
            ax.set_ylim(top=self.pk_limits[pk_key])
            ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_{pk_key}"] = f
        return figures
