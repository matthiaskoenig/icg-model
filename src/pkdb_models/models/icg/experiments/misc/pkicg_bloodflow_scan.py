"""Dependency of the pharmacokinetic parameters of indocyanine green on blood flow."""

from typing import override

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class ICGPKBloodFlowScan(ICGSimulationExperiment):
    """Dependency of the pharmacokinetic parameters on the blood flow.

    The hepatic blood flow and the cardiac output are scanned in every stage of
    cirrhosis and compared to the data of Grundmann1992, Skak1987, Pind2016,
    Moeller1998 and Moeller2018.
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Tab1", study="Grundmann1992")
        for measurement_type in df.measurement_type.unique():
            if measurement_type in ["thalf", "clearance", "hbf", "cardiac output"]:
                df_subset = df.copy()[df.measurement_type == measurement_type]
                dset = DataSet.from_df(df_subset, self.ureg)
                dsets[self.to_sid(f"Grundmann_{measurement_type}")] = dset

        df = self.load_dataframe("Tab2", study="Pind2016")
        for measurement in df.measurement.unique():
            if measurement in ["retention ratio", "clearance", "hbf", "cardiac output"]:
                df_subset = df.copy()[df.measurement == measurement]
                dset = DataSet.from_df(df_subset, self.ureg)
                dsets[self.to_sid(f"Pind_{measurement}")] = dset

        df = self.load_dataframe("Tab2", study="Moeller1998")
        for measurement in df.measurement.unique():
            df_subset = df.copy()[df.measurement == measurement]
            dset = DataSet.from_df(df_subset, self.ureg)
            dsets[self.to_sid(f"Moeller1998_{measurement}")] = dset

        df = self.load_dataframe("Tab2", study="Moeller2018")
        for measurement in df.measurement.unique():
            df_subset = df.copy()[
                (df.measurement == measurement) & (df.group != "Controls")
            ]
            dset = DataSet.from_df(df_subset, self.ureg)
            dsets[self.to_sid(f"Moeller2018_{measurement}")] = dset

        df = self.load_dataframe("Tab3", study="Skak1987")
        for measurement in df.measurement.unique():
            df_subset = df.copy()[
                (df.measurement == measurement) & (df.group == "healthy")
            ]
            dset = DataSet.from_df(df_subset, self.ureg)
            dsets[self.to_sid(f"Skak_{measurement}")] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Scan]:
        tcscans = {}
        tcs = self.cirrhosis_timecourses()
        for key in self.cirrhosis_map:
            tc = tcs[key]

            # hbf scan
            tcscans[self.to_sid(f"scan_f_bloodflow_{key}")] = Scan(
                simulation=tc,
                dimensions=[
                    Dimension("dim_hbf", values=self.bloodflow_scan_changes(stop=3)),
                ],
            )

            # cardiac output scan
            tcscans[self.to_sid(f"scan_f_cardiac_output_{key}")] = Scan(
                simulation=tc,
                dimensions=[
                    Dimension(
                        "dim_cardiacoutput",
                        values={
                            "f_cardiac_output": Q(
                                np.linspace(0.2, 3.0, num=20), "dimensionless"
                            )
                        },
                    ),
                ],
            )
        return tcscans

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()

        return {
            **self.fig_icgpk_vs_bloodflow_cirrhosis(),
            **self.fig_icgpk_vs_cardiac_output_cirrhosis(),
        }

    def plot_ctp_data(
        self,
        ax: Axes,
        xdset: DataSet,
        ydset: DataSet,
        label: str,
        marker: str,
        factor: float | None = None,
    ) -> None:
        """Plot the mean and SD of the groups CTP A-C of a study.

        Args:
            ax: the axes to plot in.
            xdset: the data of the x-axis, the groups CTP A-C.
            ydset: the data of the y-axis, the groups CTP A-C.
            label: the label of the data.
            marker: the marker of the data.
            factor: factor of the unit conversion of the y data.
        """
        x = xdset["mean"].to_numpy()
        x_sd = xdset["mean_sd"].to_numpy()
        y = ydset["mean"].to_numpy()
        y_sd = ydset["mean_sd"].to_numpy()
        if factor is not None:
            y = y * factor
            y_sd = y_sd * factor
        for k, color in enumerate(list(self.cirrhosis_colors.values())[1:]):
            ax.errorbar(
                x[k],
                y[k],
                xerr=x_sd[k],
                yerr=y_sd[k],
                marker=marker,
                color=color,
                label=label if k == 0 else None,
                linestyle="",
            )

    def plot_grundmann1992(self, ax: Axes, pk_key: str, xdset: DataSet) -> None:
        """Plot the data of the controls of Grundmann1992.

        Args:
            ax: the axes to plot in.
            pk_key: the pharmacokinetic parameter, `thalf_kel` or `cl_kel`.
            xdset: the data of the x-axis.
        """
        ydset = self._datasets[
            self.to_sid(
                {"thalf_kel": "Grundmann_thalf", "cl_kel": "Grundmann_clearance"}[
                    pk_key
                ]
            )
        ]
        ax.errorbar(
            xdset["mean"].to_numpy(),
            ydset["mean"].to_numpy(),
            xerr=xdset["mean_sd"].to_numpy(),
            yerr=ydset["mean_sd"].to_numpy(),
            marker="s",
            color="black",
            label="Grundmann 1992 controls",
            linestyle="",
        )

    def plot_pind2016(self, ax: Axes, pk_key: str, xdset: DataSet) -> None:
        """Plot the data of the groups CTP A-C of Pind2016.

        Args:
            ax: the axes to plot in.
            pk_key: the pharmacokinetic parameter, `icgr15` or `cl_kel`.
            xdset: the data of the x-axis.
        """
        if pk_key == "icgr15":
            self.plot_ctp_data(
                ax,
                xdset=xdset,
                ydset=self._datasets["Pind_retention_ratio"],
                label="Pind 2016 CTP A-C",
                marker="s",
                factor=100,  # [-] -> [%]
            )
        else:
            self.plot_ctp_data(
                ax,
                xdset=xdset,
                ydset=self._datasets["Pind_clearance"],
                label="Pind 2016 CTP A-C",
                marker="s",
            )

    def fig_icgpk_vs_bloodflow_cirrhosis(self) -> dict[str, Figure]:
        """Plot the pharmacokinetic parameters against the hepatic blood flow.

        Returns:
            The figures by id, `fig_bf_cir<parameter>`.
        """
        figures = {}

        for pk_key in self.pk_scan_parameters:
            f, ax = base_figure()
            # the hepatic blood flow of the model at 75 kg
            ax.axvline(
                x=self.hepatic_bloodflow().magnitude, color="grey", linestyle="--"
            )

            for key in self.cirrhosis_map:
                sim_key = f"scan_f_bloodflow_{key}"
                xres = self.results[self.to_sid(f"task_{sim_key}")]
                df = self.pk_dfs[self.to_sid(sim_key)]

                x = self.scan_values(xres, "Qh", self.unit_bloodflow)
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
                f"Hepatic blood flow [{self.unit_bloodflow}]", fontdict=self.font
            )
            ax.set_ylabel(
                f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                fontdict=self.font,
            )

            # add data from Grundmann1992
            if pk_key in ["thalf_kel", "cl_kel"]:
                self.plot_grundmann1992(
                    ax, pk_key=pk_key, xdset=self._datasets["Grundmann_hbf"]
                )

            # add data from Skak1987
            if pk_key == "cl_kel":
                x_skak = self._datasets["Skak_hbf"].get_quantity("value")
                y_skak = self._datasets["Skak_clearance"].get_quantity("value")
                # the hepatic plasma flow to the blood flow with the
                # hematocrit of the model
                hct = float(np.asarray(xres["HCT"].values).flat[0])
                x_mean = np.mean(x_skak) / (1 - hct)
                y_mean = np.mean(y_skak) * 1000
                x_sd = np.std(x_skak) / (1 - hct)
                y_sd = np.std(y_skak) * 1000
                ax.errorbar(
                    x_mean.magnitude,
                    y_mean.magnitude,
                    xerr=x_sd.magnitude,
                    yerr=y_sd.magnitude,
                    marker="^",
                    color="black",
                    label="Skak 1987 controls",
                    linestyle="",
                )

            # add data from Pind2016
            if pk_key in ["icgr15", "cl_kel"]:
                self.plot_pind2016(ax, pk_key=pk_key, xdset=self._datasets["Pind_hbf"])

            # add data from Moeller1998
            if pk_key == "cl_kel":
                self.plot_ctp_data(
                    ax,
                    xdset=self._datasets["Moeller1998_hbf"],
                    ydset=self._datasets["Moeller1998_clearance"],
                    label="Moeller 1998 CTP A-C",
                    marker="D",
                    factor=1000,  # l -> ml
                )

            # add data from Moeller 2018
            if pk_key in ["icgr15", "thalf_kel", "icgpdr"]:
                self.plot_ctp_data(
                    ax,
                    xdset=self._datasets["Moeller2018_hbf"],
                    ydset=self._datasets[
                        {
                            "icgr15": "Moeller2018_retention_ratio",
                            "thalf_kel": "Moeller2018_thalf",
                            "icgpdr": "Moeller2018_kel",
                        }[pk_key]
                    ],
                    label="Moeller 2018 CTP A-C",
                    marker="^",
                    # 1/min -> %/min
                    factor=100 if pk_key == "icgpdr" else None,
                )

            ax.set_ylim(top=self.pk_limits[pk_key])
            if pk_key == "cl_kel":
                ax.set_ylim(top=1550)
            ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_bf_cir{pk_key}"] = f
        return figures

    def fig_icgpk_vs_cardiac_output_cirrhosis(self) -> dict[str, Figure]:
        """Plot the pharmacokinetic parameters against the cardiac output.

        Returns:
            The figures by id, `fig_co_cir<parameter>`.
        """
        figures = {}

        for pk_key in self.pk_scan_parameters:
            f, ax = base_figure()
            # the cardiac output of the model at 75 kg (75 ml/min/kg)
            ax.axvline(
                x=self.hepatic_bloodflow().magnitude / 0.255,
                color="grey",
                linestyle="--",
            )

            for key in self.cirrhosis_map:
                # get simulation results
                sim_key = f"scan_f_cardiac_output_{key}"
                xres = self.results[self.to_sid(f"task_{sim_key}")]
                df = self.pk_dfs[self.to_sid(sim_key)]

                x = self.scan_values(xres, "QC", self.unit_bloodflow)
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
                    f"Cardiac output [{self.unit_bloodflow}]", fontdict=self.font
                )
                ax.set_ylabel(
                    f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                    fontdict=self.font,
                )

            # add data from Grundmann1992
            if pk_key in ["thalf_kel", "cl_kel"]:
                self.plot_grundmann1992(
                    ax, pk_key=pk_key, xdset=self._datasets["Grundmann_cardiac_output"]
                )

            # add data from Pind2016
            if pk_key in ["icgr15", "cl_kel"]:
                self.plot_pind2016(
                    ax, pk_key=pk_key, xdset=self._datasets["Pind_cardiac_output"]
                )

            ax.set_ylim(top=self.pk_limits[pk_key])
            ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_co_cir{pk_key}"] = f
        return figures
