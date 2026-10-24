"""Mapping of the Child-Turcotte-Pugh score on the degree of cirrhosis."""

import logging
from typing import ClassVar, override

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure, scan_color

logger = logging.getLogger(__name__)


class CTPscoremapping(ICGSimulationExperiment):
    """Mapping of the CTP score on the degree of cirrhosis (`f_cirrhosis`).

    The pharmacokinetic parameters of a scan of the degree of cirrhosis are
    compared to the parameters of the CTP classes of Figg1995, Moeller1998,
    Moeller2018 and Herold2001. The mean body weight of Moeller1998 is 73 kg.
    """

    sampling: ClassVar[np.ndarray] = np.linspace(0, 0.9, num=50)

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for sid in [
            "Figg1995_Fig1mean",
            "Figg1995_Fig1",
            "Herold2001_Fig1mean",
            "Moeller1998_Tab2",
            "Moeller2018_Tab2",
        ]:
            study, fig_id = sid.split("_", 1)
            df = self.load_dataframe(fig_id, study=study)
            if sid == "Moeller2018_Tab2":
                for measurement in df.measurement.unique():
                    df_subset = df.copy()[df.measurement == measurement]
                    # FIXME: support multiple units in datasets !
                    dset = DataSet.from_df(df_subset, self.ureg)
                    dsets[self.to_sid(f"Moeller2018_Tab2_{measurement}")] = dset
            elif sid == "Moeller1998_Tab2":
                df_subset = df.copy()[df.measurement == "clearance"]
                dset = DataSet.from_df(df_subset, self.ureg)
                # FIXME: mean and mean_sd converted from l/min to meter^3/kilogram/second
                dset.unit_conversion("mean", factor=1000 * 60 * 1000 / Q(73, "kg"))
                dsets[self.to_sid(sid)] = dset
            elif sid == "Figg1995_Fig1":
                for group in df.group.unique():
                    dset = DataSet.from_df(df, self.ureg)
                    dsets[self.to_sid(f"{sid}_{group}")] = dset[df.group == group]
            else:
                dset = DataSet.from_df(df, self.ureg)
                dsets[self.to_sid(sid)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Scan]:
        return {
            "scan_cirrhosis": Scan(
                simulation=self.normal_timecourse(),
                dimensions=[
                    Dimension(
                        "dim_cirrhosis",
                        values={
                            **self.cirrhosis_scan_changes(num=50),
                        },
                    ),
                ],
            )
        }

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()
        return self.fig_icgpk_vs_cirrhosis()

    def fig_icgpk_vs_cirrhosis(self) -> dict[str, Figure]:
        """Map the CTP classes on the degree of cirrhosis.

        Every CTP class of the data is mapped on the degree of cirrhosis of the
        scan with the closest pharmacokinetic parameter. The mean of the
        mappings of the parameters is the degree of cirrhosis of a CTP class.

        Returns:
            The figures by id, `fig_<parameter>` of the mappings of the
            parameters and `fig_mapping` of the CTP classes.
        """
        cmap = scan_color(9)
        figures = {}
        markers = ["s", "^", "D", "o"]

        figg1995_groups = [
            "Controls",
            "Mild_cirrhosis",
            "Moderate_cirrhosis",
            "Severe_cirrhosis",
        ]

        # the data, its label and the unit of every parameter
        data_map: dict[str, tuple[pd.DataFrame, str, str]] = {
            "thalf_kel": (
                self._datasets["Moeller2018_Tab2_thalf"],
                "Moeller 2018",
                "min",
            ),
            "kel": (
                pd.concat(
                    [
                        self._datasets["Moeller2018_Tab2_kel"],
                        self._datasets["Herold2001_Fig1mean"],
                    ]
                ),
                "Moeller 2018 & Herold 2001",
                "1/min",
            ),
            "icgr15": (
                self._datasets["Moeller2018_Tab2_retention_ratio"],
                "Moeller 2018",
                "%",
            ),
            "cl_kel": (
                pd.concat(
                    [
                        self._datasets["Figg1995_Fig1mean"],
                        self._datasets["Moeller1998_Tab2"],
                    ]
                ),
                "Figg 1995 & Moeller 1998",
                "ml/min/kg",
            ),
        }

        # f_cirrhosis of every CTP score or class (A, B, C, ...) by parameter
        ctp_vs_f_dict: dict[str, dict[str | float, float]] = {}

        for pk_key, (dset, data_label, unit) in data_map.items():
            f, ax = base_figure()

            # cirrhosis scan, the median and the band of the draws in a run
            # over them; the clearance per body weight of 75 kg
            df = self.pk_dfs["scan_cirrhosis"]
            x = np.asarray(self.sampling, dtype=float)
            per_kg = 75.0 if pk_key == "cl_kel" else 1.0
            y, lower, upper = self.pk_band(
                df, pk_key, "ml/min" if pk_key == "cl_kel" else unit
            )
            y = y / per_kg
            if lower is not None and upper is not None:
                ax.fill_between(
                    x,
                    lower / per_kg,
                    upper / per_kg,
                    color="black",
                    alpha=self.band_alpha,
                    linewidth=0,
                )

            # plot cirrhosis scan
            ax.scatter(
                x,
                y,
                marker="o",
                c=x,
                label="Scan: cirrhosis",
                cmap=cmap,
            )
            ax.plot(x, y, marker="o", linestyle="-", ms=0, zorder=0)
            x_stars = list(self.cirrhosis_map.values())
            y_stars = np.interp(x_stars, x, y)
            for k, (key, color) in enumerate(self.cirrhosis_colors.items()):
                ax.plot(
                    x_stars[k],
                    y_stars[k],
                    marker="*",
                    linestyle="",
                    ms=20,
                    color=color,
                    label=key,
                    zorder=100,
                )

            # ctp mapping
            x_values: dict[str | float, float] = {}

            # vectors
            pk_vec = dset["mean"]
            ctps = dset["ctp"]
            sds = dset["mean_sd"]
            counts = dset["count"]

            for ctp, value, sd, count in zip(ctps, pk_vec, sds, counts, strict=True):
                label = data_label if value == pk_vec.iloc[0] else None

                x_index = np.nanargmin(np.abs(y - value))
                x_values[ctp] = x[x_index]
                ax.errorbar(
                    x[x_index],
                    value,
                    yerr=sd,
                    capsize=3,
                    capthick=2,
                    color="black",
                    marker="s",
                    linestyle="",
                    linewidth=2,
                    label=label,
                )

                # annotation
                if pk_key in ["icgr15", "thalf_kel"]:
                    xy = (x[x_index] - 0.09, value + sd + 0.8)
                elif pk_key == "kel":
                    xy = (x[x_index] + 0.02, value + sd)
                elif ctp in [5.0, 6.0, "A", "B"]:
                    xy = (x[x_index], value + sd + 0.6)
                elif ctp == "C":
                    xy = (x[x_index] + 0.005, value - sd + 0.5)
                elif ctp == 8.5:
                    xy = (x[x_index], value + sd + 0.7)
                else:
                    xy = (x[x_index], value + sd)
                ax.annotate(f"CTP: {ctp} (n={count})", xy=xy)

                # vertical lines
                ax.plot(
                    [x[x_index], x[x_index]],
                    [0, value],
                    color="dimgrey",
                    linestyle="--",
                    linewidth=2,
                    zorder=0,
                )

            # individual data points
            if pk_key == "cl_kel":
                x_vec = np.array(list(x_values.values()))  # for subscriptability
                for k, group in enumerate(figg1995_groups):
                    individual_dset = self._datasets[
                        self.to_sid(f"Figg1995_Fig1_{group}")
                    ]
                    for y_individual in individual_dset.get_quantity("y").magnitude:
                        ax.plot(
                            x_vec[k],
                            y_individual,
                            linestyle="",
                            color="black",
                            marker="s",
                            markersize="3",
                        )

            ctp_vs_f_dict[pk_key] = x_values

            ax.set_xlabel("Degree of cirrhosis [-]", fontdict=self.font)
            ax.set_ylabel(f"{self.pk_labels[pk_key]} [{unit}]", fontdict=self.font)
            ax.set_ylim(bottom=0)
            if pk_key == "kel":
                ax.legend(fontsize=self.legend_font_size, loc="lower left")
            else:
                ax.legend(fontsize=self.legend_font_size)

            figures[f"fig_{pk_key}"] = f

        # ctp vs f_cirrhosis plot:
        f_mapping, ax = base_figure()

        mean_dict: dict[str | float, list[float]] = {
            "Controls": [],
            "A": [],
            "B": [],
            "C": [],
        }

        for (pk_key, value_dict), marker in zip(
            ctp_vs_f_dict.items(), markers, strict=True
        ):
            # only the CTP classes, not the CTP scores
            for key in list(value_dict):
                if isinstance(key, float):
                    del value_dict[key]

            for key, value in value_dict.items():
                mean_dict[key].append(value)

            ax.plot(
                list(value_dict.values()),
                list(value_dict.keys()),
                color="black",
                linestyle="",
                marker=marker,
                label=self.pk_labels[pk_key],
            )
            ax.set_xlabel("Degree of cirrhosis [-]", fontdict=self.font)
            ax.set_ylabel("CTP-class", fontdict=self.font)

        means = [np.mean(values) for values in mean_dict.values()]
        sds = [np.std(values) for values in mean_dict.values()]
        logger.info(
            "f_cirrhosis of the CTP classes (mean, SD): %s",
            {
                key: (mean, sd)
                for key, mean, sd in zip(mean_dict, means, sds, strict=True)
            },
        )
        ax.errorbar(
            means[0],
            0,
            xerr=sds[0],
            color="black",
            label=r"Mean$\pm$SD",
            marker="*",
            ms=20,
        )
        ax.plot(means, list(mean_dict.keys()), color="black", linestyle="-", zorder=0)

        # position of the label of every stage of cirrhosis
        label_offsets = {
            "Control": (-0.09, 0.1),
            "Mild cirrhosis": (-0.2, 0.1),
            "Moderate cirrhosis": (-0.31, 0.1),
        }
        for k, key in enumerate(self.cirrhosis_map):
            if k != 0:
                ax.errorbar(
                    means[k],
                    k,
                    xerr=sds[k],
                    color=self.cirrhosis_colors[key],
                    marker="*",
                    zorder=0,
                    ms=20,
                )

            # vertical lines
            ax.plot(
                [means[k], means[k]],
                [-1, k],
                color=self.cirrhosis_colors[key],
                linestyle="--",
                linewidth=2,
                zorder=0,
            )
            dx, dy = label_offsets.get(key, (-0.25, 0.08))
            ax.annotate(key, xy=(means[k] + dx, k + dy), fontsize=12)
        ax.set_ylim(bottom=-0.2)
        ax.legend(fontsize=self.legend_font_size)

        figures["fig_mapping"] = f_mapping

        return figures
