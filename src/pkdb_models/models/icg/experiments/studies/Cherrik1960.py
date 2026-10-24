"""Simulation experiment of the study Cherrik1960."""

from typing import ClassVar, override

import numpy as np
from matplotlib.ticker import FormatStrFormatter, NullFormatter
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure, scan_color
from pkdb_models.models.icg.icg_pk import summarize


class Cherrik1960(ICGSimulationExperiment):
    """Simulation experiment of Cherrik1960 (Cherrick GR et al., J Clin Invest 1960).

    Half-life and retention after 20 min of indocyanine green (Table VII) in 9
    patients without liver disease and in patients with liver disease (mostly
    cirrhosis), after a rapid intravenous injection of 0.5 mg/kg. Venous samples
    were taken at 3 and 5 min and every 5 min for at least 1 h; the half-life is
    from the exponential part of the curve, the retention R20 relative to the
    concentration at the time 0 (the back extrapolation). The model takes the
    log-linear regression of 3 - 20 min (`window`) for both, `thalf_kel` and
    `icgr20`. Simulated as a scan of the degree of cirrhosis (tissue loss and
    shunts) with the stages of `cirrhosis_map`. The controls also had 1.0 and
    2.0 mg/kg with an unchanged R20 (3.9 % and 4.0 % vs 3.85 % at 0.5 mg/kg,
    text). A validation, no fit mapping.
    """

    #: the window of the log-linear regression [min]
    window: ClassVar[tuple[float, float]] = (3.0, 20.0)

    colors: ClassVar[dict[str, str]] = {
        "Controls": "black",
        "Cirrhosis": "purple",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        tab_id = "Tab7"
        dsets = {}
        for df in self.load_dataframes_by_substance(tab_id).values():
            # the half-life and the retention ratio, each in its unit: a data
            # set holds the rows of one measurement (sbmlsim reads a column in
            # one unit)
            measurements = {
                measurement: DataSet.from_df(
                    df[df.measurement == measurement], self.ureg
                )
                for measurement in ["thalf", "retention ratio"]
            }
            for group in df.group.unique():
                for measurement, dset in measurements.items():
                    dsets[self.to_sid(f"{tab_id}_ICG{measurement}_{group}")] = dset[
                        dset.group == group
                    ]
        return dsets

    @override
    def simulations(self) -> dict[str, Scan]:
        return {
            "scan_cirrhosis": Scan(
                simulation=self.normal_timecourse(),
                dimensions=[
                    Dimension(
                        "dim_cirrhosis",
                        values=self.cirrhosis_scan_changes(stop=0.95),
                    ),
                ],
            ),
            "four_cirrhosis_scan": Scan(
                simulation=self.normal_timecourse(),
                dimensions=[
                    Dimension(
                        "dim_cirrhosis",
                        values={
                            "f_tissue_loss": Q(
                                list(self.cirrhosis_map.values()), "dimensionless"
                            ),
                            "f_shunts": Q(
                                list(self.cirrhosis_map.values()), "dimensionless"
                            ),
                        },
                    )
                ],
            ),
        }

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk(window=self.window)
        return {"fig_thalf_vs_r20": self._figure_thalf_vs_r20()}

    def _median(self, sim_key: str, key: str) -> np.ndarray:
        """Get a parameter of every point of a scan, the median over the draws."""
        df = self.pk_dfs[sim_key]
        if self.uncertainty is None:
            return df[key].to_numpy(dtype=float)
        return summarize(df, key, self.uncertainty.quantiles)["median"].to_numpy()

    def _figure_thalf_vs_r20(self) -> FigureMPL:
        """Figure of the retention after 20 min against the half-life."""
        # colormap of the scan
        cmap = scan_color(length=9)
        x = self._median("scan_cirrhosis", "thalf_kel")
        y = self._median("scan_cirrhosis", "icgr20")
        x_stars = self._median("four_cirrhosis_scan", "thalf_kel")
        y_stars = self._median("four_cirrhosis_scan", "icgr20")

        f, ax = base_figure()
        ax.scatter(
            x,
            y,
            marker="o",
            c=np.linspace(0, 1, len(x)),
            label="Scan: cirrhosis",
            cmap=cmap,
            s=16,
        )
        ax.plot(x, y, linestyle="-", color="grey", zorder=0)
        for k, (key, color) in enumerate(self.cirrhosis_colors.items()):
            ax.plot(
                x_stars[k],
                y_stars[k],
                marker="*",
                linestyle="",
                color=color,
                label=key,
                ms=20,
            )

        ax.set_xlabel(f"{self.pk_labels['thalf_kel']} [min]", fontdict=self.font)
        ax.set_ylabel(f"{self.pk_labels['icgr20']} [%]", fontdict=self.font)

        # data
        for group in ["Controls", "Cirrhosis"]:
            x = self._datasets[self.to_sid(f"Tab7_ICGthalf_{group}")].get_quantity(
                "value"
            )
            y = self._datasets[
                self.to_sid(f"Tab7_ICGretention ratio_{group}")
            ].get_quantity("value")
            ax.plot(
                x.magnitude,
                y.magnitude,
                marker="s",
                linestyle="",
                color=self.colors[group],
                label=f"Cherrick 1960 {group.lower()}",
            )

        # the legend in the order: scan, stages of cirrhosis, data
        handles, labels = ax.get_legend_handles_labels()
        order = [
            "Scan: cirrhosis",
            *self.cirrhosis_colors,
            "Cherrick 1960 controls",
            "Cherrick 1960 cirrhosis",
        ]
        by_label = dict(zip(labels, handles, strict=True))
        ax.legend(
            [by_label[label] for label in order],
            order,
            fontsize=self.legend_font_size,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ticks = [2, 5, 10, 20, 50, 100]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(FormatStrFormatter("%g"))
            axis.set_minor_formatter(NullFormatter())
        return f


if __name__ == "__main__":
    MODEL.run_experiments(Cherrik1960, output_dir=Cherrik1960.__name__)
