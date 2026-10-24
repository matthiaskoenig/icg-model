"""Simulation experiment of the study Ohwada2006."""

from typing import ClassVar, override

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from pkpdsim.scan import start_values
from sbmlsim import Q
from sbmlsim.data import Data, DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan, Simulation

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure
from pkdb_models.models.icg.icg_pk import PK_SELECTIONS, summarize
from pkdb_models.models.icg.nhanes.sampling import sampled_changes


class Ohwada2006(ICGSimulationExperiment):
    """Simulation experiment of Ohwada2006.

    Estimated remnant vs measured elimination rate of indocyanine green after
    hepatectomy, 20 mg intravenous, pulse dye densitometry (with the data of
    Okochi2002 and Sunagawa2021). Ohwada2006: 75 patients with hepatocellular
    carcinoma (63 with hepatitis C, 72 Child-Pugh A), the rate measured just
    after the hepatectomy (Fig 1); Okochi2002 (Fig 2) and Sunagawa2021 (Fig 2A)
    on the first postoperative day. The estimated remnant rate is the
    preoperative rate times 1 - the resection rate. The rate of the model is
    the regression of 2 - 15 min (`kel`, the densitometer fits the first
    minutes after the mixing); the model has no regeneration or other
    postoperative change.
    """

    #: the study, its figure and its label
    studies: ClassVar[dict[str, tuple[str, str]]] = {
        "Ohwada2006": ("Fig1", "Ohwada 2006 (after resection)"),
        "Okochi2002": ("Fig2", "Okochi 2002 (POD1)"),
        "Sunagawa2021": ("Fig2", "Sunagawa 2021 (POD1)"),
    }
    dose: ClassVar[float] = 20.0  # [mg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for study, (fig_id, _) in self.studies.items():
            df = self.load_dataframe(fig_id, study=study)
            dsets[self.to_sid(study)] = DataSet.from_df(df, self.ureg)
        return dsets

    @override
    def data(self) -> dict[str, Data]:
        # the scans record only what the pharmacokinetic parameters and the
        # figures (the sampled resection rates) read
        self.add_selections_data(selections=[*PK_SELECTIONS, "resection_rate"])
        return {}

    @override
    def simulations(self) -> dict[str, Simulation | Scan]:
        tcscans: dict[str, Simulation | Scan] = {}
        tcs = self.cirrhosis_timecourses(dose=self.dose)
        for key in self.cirrhosis_map:
            tc = tcs[key]
            # scan of liver volume
            tcscans[self.to_sid(f"scan_rr_{key}")] = Scan(
                simulation=tc,
                dimensions=[
                    Dimension(
                        "dim_rr",
                        values={
                            "resection_rate": Q(self.resection_rates, "dimensionless"),
                        },
                    ),
                ],
            )

            # reference simulation pre-op
            tcscans[self.to_sid(f"scan_reference_{key}")] = tc

            # bloodflow variation on top of resection
            tcscans[self.to_sid(f"sampled_rr_bf_{key}")] = Scan(
                simulation=tc,
                dimensions=[
                    Dimension(
                        "dim_rr_bfvariation",
                        values=sampled_changes(
                            ["f_bloodflow", "resection_rate"],
                            size=400,
                            random_state=42,
                            hepatectomy=True,
                        )
                        # the samples run with the fitted values, not over the draws
                        | self.fitted_values(400),
                    ),
                ],
            )

        return tcscans

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk()
        return {
            **self._figure_icg_resection(),
            **self._figure_icg_resection_bfvariation(),
        }

    def _kel(self, sim_key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get the elimination rate of every point of a scan [1/min].

        Returns:
            The rate (the median of the draws in a run over them), and the
            lower and the upper quantile (the rate without draws).
        """
        df: pd.DataFrame = self.pk_dfs[self.to_sid(sim_key)]
        unit = df["kel_unit"].iloc[0]
        run = self.uncertainty
        if run is None:
            kel = Q(df["kel"].to_numpy(dtype=float), unit).to("1/min").magnitude
            return kel, kel, kel
        summary = summarize(df, "kel", run.quantiles)
        median, lower, upper = (
            Q(summary[column].to_numpy(), unit).to("1/min").magnitude
            for column in ["median", "lower", "upper"]
        )
        return median, lower, upper

    def _figure_template(self) -> tuple[FigureMPL, Axes]:
        """Get the figure with the data of the studies."""
        f, ax = base_figure()
        ax.set_xlabel("Estimated remnant ICG-kel [1/min]", fontdict=self.font)
        ax.set_ylabel("Postoperative ICG-kel [1/min]", fontdict=self.font)

        # data
        markers = ["s", "*", "D"]
        for (study_id, (_, label)), marker in zip(
            self.studies.items(), markers, strict=True
        ):
            x = self._datasets[self.to_sid(study_id)].get_quantity("x")
            y = self._datasets[self.to_sid(study_id)].get_quantity("y")
            ax.plot(
                x.magnitude,
                y.magnitude,
                marker=marker,
                ms=3.5,
                linestyle="",
                color="white",
                markeredgecolor="black",
                label=label,
                alpha=1,
                zorder=2,
            )
        return f, ax

    def _legend(self, ax: Axes) -> None:
        """Draw the legend right of the axes, the simulations first."""
        n_data = len(self.studies)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles[n_data:] + handles[:n_data],
            labels[n_data:] + labels[:n_data],
            fontsize=self.legend_font_size,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
        )

    def _figure_icg_resection(self) -> dict[str, FigureMPL]:
        """Estimated remnant kel vs measured postoperative kel."""
        f, ax = self._figure_template()

        for key in self.cirrhosis_map:
            kel, lower, upper = self._kel(f"scan_rr_{key}")
            # the estimated remnant kel with the preoperative kel (index 0)
            erkel = (1 - self.resection_rates) * kel[0]
            if self.uncertainty is not None:
                ax.fill_between(
                    erkel, lower, upper, color=self.cirrhosis_colors[key], alpha=0.2
                )
            ax.plot(
                erkel,
                kel,
                marker="o",
                ms=3.5,
                color=self.cirrhosis_colors[key],
                label=f"Sim: {key}",
            )

        self._legend(ax)
        return {"fig1": f}

    def _figure_icg_resection_bfvariation(self) -> dict[str, FigureMPL]:
        """Estimated remnant kel vs measured postoperative kel.

        With the sampling of the blood flow after the resection.
        """
        f, ax = self._figure_template()
        alpha = {
            "Control": 0.7,
            "Mild cirrhosis": 0.7,
            "Moderate cirrhosis": 0.2,
            "Severe cirrhosis": 0.1,
        }

        for key in self.cirrhosis_map:
            sim_key = f"sampled_rr_bf_{key}"
            xres = self.results[self.to_sid(f"task_{sim_key}")]
            resection_rates = (
                start_values(xres, "resection_rate")
                .transpose("dim_rr_bfvariation", ...)
                .to_numpy()
                .reshape(xres.ds.sizes["dim_rr_bfvariation"], -1)[:, 0]
            )

            # the preoperative kel and the kel of the samples
            kel_ref, _, _ = self._kel(f"scan_reference_{key}")
            kel, _, _ = self._kel(sim_key)

            # estimated remnant kel
            erkel = (1 - resection_rates) * kel_ref[0]

            ax.plot(
                erkel,
                kel,
                marker="o",
                ms=3.5,
                color=self.cirrhosis_colors[key],
                label=f"Sim: {key}",
                linestyle="",
                alpha=alpha[key],
                zorder=1,
            )

        self._legend(ax)
        return {"fig1_bfvariation": f}
