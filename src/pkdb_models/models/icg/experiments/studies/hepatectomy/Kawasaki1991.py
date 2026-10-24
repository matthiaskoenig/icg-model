"""Simulation experiment of the study Kawasaki1991."""

from typing import override

import numpy as np
from sbmlsim.data import DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure
from pkdb_models.models.icg.nhanes.sampling import postop_bloodflow_distribution


class Kawasaki1991(ICGSimulationExperiment):
    """Simulation experiment of Kawasaki1991.

    Change of the portal blood flow (Doppler ultrasound) by partial
    hepatectomy, the data of the distribution of the hepatic blood flow after
    the resection (`f_bloodflow`) which is sampled in the simulations of
    hepatectomy (`nhanes.sampling.postop_bloodflow_distribution`). The data
    are the ratio after/before of the 17 patients measured intraoperatively,
    before and right after the resection (Fig 2a; 15.3 +- 4.3 to 14.3 +- 5.4
    ml/min/kg, mean +- SD); the postoperative measurements (10 - 12 days, n=56,
    -9 %) are group means only. The portal flow stands in for the hepatic
    blood flow (`f_bloodflow` scales the hepatic and the gastrointestinal
    flow), an assumption. The experiment has no simulations.
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        df = self.load_dataframe("Fig2")
        return {"bf_change": DataSet.from_df(df, self.ureg)}

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        return {**self._figure_histogram()}

    def _figure_histogram(self) -> dict[str, FigureMPL]:
        """Distribution of the change of the portal blood flow by the resection."""
        f, ax = base_figure()
        ax.set_xlabel("Intraoperative PBF change [%]", fontdict=self.font)
        ax.set_ylabel("Density [-]", fontdict=self.font)

        dist, _, pars = postop_bloodflow_distribution()
        left_bound = pars["loc"]
        right_bound = pars["loc"] + pars["scale"]
        x = np.linspace(left_bound - 0.001, right_bound, 1000)
        y = dist.pdf(x)
        ax.plot(x, y, color="red", label="Uniform PDF")

        n_samples = 10000
        samples = dist.rvs(random_state=42, size=n_samples)
        ax.hist(
            x=samples,
            bins=50,
            label="Samples",
            density=True,
            color="darkgrey",
            alpha=1,
        )

        data = self._datasets["bf_change"]["change"]
        ax.hist(x=data, bins=15, label="Kawasaki 1991 (n=17)", density=True, alpha=0.8)
        ax.set_xticks(ticks=[0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6])
        ax.set_xticklabels(["-60", "-40", "-20", "0", "+20", "+40", "+60"])

        # the legend above the histograms
        ax.set_ylim(top=2.6)
        ax.legend(fontsize=self.legend_font_size, loc="upper center")

        return {"fig_hist": f}
