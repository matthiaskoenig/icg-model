"""Simulation experiment of the study DOnofrio2014."""

from typing import override

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure
from pkdb_models.models.icg.icg_pk import summarize


class DOnofrio2014(ICGSimulationExperiment):
    """Simulation experiment of DOnofrio2014.

    Retention of indocyanine green after 15 min against the serum bilirubin in
    31 candidates for liver surgery (liver tumours: hepatocellular carcinoma,
    cholangiocarcinoma, one metastasis; 7 with chronic hepatitis B, 1 with
    hepatitis C; bilirubin < 3 mg/dl), fasted from midnight, after an
    intravenous bolus of 0.5 mg/kg measured by pulse dye densitometry (LiMON).
    The device computes the retention from the plasma disappearance rate,
    R15 = exp(-15 min * PDR), so the PDR of Table 1 is the same measurement and
    is not shown again; the model's R15 is relative to the back extrapolation of
    the log-linear regression (`icg_pk`). Simulated as a scan of the bilirubin
    for every stage of cirrhosis: the patients are mostly not cirrhotic, the
    figure shows that the retention does not depend on the bilirubin below
    3 mg/dl. A validation, no fit mapping.
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        tab_id = "Tab1"
        df = self.load_dataframe(tab_id)
        # the retention and the bilirubin of every patient in one row
        columns = {
            "r15": (df.substance == "indocyanine green")
            & (df.measurement_type == "retention ratio"),
            "bilirubin": (df.substance == "bilirubin")
            & (df.measurement_type == "concentration_ss"),
        }
        frames = [
            df[mask][["individual", "value", "unit"]]
            .rename(columns={"value": key, "unit": f"{key}_unit"})
            .set_index("individual")
            for key, mask in columns.items()
        ]
        df_patients = pd.concat(frames, axis=1, join="inner").reset_index()
        return {f"{tab_id}_R15_bilirubin": DataSet.from_df(df_patients, self.ureg)}

    @override
    def simulations(self) -> dict[str, Scan]:
        tcs = self.cirrhosis_timecourses()
        return {
            self.to_sid(f"scan_bilirubin_{key}"): Scan(
                simulation=tcs[key],
                dimensions=[
                    Dimension("dim_bil", values=self.bilirubin_scan_changes()),
                ],
            )
            for key in self.cirrhosis_map
        }

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk()
        return {"fig_icgr15": self._figure_r15_vs_bilirubin()}

    def _plot_scan(
        self, ax: Axes, sim_key: str, key: str, color: str, label: str
    ) -> None:
        """Plot a parameter of a scan of the bilirubin, with the band of the draws."""
        x = self.bilirubin_scan_changes()["[LI__bil_ext]"].to(self.unit_bilirubin)
        df = self.pk_dfs[sim_key]
        unit_pk = df[f"{key}_unit"].iloc[0]

        def convert(values: np.ndarray) -> np.ndarray:
            """Get the values in the unit of the figure."""
            return Q(values, unit_pk).to(self.pk_units[key]).magnitude

        if self.uncertainty is None:
            y = convert(df[key].to_numpy(dtype=float))
        else:
            summary = summarize(df, key, self.uncertainty.quantiles)
            y, lower, upper = (
                convert(summary[column].to_numpy())
                for column in ["median", "lower", "upper"]
            )
            ax.fill_between(
                x.magnitude, lower, upper, color=color, alpha=0.2, linewidth=0
            )
        ax.plot(x.magnitude, y, marker="o", linestyle="-", color=color, label=label)

    def _figure_r15_vs_bilirubin(self) -> FigureMPL:
        """Figure of the retention after 15 min against the plasma bilirubin."""
        f, ax = base_figure()
        ax.set_xlabel("Plasma bilirubin [mmole/l]", fontdict=self.font)
        ax.set_ylabel(
            f"{self.pk_labels['icgr15']} [{self.pk_units['icgr15']}]",
            fontdict=self.font,
        )
        # the upper limit of normal bilirubin, 1.2 mg/dl
        ax.axvline(
            x=self.bilirubin_mM(1.2, "mg/dl").magnitude, color="grey", linestyle="--"
        )

        for key in self.cirrhosis_map:
            self._plot_scan(
                ax,
                sim_key=self.to_sid(f"scan_bilirubin_{key}"),
                key="icgr15",
                color=self.cirrhosis_colors[key],
                label=key,
            )

        # data
        dset = self._datasets["Tab1_R15_bilirubin"]
        x = dset.get_quantity("bilirubin")
        y = dset.get_quantity("r15")
        ax.plot(
            self.bilirubin_mM(x.magnitude, str(x.units)).magnitude,
            y.to(self.pk_units["icgr15"]).magnitude,
            marker="s",
            linestyle="",
            markeredgecolor="black",
            color="white",
            label="D'Onofrio 2014 liver tumours",
        )

        ax.set_xscale("log")
        ax.legend(fontsize=self.legend_font_size)
        return f


if __name__ == "__main__":
    MODEL.run_experiments(DOnofrio2014, output_dir=DOnofrio2014.__name__)
