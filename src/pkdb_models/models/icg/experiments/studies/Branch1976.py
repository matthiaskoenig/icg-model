"""Simulation experiment of the study Branch1976."""

from typing import ClassVar, override

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


class Branch1976(ICGSimulationExperiment):
    """Simulation experiment of Branch1976.

    Blood clearance of indocyanine green against the serum bilirubin (Table I)
    in 6 normal subjects and 20 patients with stable chronic liver disease
    (cirrhosis of several causes, chronic active hepatitis, scleroderma,
    nodular transformation; several with ascites or a surgical portacaval
    shunt) after a rapid intravenous bolus of 0.5 mg/kg. The clearance is the
    dose over the area of a biexponential fit (venous samples every 3 min for
    21 min in the normals, up to 90 min in the patients), converted to a whole
    blood clearance with the hematocrit: the model's blood clearance
    `cl_blood` (the plasma clearance over 1 - HCT) of a simulation of 60 min
    is compared. Simulated as a scan of the bilirubin for every stage of
    cirrhosis; bilirubin is one of many determinants of the clearance in these
    patients, so the figure is a qualitative validation, no fit mapping.

    Patient 19 of Table I (bilirubin 8.6 mg/dl, albumin 2.8 g/dl) has a
    clearance of 678 ml/min, which is no point of Fig 1B and does not give the
    mean of the patients with albumin < 3 g/dl in Table II (182 +- 55 ml/min,
    with 678 the mean is 248 ml/min): a probable misprint, not shown.
    """

    #: the patients of Table I which are not used, see the docstring
    excluded: ClassVar[set[str]] = {"P19"}
    #: the end of the simulations [min]
    end: ClassVar[float] = 60.0

    @override
    def datasets(self) -> dict[str, DataSet]:
        tab_id = "Tab1"
        df = self.load_dataframe(tab_id)
        df = df[~df.individual.isin(self.excluded)]
        # the clearance and the bilirubin of every subject in one row
        columns = {
            "cl": (df.substance == "indocyanine green")
            & (df.measurement_type == "clearance"),
            "bilirubin": (df.substance == "bilirubin")
            & (df.measurement_type == "concentration_ss"),
        }
        frames = [
            df[mask][["individual", "group", "value", "unit"]]
            .rename(columns={"value": key, "unit": f"{key}_unit"})
            .set_index(["individual", "group"])
            for key, mask in columns.items()
        ]
        df_subjects = pd.concat(frames, axis=1, join="inner").reset_index()
        return {
            f"{tab_id}_{group}": DataSet.from_df(
                df_subjects[df_subjects.group == group], self.ureg
            )
            for group in ["control", "patients"]
        }

    @override
    def simulations(self) -> dict[str, Scan]:
        tcs = self.cirrhosis_timecourses(end=self.end)
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
        return {"fig_cl_vs_bil": self._figure_clearance_vs_bilirubin()}

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

    def _figure_clearance_vs_bilirubin(self) -> FigureMPL:
        """Figure of the blood clearance of indocyanine green against the bilirubin."""
        f, ax = base_figure()
        ax.set_xlabel("Plasma bilirubin [mmole/l]", fontdict=self.font)
        ax.set_ylabel(
            f"{self.pk_labels['cl_blood']} [{self.pk_units['cl_blood']}]",
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
                key="cl_blood",
                color=self.cirrhosis_colors[key],
                label=key,
            )

        # data
        for group, color, label in [
            ("control", "black", "normal"),
            ("patients", "purple", "chronic liver disease"),
        ]:
            dset = self._datasets[f"Tab1_{group}"]
            x = dset.get_quantity("bilirubin")
            y = dset.get_quantity("cl")
            ax.plot(
                self.bilirubin_mM(x.magnitude, str(x.units)).magnitude,
                y.to(self.pk_units["cl_blood"]).magnitude,
                marker="s",
                linestyle="",
                color=color,
                label=f"Branch 1976 {label}",
            )

        ax.set_xscale("log")
        # room for the legend above the data
        ax.set_ylim(bottom=0, top=3600)
        ax.legend(fontsize=self.legend_font_size - 2, loc="upper left")
        return f


if __name__ == "__main__":
    MODEL.run_experiments(Branch1976, output_dir=Branch1976.__name__)
