"""Simulation experiment of the study Kawasaki1988."""

from typing import ClassVar, override

import numpy as np
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure
from pkdb_models.models.icg.icg_pk import summarize


class Kawasaki1988(ICGSimulationExperiment):
    """Simulation experiment of Kawasaki1988.

    Dependency of the blood clearance of indocyanine green on the plasma
    bilirubin in healthy subjects (5 men) and in cirrhosis (7 patients) after
    an intravenous bolus of 0.5 mg/kg. The clearance of the paper is the blood
    clearance `CL_B,ICG = dose/AUC * 1/(1 - Hct)` with the area under the
    plasma curve to infinity of a biexponential fit of the samples from 2 to
    80 min; the individual hematocrit is not reported. The model's blood
    clearance is `cl_blood` (the dose over the area to infinity of an 80 min
    simulation over `1 - HCT` of the model) per body weight. The body weight
    is not reported, the model's 75 kg are used. The bilirubin of the subjects
    (0.3 - 2.2 mg/dl, Table I) is the x value of the data; the scan of the
    bilirubin is an extrapolation beyond it.
    """

    #: the groups of the data with their color and label; the patient with
    #: liver fibrosis (the group EPO of the curation) is not shown
    groups: ClassVar[dict[str, tuple[str, str]]] = {
        "normal": ("black", "healthy (n=5)"),
        "cirrhosis": ("purple", "cirrhosis (n=7)"),
    }
    #: the body weight of the simulations [kg]
    bodyweight: ClassVar[float] = 75.0

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Tab3")
        df_icg = df[
            (df.substance == "indocyanine green") & (df.measurement_type == "clearance")
        ]
        # the serum bilirubin of Table I [mg/dl] of every subject
        df_bil = df[df.substance == "bilirubin"][["individual", "value"]].rename(
            columns={"value": "bilirubin"}
        )
        df_icg = df_icg.merge(df_bil, on="individual", how="left")
        for group in self.groups:
            # the blood clearance [ml/min/kg]
            dsets[self.to_sid(f"Tab3_ICGclb_{group}")] = DataSet.from_df(
                df_icg[df_icg.group == group], self.ureg
            )
        return dsets

    @override
    def simulations(self) -> dict[str, Scan]:
        # 0.5 mg/kg, the samples of the study to 80 min
        tcs = self.cirrhosis_timecourses(dose=0.5 * self.bodyweight, end=80)
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
        return {**self._figure_icgpk_vs_bilirubin()}

    def _figure_icgpk_vs_bilirubin(self) -> dict[str, FigureMPL]:
        """Blood clearance of indocyanine green per body weight vs plasma bilirubin."""
        f, ax = base_figure()
        ax.set_xlabel("Plasma bilirubin [mmole/l]", fontdict=self.font)
        ax.set_ylabel(
            f"{self.pk_labels['cl_blood']}/BW\n[ml/min/kg]", fontdict=self.font
        )
        # the bilirubin of the model
        ax.axvline(x=0.01, color="grey", linestyle="--")

        x = (
            self.bilirubin_scan_changes()["[LI__bil_ext]"]
            .to(self.unit_bilirubin)
            .magnitude
        )
        bodyweight = Q(self.bodyweight, "kg")
        for key in self.cirrhosis_map:
            df = self.pk_dfs[self.to_sid(f"scan_bilirubin_{key}")]
            unit = df["cl_blood_unit"].iloc[0]

            def convert(values: np.ndarray, unit: str = unit) -> np.ndarray:
                """Get the blood clearance per body weight in ml/min/kg."""
                return (Q(values, unit) / bodyweight).to("ml/min/kg").magnitude

            if self.uncertainty is None:
                y = convert(df["cl_blood"].to_numpy(dtype=float))
            else:
                summary = summarize(df, "cl_blood", self.uncertainty.quantiles)
                y, lower, upper = (
                    convert(summary[column].to_numpy())
                    for column in ["median", "lower", "upper"]
                )
                ax.fill_between(
                    x,
                    lower,
                    upper,
                    color=self.cirrhosis_colors[key],
                    alpha=0.2,
                    linewidth=0,
                )
            ax.plot(
                x,
                y,
                marker="o",
                linestyle="-",
                color=self.cirrhosis_colors[key],
                label=key,
            )

        # data: the bilirubin of every subject against its blood clearance
        for group, (color, label) in self.groups.items():
            dset = self._datasets[self.to_sid(f"Tab3_ICGclb_{group}")]
            x_data = self.bilirubin_mM(dset["bilirubin"].to_numpy(), "mg/dl")
            y_data = dset.get_quantity("value").to("ml/min/kg")
            ax.plot(
                x_data.to(self.unit_bilirubin).magnitude,
                y_data.magnitude,
                marker="s",
                linestyle="",
                color=color,
                label=f"Kawasaki 1988 {label}",
            )

        ax.set_xscale("log")
        # the legend above the data
        ax.set_ylim(bottom=0, top=26)
        ax.legend(fontsize=self.legend_font_size - 1, loc="upper right")

        return {"fig_cl_blood_perkg": f}
