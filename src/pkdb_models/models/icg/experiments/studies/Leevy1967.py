"""Simulation experiment of the study Leevy1967."""

from typing import ClassVar, override

import numpy as np
import pandas as pd
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan, Simulation

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure, scan_color
from pkdb_models.models.icg.icg_pk import summarize


class Leevy1967(ICGSimulationExperiment):
    """Simulation experiment of Leevy1967.

    Dose dependency of the plasma disappearance rate of indocyanine green in
    healthy subjects (n=10) and in liver disease (Fig 6): every subject
    received 0.5 mg/kg and one hour later 5 mg/kg, both injected in less than
    10 s; the curves were recorded by ear densitometry. The initial PDR of the
    paper is the semilogarithmic slope of the points every 2 min from the
    third minute, before the second phase which began between the seventh and
    the tenth minute (`window`, 3 - 9 min after each injection). The
    simulations give the two doses in sequence; the body weight is not
    reported, the model's 75 kg are used. The scans of the degree of
    cirrhosis (`f_tissue_loss` = `f_shunts`) are compared with the data of
    the controls and the cirrhosis patients.
    """

    #: the window of the initial PDR after an injection [min]
    window: ClassVar[tuple[float, float]] = (3.0, 9.0)
    #: the time of the second injection [min]
    t_second: ClassVar[float] = 60.0
    #: the body weight of the simulations [kg]
    bodyweight: ClassVar[float] = 75.0

    #: the PDR of the first (0.5 mg/kg) and of the second injection (5 mg/kg)
    #: by simulation key
    pk_first: dict[str, pd.DataFrame]
    pk_second: dict[str, pd.DataFrame]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for df in self.load_dataframes_by_substance("Fig6").values():
            dset = DataSet.from_df(df, self.ureg)
            for group in df.group.unique():
                dsets[self.to_sid(f"Fig6_ICGpdr_{group}")] = dset[dset.group == group]
        return dsets

    def _sequence(self) -> Simulation:
        """Get the simulation of 0.5 mg/kg and 5 mg/kg one hour later."""
        return Simulation(
            time_unit="min",
            end=self.t_second + 30,
            preinit_changes={
                **self.default_changes(),
                "IVDOSE_icg": Q(0.5 * self.bodyweight, "mg"),
            },
            changes=[
                self.add_doses(
                    self.t_second,
                    {"IVDOSE_icg": Q(5.0 * self.bodyweight, "mg")},
                )
            ],
        )

    @override
    def simulations(self) -> dict[str, Scan]:
        # the stages of cirrhosis of `cirrhosis_map`
        stages = {
            "f_tissue_loss": Q(list(self.cirrhosis_map.values()), "dimensionless"),
            "f_shunts": Q(list(self.cirrhosis_map.values()), "dimensionless"),
        }
        return {
            "scan_cirrhosis": Scan(
                simulation=self._sequence(),
                dimensions=[
                    Dimension("dim_cirrhosis", values=self.cirrhosis_scan_changes()),
                ],
            ),
            "scan_cirrhosis_stages": Scan(
                simulation=self._sequence(),
                dimensions=[Dimension("dim_cirrhosis", values=stages)],
            ),
        }

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_first = self.calculate_icg_pk(window=self.window)
        self.pk_second = self.calculate_icg_pk(
            window=(self.window[0] + self.t_second, self.window[1] + self.t_second)
        )
        return {**self._figure_icgpdr()}

    def _pdr(self, df: pd.DataFrame) -> np.ndarray:
        """Get the PDR [%/min] of the points of a scan, the median of the draws."""
        if self.uncertainty is None:
            values = df["icgpdr"].to_numpy(dtype=float)
        else:
            values = summarize(df, "icgpdr", self.uncertainty.quantiles)[
                "median"
            ].to_numpy()
        return Q(values, df["icgpdr_unit"].iloc[0]).to("percent/min").magnitude

    def _figure_icgpdr(self) -> dict[str, FigureMPL]:
        """Plasma disappearance rate with 0.5 mg/kg vs with 5.0 mg/kg."""
        f, ax = base_figure()
        ax.set_xlabel("ICG-PDR [%/min] (Dose: 5.0 mg/kg)", fontdict=self.font)
        ax.set_ylabel("ICG-PDR [%/min] (Dose: 0.5 mg/kg)", fontdict=self.font)

        x = self._pdr(self.pk_second["scan_cirrhosis"])
        y = self._pdr(self.pk_first["scan_cirrhosis"])
        x_stars = self._pdr(self.pk_second["scan_cirrhosis_stages"])
        y_stars = self._pdr(self.pk_first["scan_cirrhosis_stages"])

        # simulation
        cmap = scan_color(9, reverse=True)
        ax.scatter(x, y, marker="o", c=y, label="Scan: cirrhosis", cmap=cmap, s=16)
        ax.plot(x, y, marker="o", linestyle="-", ms=0, zorder=0)
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

        # data: color, marker and label of the datasets
        data_map = {
            "Fig6_ICGpdr_Controls": ("black", "s", "controls"),
            "Fig6_ICGpdr_Cirrhosis": ("purple", "s", "cirrhosis"),
        }
        for dset_key, (color, marker, label) in data_map.items():
            dset = self._datasets[self.to_sid(dset_key)]
            ax.plot(
                dset.get_quantity("x").magnitude,
                dset.get_quantity("y").magnitude,
                marker=marker,
                linestyle="",
                color=color,
                label=f"Leevy 1967 {label}",
            )

        ax.set_xlim(left=0, right=37)
        ax.set_ylim(bottom=0, top=37)
        ax.legend(fontsize=self.legend_font_size - 2, loc="upper left")

        return {"fig_icgpdr": f}
