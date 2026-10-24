"""Simulation experiment of the study Sakka2006."""

from typing import override

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pkpdsim.scan import scan_dimensions, start_values
from sbmlsim import Q
from sbmlsim.data import Data, DataSet
from sbmlsim.model import RoadrunnerSBMLModel
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.icg_pk import PK_SELECTIONS
from pkdb_models.models.icg.nhanes.sampling import sampled_changes


class Sakka2006(ICGSimulationExperiment):
    """Simulation experiment of Sakka2006.

    Plasma disappearance rate vs blood clearance of indocyanine green in 209
    critically ill, ventilated patients (sepsis, ARDS, head injury, shock),
    0.3 mg/kg into a central vein, measured with a fiber-optic catheter in the
    aorta (COLD system, Fig 1). The blood clearance of the system is the total
    blood volume (dose over the blood concentration at the time 0, the back
    extrapolation) times the disappearance rate, per body surface area
    [ml/min/m^2]; normal 500 - 750 ml/min/m^2 and a PDR above 18 %/min. The
    model gives it as `cl_kel` / (1 - HCT) over the body surface area of the
    model (`BSA`, from the body weight and the height of 170 cm). The patients
    are not described by the model: a qualitative comparison of the relation
    of the PDR and the clearance, with a sampling of the body weight and of
    the amount of the transport protein.
    """

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for df in self.load_dataframes_by_substance("Fig1").values():
            # the clearance and the plasma disappearance rate, each in its
            # unit: a data set holds the rows of one measurement (sbmlsim reads a
            # column in one unit)
            for measurement in df.measurement.unique():
                dsets[f"Fig1_ICG{measurement}"] = DataSet.from_df(
                    df[df.measurement == measurement], self.ureg
                )
        return dsets

    @override
    def data(self) -> dict[str, Data]:
        # the scans record only what the pharmacokinetic parameters and the
        # clearance per body surface area read
        self.add_selections_data(selections=[*PK_SELECTIONS, "BSA"])
        return {}

    @override
    def simulations(self) -> dict[str, Scan]:
        n_samples = 20
        dose_per_bw = Q(0.3, "mg/kg")
        bodyweights = Q(np.linspace(45, 90, num=n_samples), "kg")

        # the fitted uptake of the model
        model = self._models["model"]
        if not isinstance(model, RoadrunnerSBMLModel) or model.r is None:
            raise TypeError(f"The model of '{self.sid}' is not loaded: '{model}'.")
        vmax_ref = Q(model.r["LI__ICGIM_Vmax"], model.uinfo["LI__ICGIM_Vmax"])
        vmax_scan = Q(
            np.linspace(vmax_ref.magnitude * 0.2, vmax_ref.magnitude * 4, n_samples),
            vmax_ref.units,
        )

        # sampling of variability, the dose of 0.3 mg/kg of every body weight
        changes = sampled_changes(pars=["LI__f_oatp1b3", "BW"], random_state=42)
        changes["IVDOSE_icg"] = (dose_per_bw * changes["BW"]).to("mg")
        # the samples run with the fitted values, not over the draws
        changes.update(self.fitted_values(np.size(changes["BW"].magnitude)))

        return {
            "scan_bw_protein_sampling": Scan(
                # the scan sets the dose
                simulation=self.normal_timecourse(dose=0),
                dimensions=[Dimension("dim_sampling", values=changes)],
            ),
            "scan_bw": Scan(
                # the scan sets the dose
                simulation=self.normal_timecourse(dose=0),
                dimensions=[
                    Dimension(
                        "dim_bodyweight",
                        values={
                            "BW": bodyweights,
                            "IVDOSE_icg": (dose_per_bw * bodyweights).to("mg"),
                        },
                    ),
                ],
            ),
            "scan_protein": Scan(
                simulation=self.normal_timecourse(dose=0.3 * 75),
                dimensions=[
                    Dimension("dim_protein", values={"LI__ICGIM_Vmax": vmax_scan})
                ],
            ),
        }

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk()
        return {**self._figure_icgpdr()}

    def _start_values(self, sim_key: str, df: pd.DataFrame, key: str) -> np.ndarray:
        """Get a variable at the start of the simulation of every row of the parameters.

        Args:
            sim_key: the key of the scan.
            df: the parameters of the scan, a row per point (and draw).
            key: the variable.
        """
        res = self.results[f"task_{sim_key}"]
        across = None if self.uncertainty is None else self.uncertainty.id
        if across is not None and across not in res.dims:
            # a scan which runs with the fitted values, not over the draws
            across = None
        dims = scan_dimensions(res, across)
        start = start_values(res, key)
        values = []
        for _, row in df.iterrows():
            position = {dim: int(row["point"]) for dim in dims}
            if across is not None:
                position[across] = int(row["draw"])
            values.append(start.isel(position).item())
        return Q(np.asarray(values), res.units[key]).magnitude

    def _blood_clearance_per_bsa(self, sim_key: str, df: pd.DataFrame) -> np.ndarray:
        """Get the blood clearance of the COLD system per body surface area.

        The total blood volume times the disappearance rate,
        `cl_kel` / (1 - HCT), over the body surface area [ml/min/m^2].
        """
        cl_kel = Q(df["cl_kel"].to_numpy(dtype=float), df["cl_kel_unit"].iloc[0])
        hct = self._start_values(sim_key, df, "HCT")
        bsa = Q(self._start_values(sim_key, df, "BSA"), "m**2")
        return (cl_kel / (1 - hct) / bsa).to("ml/min/m**2").magnitude

    def _figure_icgpdr(self) -> dict[str, FigureMPL]:
        """Plasma disappearance rate vs blood clearance per body surface area.

        In a run over the draws of the fitted parameters every point is the
        median of its draws.
        """
        f, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlabel("ICG-blood clearance [ml/min/m²]", fontdict=self.font)
        ax.set_ylabel(f"{self.pk_labels['icgpdr']} [%/min]", fontdict=self.font)
        ax.tick_params(axis="both", labelsize=self.tick_font_size)

        # data
        y = self._datasets["Fig1_ICGpdr"].get_quantity("value")
        x = self._datasets["Fig1_ICGclearance"].get_quantity("value")
        ax.plot(
            x.to("ml/min/m**2").magnitude,
            y.to("percent/min").magnitude,
            marker="s",
            linestyle="",
            color="black",
            label="Sakka 2006",
            ms=3,
        )

        # label and color of the scans
        styles = {
            "scan_bw_protein_sampling": ("Sim: body weight & protein", "orange"),
            "scan_protein": ("Sim: protein scan", "green"),
            "scan_bw": ("Sim: body weight scan", "blue"),
        }
        for key, (label, color) in styles.items():
            df = self.pk_dfs[key].copy()
            df["clb_bsa"] = self._blood_clearance_per_bsa(key, df)
            df["pdr"] = (
                Q(df["icgpdr"].to_numpy(dtype=float), df["icgpdr_unit"].iloc[0])
                .to("percent/min")
                .magnitude
            )
            points = df.groupby("point")[["clb_bsa", "pdr"]].median()
            ax.plot(
                points["clb_bsa"].to_numpy(),
                points["pdr"].to_numpy(),
                marker="o",
                linestyle="",
                color=color,
                ms=4,
                label=label,
            )

        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.legend(
            fontsize=self.legend_font_size, loc="upper left", bbox_to_anchor=(1.02, 1.0)
        )

        return {"fig_icgpdr_vs_clearance": f}
