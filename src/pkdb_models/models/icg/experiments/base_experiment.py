"""Base class of the indocyanine green simulation experiments."""

from collections.abc import Mapping, Sequence
from typing import ClassVar, NamedTuple

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.scan import end_values, start_values
from sbmlsim import Q
from sbmlsim.model import RoadrunnerSBMLModel
from sbmlsim.result import ScanResult
from sbmlsim.simulation import Scan, Simulation
from sbmlsim.simulation.sampling import local, parameters_of
from sbmlsim.units import Quantity

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.icg_pk import (
    PDR_WINDOW,
    calculate_icg_parameters,
    summarize,
)


class MolecularWeights(NamedTuple):
    """Molecular weights for the conversion of the data."""

    icg: Quantity


def exclude_parameters_icg(pid: str) -> bool:
    """Filter for excluding parameter ids in sensitivity.

    For instance constants like molecular weights or conversion factors
    are not changed in sensitivity.
    """
    if pid.startswith(("Mr_", "ti_")):
        return True
    if pid in {
        "f_bloodflow",
        "f_cardiac_output",
        "f_shunts",
        "f_tissue_loss",
        "LI__f_oatp1b3",
        "resection_rate",
    }:
        return True
    # must always be 1.0
    return pid.startswith("FQlu")


class ICGSimulationExperiment(PKPDSimulationExperiment):
    """Base class of the indocyanine green simulation experiments."""

    pkpd_model = MODEL
    selections = (
        "time",
        # whole-body parameters
        "BW",
        "HCT",
        # dosing
        "IVDOSE_icg",
        "Ri_icg",
        "CLinfusion_icg",
        # extraction ratio
        "ER_icg",
        # blood flows
        "Qh",
        "Qha",
        "Qgi",
        "Qpo",
        "Qre",
        "Qlu",
        "QC",
        "FQh",
        # plasma flows, the convective flows of indocyanine green
        "Qh_plasma",
        # volumes
        "FVli",
        "Vli",
        # model settings
        "f_shunts",
        "f_tissue_loss",
        "f_bloodflow",
        "f_cardiac_output",
        "resection_rate",
        # icg concentrations/amounts
        "[Cve_icg]",
        "[Car_icg]",
        "[Chv_icg]",
        "[Cpo_icg]",
        "[LI__icg]",
        "[LI__icg_bi]",
        "[Cre_plasma_icg]",
        "[Cli_plasma_icg]",
        "Afeces_icg",
        # bilirubin concentrations
        "[LI__bil_ext]",
        # kinetic parameters
        "LI__ICGIM_Vmax",
        "LI__f_oatp1b3",
        # transport rates: the uptake, the reflux into plasma, the export into
        # bile and the biliary excretion
        "LI__ICGIM",
        "LI__ICGEX",
        "LI__ICGLI2CA",
        "LI__ICGLI2BI",
    )

    font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 18}
    #: the opacity of the bands of the draws of the fitted parameters
    band_alpha = 0.2
    tick_font_size = 15
    legend_font_size = 13

    unit_bile_rate = "µmole/min"
    unit_bilirubin = "mmole/l"
    unit_bloodflow = "l/min"
    unit_clearance = "ml/min"
    unit_dose = "mg"
    unit_er = "dimensionless"
    unit_f = "dimensionless"
    unit_icg = "µmol/l"
    unit_time = "min"

    # hepatectomy
    resection_rates: ClassVar[np.ndarray] = np.linspace(0, 0.9, num=20)

    cirrhosis_colors: ClassVar[Mapping[str, str]] = {
        "Control": "black",
        "Mild cirrhosis": "blue",
        "Moderate cirrhosis": "purple",
        "Severe cirrhosis": "red",
    }

    pk_labels: ClassVar[Mapping[str, str]] = {
        "kel": "ICG-kel",
        "lambda_z": "ICG-λz",
        "icgpdr": "ICG-PDR",
        "icgr15": "ICG-R15",
        "icgr20": "ICG-R20",
        "cl": "ICG-clearance",
        "cl_kel": "ICG-clearance",
        "clperbw": "ICG-clearance/BW",
        "cl_auc": "ICG-clearance",
        "cl_blood": "ICG-blood clearance",
        "thalf": "ICG-thalf",
        "thalf_kel": "ICG-thalf",
        "v0": "ICG-V0",
    }

    pk_units: ClassVar[Mapping[str, str]] = {
        "kel": "1/min",
        "icgpdr": "%/min",
        "icgr15": "%",
        "icgr20": "%",
        "cl": "ml/min",
        "cl_kel": "ml/min",
        "cl_auc": "ml/min",
        "cl_blood": "ml/min",
        "clperbw": "ml/min/kg",
        "thalf": "min",
        "thalf_kel": "min",
        "lambda_z": "1/min",
        "v0": "l",
    }

    pk_limits: ClassVar[dict[str, float]] = {
        "icgpdr": 45,
        "icgr15": 110,
        "cl": 1000,
        "cl_kel": 1000,
        "thalf": 125,
        "thalf_kel": 125,
    }

    #: the parameters of the liver function test which the scans show: the
    #: half-life, the clearance (`kel` * `v0`), the plasma disappearance rate
    #: and the retention after 15 min, the regression of 2 - 15 min
    pk_scan_parameters: ClassVar[tuple[str, ...]] = (
        "thalf_kel",
        "cl_kel",
        "icgpdr",
        "icgr15",
    )

    Mr: ClassVar[MolecularWeights] = MolecularWeights(icg=Q(774.96493, "g/mole"))

    #: molecular weight of bilirubin, for the conversion of serum bilirubin
    Mr_bilirubin: ClassVar[Quantity] = Q(584.6623, "g/mole")

    @classmethod
    def bilirubin_mM(cls, value: float | np.ndarray, unit: str = "mg/dl") -> Quantity:
        """Get a serum bilirubin in mmole/l, the unit of `[LI__bil_ext]`.

        Args:
            value: the bilirubin, e.g. in mg/dl (1 mg/dl = 0.0171 mmole/l).
            unit: the unit of the value, a mass or an amount concentration.

        Returns:
            The bilirubin in mmole/l.
        """
        bilirubin = Q(value, unit)
        if bilirubin.check("[mass]/[length]**3"):
            bilirubin = bilirubin / cls.Mr_bilirubin
        return bilirubin.to("mmole/l")

    @staticmethod
    def hepatic_bloodflow(
        bodyweight: float = 75.0,
        f_bloodflow: float = 1.0,
        f_cardiac_output: float = 1.0,
    ) -> Quantity:
        """Get the hepatic blood flow `Qh` of the model for a body weight.

        `Qh = BW * COBW * f_cardiac_output * FQh * f_bloodflow` with the cardiac
        output per body weight `COBW` 75 ml/min/kg and the fraction of the
        hepatic blood flow `FQh` 0.255 of the model, 1.434 l/min at 75 kg; the
        reference lines of the figures of the blood flow are drawn at it. The
        hepatic plasma flow is `Qh * (1 - HCT)` (`Qh_plasma`).

        Args:
            bodyweight: the body weight [kg].
            f_bloodflow: the factor of the hepatic blood flow.
            f_cardiac_output: the factor of the cardiac output.

        Returns:
            The hepatic blood flow in l/min.
        """
        cobw = Q(75, "ml/min/kg")
        fqh = 0.255
        return (Q(bodyweight, "kg") * cobw * f_cardiac_output * fqh * f_bloodflow).to(
            "l/min"
        )

    #: the pharmacokinetic parameters of the scans by simulation key, which the
    #: figures of an experiment calculate with `calculate_icg_pk`
    pk_dfs: dict[str, pd.DataFrame]

    def cirrhosis_timecourses(
        self, dose: float | None = 37.5, end: float = 30
    ) -> dict[str, Simulation]:
        """Get the simulations of an intravenous bolus for the stages of cirrhosis.

        Args:
            dose: intravenous dose of indocyanine green [mg], no dose for `None`
                (the simulations of a scan of the dose, which sets it).
            end: end of the simulations [min].

        Returns:
            The simulations by stage of cirrhosis, the keys of `cirrhosis_map`.
        """
        doses = {} if dose is None else {"IVDOSE_icg": Q(dose, "mg")}
        return {
            key: Simulation(
                time_unit="min",
                end=end,
                preinit_changes={
                    **self.default_changes(),
                    "f_tissue_loss": Q(value, "dimensionless"),
                    "f_shunts": Q(value, "dimensionless"),
                    **doses,
                },
            )
            for key, value in self.cirrhosis_map.items()
        }

    def normal_timecourse(self, dose: float = 37.5, end: float = 30) -> Simulation:
        """Get the simulation of an intravenous bolus of indocyanine green.

        Args:
            dose: intravenous dose of indocyanine green [mg], no dose for 0 (the
                simulation of a scan of the dose, which sets it).
            end: end of the simulation [min].

        Returns:
            The simulation.
        """
        changes: dict[str, Quantity | float] = dict(self.default_changes())
        if dose != 0:
            changes["IVDOSE_icg"] = Q(dose, "mg")
        return Simulation(time_unit="min", end=end, preinit_changes=changes)

    def cirrhosis_scan_changes(
        self, start: float = 0, stop: float = 0.9, num: int = 20
    ) -> dict[str, Quantity]:
        """Get the changes of a scan of the degree of cirrhosis.

        The tissue loss and the shunts change together. The changes are vectors,
        the values of a dimension of a scan.
        """
        return {
            "f_tissue_loss": Q(
                np.linspace(start=start, stop=stop, num=num), "dimensionless"
            ),
            "f_shunts": Q(
                np.linspace(start=start, stop=stop, num=num), "dimensionless"
            ),
        }

    def bloodflow_scan_changes(
        self, start: float = 0.1, stop: float = 2.0, num: int = 20
    ) -> dict[str, Quantity]:
        """Get the changes of a scan of the hepatic blood flow (`f_bloodflow`)."""
        return {
            "f_bloodflow": Q(
                np.linspace(start=start, stop=stop, num=num), "dimensionless"
            ),
        }

    def bilirubin_scan_changes(
        self, start: float = -3, stop: float = -0.5, num: int = 20, log: bool = True
    ) -> dict[str, Quantity]:
        """Get the changes of a scan of the bilirubin concentration.

        Scan from start to stop with num points. The log flag makes either a
        linear (False) or logarithmic (True) scan.
        """
        if log:
            bil_array = Q(np.logspace(start, stop, num=num), "mmole/l")
        else:
            bil_array = Q(np.linspace(start, stop, num=num), "mmole/l")

        return {
            "[LI__bil_ext]": bil_array,
        }

    def calculate_icg_pk(
        self,
        scans: Sequence[str] | None = None,
        window: tuple[float, float] = PDR_WINDOW,
    ) -> dict[str, pd.DataFrame]:
        """Calculate the ICG parameters of the scans, see `icg_pk`.

        In a run over the draws of the fitted parameters (`uncertainty`) the
        parameters are calculated for every draw (the column `draw`).

        Args:
            scans: keys of the simulations, all simulations by default.
            window: the window of the log-linear regression of the liver
                function test [min] (`kel`, `icgpdr`, `c0`, `icgr15`).

        Returns:
            The parameters of every scan by simulation key.
        """
        keys = scans if scans else list(self._simulations)
        across = None if self.uncertainty is None else self.uncertainty.id
        return {
            sim_key: calculate_icg_parameters(
                experiment=self,
                xres=self.results[f"task_{sim_key}"],
                across=across,
                window=window,
            )
            for sim_key in keys
        }

    def _across(self) -> str | None:
        """Get the dimension of the draws of the fitted parameters, `None` without."""
        return None if self.uncertainty is None else self.uncertainty.id

    def fitted_values(self, size: int) -> dict[str, Quantity]:
        """Get the fitted parameters at their fitted values as a dimension of a scan.

        A scan of a sampled population (the variability of the blood flow or
        the body weight) shows the spread of the samples; in a run over the
        draws of the fitted parameters it would run every sample for every
        draw. The fitted values, set by the dimension of the samples, keep the
        scan out of the draws (`pkpdsim.uncertainty.apply_uncertainty`), so it
        runs once with the fitted model.

        Args:
            size: the number of points of the dimension.

        Returns:
            In a run over the draws the values of the model of the drawn
            parameters, `size` times; nothing otherwise.

        Raises:
            TypeError: if the model of the experiment is not loaded.
        """
        run = self.uncertainty
        if run is None:
            return {}
        model = self._models["model"]
        if not isinstance(model, RoadrunnerSBMLModel) or model.r is None:
            raise TypeError(f"The model of '{self.sid}' is not loaded: '{model}'.")
        return {
            pid: Q(np.full(size, float(model.r[pid])), model.uinfo[pid])
            for pid in run.dimension.values
        }

    def scan_values(
        self, res: ScanResult, sid: str, unit: str | None = None
    ) -> np.ndarray:
        """Get a variable at the start of every point of a scan of one dimension.

        In a run over the draws of the fitted parameters the values of the first
        draw: a scanned parameter is the same in every draw.

        Args:
            res: the result of the scan.
            sid: the variable, e.g. the scanned parameter.
            unit: the unit of the values, the unit of the result by default.

        Returns:
            The values at the points of the scan.
        """
        values = start_values(res, sid)
        across = self._across()
        if across is not None and across in values.dims:
            values = values.isel({across: 0})
        quantity = Q(np.asarray(values.values, dtype=float), res.units[sid])
        return quantity.to(unit).magnitude if unit else quantity.magnitude

    def pk_band(
        self, df: pd.DataFrame, key: str, unit: str | None = None
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
        """Get a parameter of `calculate_icg_pk` at the points of a scan.

        Args:
            df: the parameters of the scan.
            key: the parameter.
            unit: the unit of the values, the unit of the parameter by default.

        Returns:
            The values at the points (the median of the draws in a run over
            them), and the lower and the upper quantile, `None` without draws.
        """
        source = str(df[f"{key}_unit"].iloc[0])

        def convert(values: np.ndarray) -> np.ndarray:
            quantity = Q(np.asarray(values, dtype=float), source)
            return quantity.to(unit).magnitude if unit else quantity.magnitude

        run = self.uncertainty
        if run is None:
            return convert(df.sort_values("point")[key].to_numpy()), None, None
        summary = summarize(df, key, run.quantiles)
        return (
            convert(summary["median"].to_numpy()),
            convert(summary["lower"].to_numpy()),
            convert(summary["upper"].to_numpy()),
        )

    def scan_at_end(
        self, res: ScanResult, sid: str, unit: str | None = None
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
        """Get a variable at the last time of every point of a scan of one dimension.

        E.g. the steady state of an infusion.

        Args:
            res: the result of the scan.
            sid: the variable.
            unit: the unit of the values, the unit of the result by default.

        Returns:
            The values at the points (the median of the draws in a run over
            them), and the lower and the upper quantile, `None` without draws.
        """
        values = end_values(res, sid)

        def convert(array: np.ndarray) -> np.ndarray:
            quantity = Q(np.asarray(array, dtype=float), res.units[sid])
            return quantity.to(unit).magnitude if unit else quantity.magnitude

        run = self.uncertainty
        if run is None or run.id not in values.dims:
            return convert(values.values), None, None
        axis = values.dims.index(run.id)
        array = np.asarray(values.values, dtype=float)
        low, high = run.quantiles
        return (
            convert(np.nanmedian(array, axis=axis)),
            convert(np.nanquantile(array, low, axis=axis)),
            convert(np.nanquantile(array, high, axis=axis)),
        )

    def plot_band(
        self,
        ax: Axes,
        x: np.ndarray,
        values: np.ndarray,
        lower: np.ndarray | None,
        upper: np.ndarray | None,
        *,
        color: str,
        label: str | None = None,
        marker: str = "o",
        linestyle: str = "-",
        linewidth: float = 1.5,
        markersize: float = 6.0,
    ) -> np.ndarray:
        """Plot values against `x` with the band of the draws, if any.

        Args:
            ax: the axes.
            x: the values of the scanned variable at the points.
            values: the values at the points (the median of the draws).
            lower: the lower quantile of the draws, `None` without draws.
            upper: the upper quantile of the draws, `None` without draws.
            color: the color of the line and the band.
            label: the label of the line.
            marker: the marker of the points.
            linestyle: the style of the line.
            linewidth: the width of the line.
            markersize: the size of the markers.

        Returns:
            The values.
        """
        if lower is not None and upper is not None:
            ax.fill_between(
                x, lower, upper, color=color, alpha=self.band_alpha, linewidth=0
            )
        ax.plot(
            x,
            values,
            color=color,
            label=label,
            marker=marker,
            linestyle=linestyle,
            linewidth=linewidth,
            markersize=markersize,
        )
        return values

    def plot_pk(
        self,
        ax: Axes,
        x: np.ndarray,
        df: pd.DataFrame,
        key: str,
        unit: str | None = None,
        *,
        color: str,
        label: str | None = None,
        marker: str = "o",
        linestyle: str = "-",
        linewidth: float = 1.5,
        markersize: float = 6.0,
    ) -> np.ndarray:
        """Plot a parameter of a scan against `x`, with the band of the draws.

        Args:
            ax: the axes.
            x: the values of the scanned variable at the points.
            df: the parameters of the scan.
            key: the parameter.
            unit: the unit of the values, the unit of the parameter by default.
            color: the color of the line and the band.
            label: the label of the line.
            marker: the marker of the points.
            linestyle: the style of the line.
            linewidth: the width of the line.
            markersize: the size of the markers.

        Returns:
            The values at the points, the median of the draws in a run over them.
        """
        values, lower, upper = self.pk_band(df, key, unit)
        return self.plot_band(
            ax,
            x,
            values,
            lower,
            upper,
            color=color,
            label=label,
            marker=marker,
            linestyle=linestyle,
            linewidth=linewidth,
            markersize=markersize,
        )

    def with_sensitivity(
        self, simulations: Mapping[str, Simulation | Scan]
    ) -> dict[str, Simulation | Scan]:
        """Add a parameter sensitivity scan to every simulation which is not a scan.

        The scan `<key>_sensitivity` of a simulation `<key>` has the dimension
        `dim_sens`, the local design of the parameters of the model
        (`sbmlsim.simulation.sampling.local`): the reference, then every
        parameter alone changed by +25 % and by -25 %, in the order of the ids.
        The parameters are the constant parameters of the model which
        `exclude_parameters_icg` does not leave out, which the simulation does
        not change at a time and which are not zero after the
        pre-initialization of the simulation
        (`sbmlsim.simulation.sampling.parameters_of`).

        Args:
            simulations: the simulations by key.

        Returns:
            The simulations and their sensitivity scans.
        """
        model = self._models["model"]
        if not isinstance(model, RoadrunnerSBMLModel):
            raise TypeError(f"The model of '{self.sid}' is not loaded: '{model}'.")
        scans: dict[str, Simulation | Scan] = dict(simulations)
        for sim_key, sim in simulations.items():
            if not isinstance(sim, Simulation):
                continue
            # a scan sets its value wherever the simulation sets the target,
            # so a target which the simulation changes at a time (the
            # infusion rates of Soons1991) is the protocol and not scanned
            changed = {t for change in sim.changes for t in change.values}
            targets = parameters_of(
                model,
                exclude=lambda pid, changed=changed: (
                    exclude_parameters_icg(pid) or pid in changed
                ),
                simulation=sim,
            )
            scans[self.to_sid(f"{sim_key}_sensitivity")] = Scan(
                sim,
                # 0.75 * value and 1.25 * value
                [
                    local(
                        targets, delta=0.25, model=model, simulation=sim, id="dim_sens"
                    )
                ],
            )
        return scans
