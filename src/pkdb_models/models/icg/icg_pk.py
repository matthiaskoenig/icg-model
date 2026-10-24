"""Pharmacokinetic parameters of indocyanine green from the timecourses of a scan.

The parameters are those of indocyanine green in venous plasma after an
intravenous bolus. They are of two kinds:

- the non-compartmental analysis of `pkpdutils` (e.g. `auc_last`,
  `auc_inf_obs`, `cmax`, `thalf`, `lambda_z`, `cl`, `vz`), named like the
  variables of its `NCAResult`, with its terminal phase, the best fit of the
  last points; with the reflux of indocyanine green from the liver into the
  plasma this is the slow secondary phase of a long simulation. The
  concentrations below `LLOQ_FRACTION` of the maximum are below the limit of
  quantification, so that the noise of the integrator is not the terminal
  phase;
- the parameters of the liver function test, from the log-linear regression of
  the concentrations in a window after the bolus (`PDR_WINDOW`, 2 - 15 min by
  default, a study with another window passes its own):
  - `kel` [1/min]: the elimination rate, the negative slope of the regression;
  - `icgpdr` [%/min]: the plasma disappearance rate, `kel` in %/min;
  - `thalf_kel` [min]: the half-life of the plasma disappearance, ln 2 / `kel`;
  - `c0` [mM]: the concentration at the time of the bolus, the back
    extrapolation of the regression;
  - `v0` [l]: the volume of the dose over `c0`, the plasma volume;
  - `cl_kel` [l/min]: the clearance of the liver function test, `kel` * `v0`;
  - `icgr15`, `icgr20` [%]: the retention, the concentration at 15 min and
    20 min relative to `c0`.

The clearances from the area under the curve are `cl_auc` (the dose over
`auc_last`, i.e. to the end of the simulation), `cl_auc_perkg` (per body weight
of the simulation) and `cl_blood` (the blood clearance, `cl` over `1 - HCT` of
the simulation). The unit of every parameter is in `<name>_unit`.

A scan of one dimension, or a simulation, gives one row per point; a run over
the draws of the fitted parameters (`run_simulation pkdb_models.models.icg
--uncertainty`) adds the dimension of the draws. The rows have the position of
their point along the dimension of the scan (`point`, 0 for a simulation) and,
with draws, of the draw (`draw`); `summarize` gives the median and the
quantiles over the draws.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from pkpdsim.pk import nca_frame
from pkpdsim.scan import scan_dimensions, start_values
from pkpdutils import NCAOptions, Route, Timecourses, nca, ureg
from sbmlsim import Q
from sbmlsim.result import ScanResult

if TYPE_CHECKING:
    from pkdb_models.models.icg.experiments.base_experiment import (
        ICGSimulationExperiment,
    )

logger = logging.getLogger(__name__)

#: the intravenous dose
DOSE = "IVDOSE_icg"
#: the concentration of the analysis, indocyanine green in venous plasma
CONCENTRATION = "[Cve_icg]"
#: the limit of quantification relative to the maximal concentration
LLOQ_FRACTION = 1e-4
#: the window of the log-linear regression of the liver function test [min]
PDR_WINDOW: tuple[float, float] = (2.0, 15.0)
#: the times of the retentions [min]
RETENTIONS: dict[str, float] = {"icgr15": 15.0, "icgr20": 20.0}

#: options of the non-compartmental analysis: the terminal phase is the best
#: fit of the last points (the default of pkpdutils)
NCA_OPTIONS = NCAOptions()

#: units of the clearances and volumes of the non-compartmental analysis
UNITS: dict[str, str] = {
    "cl": "l/min",
    "cl_pred": "l/min",
    "vz": "l",
    "vz_pred": "l",
    "vss": "l",
    "vss_pred": "l",
}

#: the selections `calculate_icg_parameters` reads from the result of a scan; an
#: experiment with large scans which uses nothing else records only these
PK_SELECTIONS: tuple[str, ...] = ("time", DOSE, CONCENTRATION, "BW", "HCT")


def _rows(df: pd.DataFrame, res: ScanResult, across: str | None) -> pd.DataFrame:
    """Get the rows of the parameters of a result, in the order of the points.

    The result is a simulation (`point` 0) or a scan of one dimension, besides
    the draws `across`; its position is the column `point`, the position of
    the draw the column `draw`; the rows are ordered by the point and then the
    draw.

    Raises:
        ValueError: if the result has more than one dimension besides the draws.
    """
    dims = scan_dimensions(res, across)
    if len(dims) > 1:
        raise ValueError(
            f"A scan of one dimension is required, not of {dims} (draws: {across})."
        )
    names = {
        **dict.fromkeys(dims, "point"),
        **({across: "draw"} if across is not None else {}),
    }
    df = df.rename(columns=names)
    if not dims:
        df.insert(0, "point", 0)
    positions = ["point", *(["draw"] if across is not None else [])]
    df = df.sort_values(positions, kind="stable", ignore_index=True)
    return df[[*positions, *(c for c in df.columns if c not in positions)]]


def log_linear_regression(
    time: np.ndarray, values: np.ndarray, window: tuple[float, float]
) -> tuple[float, float]:
    """Get the elimination rate and the back extrapolation of a window.

    Args:
        time: the time points.
        values: the concentrations at the time points.
        window: the first and the last time of the regression.

    Returns:
        The elimination rate (the negative slope of the log concentrations) and
        the concentration extrapolated to the time 0; NaN for less than three
        positive concentrations in the window.
    """
    mask = (time >= window[0]) & (time <= window[1]) & (values > 0)
    if mask.sum() < 3:
        return np.nan, np.nan
    slope, intercept = np.polyfit(time[mask], np.log(values[mask]), 1)
    return float(-slope), float(np.exp(intercept))


def calculate_icg_parameters(
    experiment: ICGSimulationExperiment,
    xres: ScanResult,
    across: str | None = None,
    window: tuple[float, float] = PDR_WINDOW,
) -> pd.DataFrame:
    """Calculate the ICG parameters (PDR, R15, ...) of every simulation of a scan.

    The non-compartmental analysis of all simulations at once.

    Args:
        experiment: experiment of the scan, for the molecular weight.
        xres: result of a simulation or of a scan of one dimension.
        across: the dimension of the draws of the fitted parameters, `None`
            for a result without draws (or a result which does not run over
            them, a simulation which sets the drawn parameters itself).
        window: the window of the regression of the liver function test [min].

    Returns:
        One row per simulation, with its position (`point`, and `draw` with
        draws), the parameters and their units (`<name>_unit`).
    """
    if across is not None and across not in xres.dims:
        # a simulation which sets every drawn parameter itself runs once
        across = None
    units = xres.units
    dims = tuple(xres.dims)
    # the dose, the body weight and the hematocrit at the start of the
    # simulation of every point
    dose_mmole = (
        Q(start_values(xres, DOSE).values, units[DOSE]) / experiment.Mr.icg
    ).to("mmole")
    ds = xres.to_timecourses(CONCENTRATION).ds
    time = (
        Q(ds["times"].values if "times" in ds else ds["time"].values, units["time"])
        .to("min")
        .magnitude
    )
    c_vec = Q(ds["value"].values, units[CONCENTRATION]).to("mM").magnitude
    timecourses = Timecourses.from_arrays(
        time,
        c_vec,
        time_unit="min",
        unit="mM",
        dims=dims,
        coords={dim: ds[dim].values for dim in dims} or None,
        dose={"amount": dose_mmole.magnitude, "unit": "mmol"},
        route=Route.IV_BOLUS,
        substance="indocyanine green",
    )
    values = timecourses.ds["value"]
    cmax = values.max(dim=[d for d in values.dims if d not in dims])
    lloq = (LLOQ_FRACTION * cmax).where(cmax > 0)
    result = nca(Timecourses(timecourses.ds.assign(lloq=lloq)), options=NCA_OPTIONS)
    nca_df = nca_frame(result)
    position = tuple(nca_df[dim].to_numpy() for dim in dims)

    def at(array: np.ndarray) -> np.ndarray:
        """Get the values of an array over the dimensions in the rows of the analysis."""
        return np.reshape(np.asarray(array)[position], -1)

    # the samples with a value of every simulation, the times of the others NaN
    sampled = np.isfinite(c_vec)
    t_vec = np.where(sampled, np.broadcast_to(time, c_vec.shape), np.nan)

    # the liver function test: the log-linear regression of the window, of
    # every simulation (cheap, a regression of a few points)
    kel_values = np.empty(c_vec.shape[:-1])
    c0_values = np.empty(c_vec.shape[:-1])
    for index in np.ndindex(c_vec.shape[:-1]):
        kel_values[index], c0_values[index] = log_linear_regression(
            t_vec[index][sampled[index]], c_vec[index][sampled[index]], window
        )
    kel = ureg.Quantity(at(kel_values), "1/min")
    c0 = ureg.Quantity(at(c0_values), "mM")

    # the quantities which are combined with the parameters are quantities of
    # the unit registry of pkpdutils
    dose_pk = ureg.Quantity(at(dose_mmole.magnitude), "mmole")
    bodyweight = ureg.Quantity(at(start_values(xres, "BW").values), str(units["BW"]))
    hematocrit = ureg.Quantity(at(start_values(xres, "HCT").values), str(units["HCT"]))
    names = [name for name in nca_df.columns if f"{name}_unit" in nca_df.columns]
    pk = {
        name: ureg.Quantity(nca_df[name].to_numpy(), nca_df[f"{name}_unit"].iloc[0])
        for name in names
    }
    pk = {
        key: value.to(UNITS[key]) if key in UNITS else value
        for key, value in pk.items()
    }

    v0 = (dose_pk / c0).to("l")
    retentions = {}
    for key, t_retention in RETENTIONS.items():
        index = np.nanargmin(np.abs(t_vec - t_retention), axis=-1)[..., None]
        t_closest = at(np.take_along_axis(t_vec, index, -1)[..., 0])
        for _ in range(int(np.sum(np.abs(t_closest - t_retention) >= 0.1))):
            logger.warning("%s: no time point at %s min", key, t_retention)
        c_retention = at(np.take_along_axis(c_vec, index, -1)[..., 0])
        retentions[key] = (ureg.Quantity(c_retention, "mM") / c0).to("percent")

    cl_auc = (dose_pk / pk["auc_last"]).to("l/min")
    parameters = {
        "dose": dose_pk,
        **pk,
        "kel": kel,
        "icgpdr": kel.to("percent/min"),
        "thalf_kel": (np.log(2) / kel).to("min"),
        "c0": c0,
        "v0": v0,
        "cl_kel": (kel * v0).to("l/min"),
        **retentions,
        "cl_auc": cl_auc,
        "cl_auc_perkg": cl_auc / bodyweight,
        "cl_blood": pk["cl"] / (1 - hematocrit),
    }
    columns: dict[str, object] = {dim: nca_df[dim].to_numpy() for dim in dims}
    for key, value in parameters.items():
        columns[key] = value.magnitude
        columns[f"{key}_unit"] = str(value.units)
    return _rows(pd.DataFrame(columns), xres, across)


def summarize(
    df: pd.DataFrame, key: str, quantiles: tuple[float, float]
) -> pd.DataFrame:
    """Get the median and the quantiles of a parameter over the draws of every point.

    Args:
        df: the parameters with the columns `point` and `draw`.
        key: the parameter.
        quantiles: the lower and the upper quantile.

    Returns:
        One row per point (the index), the columns `lower`, `median`, `upper`.
    """
    low, high = quantiles
    grouped = df.groupby("point")[key]
    return pd.DataFrame(
        {
            "lower": grouped.quantile(low),
            "median": grouped.median(),
            "upper": grouped.quantile(high),
        }
    ).astype(np.float64)
