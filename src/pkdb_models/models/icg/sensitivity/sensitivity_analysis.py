"""Sensitivity analysis of the indocyanine green model.

The outputs are the pharmacokinetic parameters of an intravenous bolus of
indocyanine green in venous plasma: the non-compartmental analysis of
`pkpdutils` (AUC, Cmax, half-life, Vz and clearance, the terminal phase of the
60 min of the simulation) and the parameters of the liver function test, the
log-linear regression of the concentrations from 2 to 15 min after the bolus
(`icg_pk.log_linear_regression`): the plasma disappearance rate (ICG-PDR), the
retention after 15 min relative to the back extrapolated concentration
(ICG-R15) and the clearance of the test (the elimination rate times the dose
over the back extrapolated concentration). The parameters are the parameters
of the model with bounds of +-15 % around their value, the fitted parameters
and the scaling factors of the hepatic blood flow, the cardiac output and the
uptake included. The analyses are run for the control, the stages of
cirrhosis and resections of the liver. They are run as a script, the local
analysis by default:

    python -m pkdb_models.models.icg.sensitivity.sensitivity_analysis
    python -m pkdb_models.models.icg.sensitivity.sensitivity_analysis local sampling

The analyses run on `pkpdsim.sensitivity.SensitivityRunner`, which writes the
table of the parameters into `results/sensitivity` of the model and the tables
and figures of an analysis into `results/sensitivity/<analysis>`.
An analysis does not vary a parameter which its simulation or one of its groups
sets (`pkpdsim.sensitivity.free_parameters`; the dose `IVDOSE_icg` and the
changes of the groups are no parameters of the analyses). The environment
variable `PKPDSIM_SENSITIVITY_N` replaces the numbers of samples of the global
analyses for a quick run, see `pkpdsim.sensitivity.sample_size`.
"""

from collections.abc import Mapping
from typing import Any

import numpy as np
from pkpdsim.sensitivity import (
    FastAnalysis,
    LocalAnalysis,
    ParameterType,
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityParameter,
    SensitivityRunner,
    SobolAnalysis,
    parameters_from_sbml,
    set_bounds,
    set_relative_bounds,
)
from pkpdutils import NCAOptions
from sbmlsim import Q
from sbmlsim.simulation import PK, Custom, Formula, Observable, Simulation

from pkdb_models.models.icg import MODEL, MODEL_PATH
from pkdb_models.models.icg.experiments.base_experiment import (
    ICGSimulationExperiment,
    exclude_parameters_icg,
)
from pkdb_models.models.icg.fitting.parameters import parameters_pk as fit_parameters
from pkdb_models.models.icg.icg_pk import (
    CONCENTRATION,
    PDR_WINDOW,
    RETENTIONS,
    log_linear_regression,
)

#: intravenous dose of the analysis [mg], 0.5 mg/kg of 75 kg
DOSE_ICG = 37.5
#: molecular weight of indocyanine green [g/mole], see `ICGSimulationExperiment.Mr`
MR_ICG = 774.96493
#: the dose of the analysis [mmole]
DOSE_ICG_MMOLE = DOSE_ICG / MR_ICG
#: fraction of the bounds of a parameter relative to its value
BOUNDS_FRACTION = 0.15
#: the scaling factors of the functions which the analyses vary
SCALING_FACTORS: tuple[str, ...] = ("f_bloodflow", "f_cardiac_output", "LI__f_oatp1b3")
#: the end [min] and the steps of the simulation (output every 0.05 min)
END = 60
STEPS = 1200
#: the limit of quantification of the non-compartmental analysis [mM], 1e-4 of
#: the maximal venous concentration after 37.5 mg (about 0.015 mM), so that the
#: noise of the integrator is not the terminal phase
LLOQ = 1.5e-6

#: outputs of the non-compartmental analysis of pkpdutils -> unit of the output
PK_OUTPUTS: dict[str, str] = {
    "auc_inf_obs": "mM*min",
    "cmax": "mM",
    "thalf": "min",
    "vz": "l",
    "cl": "ml/min",
}

#: the outputs with the ids and names of the analyses; the observable of an
#: output is `icg_<id>`
OUTPUTS: list[SensitivityOutput] = [
    SensitivityOutput(uid, name, f"icg_{uid}")
    for uid, name in [
        ("auc_inf_obs", "ICG AUC∞"),
        ("cmax", "ICG Cmax"),
        ("thalf", "ICG Half-life"),
        ("vz", "ICG Vz"),
        ("cl", "ICG CL"),
        ("icgpdr", "ICG-PDR"),
        ("icgr15", "ICG-R15"),
        ("cl_kel", "ICG CL (PDR)"),
    ]
]


def _regression(time: np.ndarray, values: dict[str, Any]) -> tuple[float, float]:
    """Get the elimination rate [1/min] and the back extrapolation [mM] of the test."""
    return log_linear_regression(
        np.asarray(time), np.asarray(values[CONCENTRATION]), PDR_WINDOW
    )


def icgpdr(time: np.ndarray, values: dict[str, Any]) -> float:
    """Get the plasma disappearance rate [%/min], the regression of 2 to 15 min."""
    kel, _ = _regression(time, values)
    return 100 * kel


def icgr15(time: np.ndarray, values: dict[str, Any]) -> float:
    """Get the retention after 15 min [%], relative to the back extrapolation."""
    _, c0 = _regression(time, values)
    concentration = np.interp(
        RETENTIONS["icgr15"], np.asarray(time), np.asarray(values[CONCENTRATION])
    )
    return 100 * float(concentration) / c0


def cl_kel(time: np.ndarray, values: dict[str, Any]) -> float:
    """Get the clearance of the test [ml/min], the rate times dose over c0."""
    kel, c0 = _regression(time, values)
    return 1000 * kel * DOSE_ICG_MMOLE / c0


#: the observables of a simulation, the outputs with their units: the analysis
#: of indocyanine green in venous plasma with the dose in mmole (the dose of the
#: simulations over the molecular weight) and the parameters of the liver
#: function test, functions of the timecourse of every simulation (`Custom`)
OBSERVABLES: list[Observable] = [
    PK(
        "pk_icg",
        CONCENTRATION,
        dose=Q(DOSE_ICG_MMOLE, "mmole"),
        route="iv_bolus",
        options=NCAOptions(lloq=LLOQ),
        parameters=list(PK_OUTPUTS),
    ),
    *(
        Formula(f"icg_{key}", f"pk_icg.{key}", unit=unit)
        for key, unit in PK_OUTPUTS.items()
    ),
    Custom("icg_icgpdr", icgpdr, unit="percent/min", symbols=[CONCENTRATION]),
    Custom("icg_icgr15", icgr15, unit="percent", symbols=[CONCENTRATION]),
    Custom("icg_cl_kel", cl_kel, unit="ml/min", symbols=[CONCENTRATION]),
]


def sensitivity_simulation(changes: Mapping[str, float]) -> Simulation:
    """Get the simulation of a group: an intravenous bolus of indocyanine green."""
    return Simulation(
        time_unit="min",
        end=END,
        preinit_changes={"IVDOSE_icg": Q(DOSE_ICG, "mg"), **changes},
    )


#: the groups the analyses are run for: the control, the stages of cirrhosis of
#: the simulation experiments and resections of the liver
SENSITIVITY_GROUPS: list[SensitivityGroup] = [
    SensitivityGroup("control", "Control", sensitivity_simulation({}), "dimgrey"),
    *(
        SensitivityGroup(
            uid,
            name,
            sensitivity_simulation(
                {
                    "f_tissue_loss": ICGSimulationExperiment.cirrhosis_map[name],
                    "f_shunts": ICGSimulationExperiment.cirrhosis_map[name],
                }
            ),
            ICGSimulationExperiment.cirrhosis_colors[name],
        )
        for uid, name in [
            ("mildCirrhosis", "Mild cirrhosis"),
            ("moderateCirrhosis", "Moderate cirrhosis"),
            ("severeCirrhosis", "Severe cirrhosis"),
        ]
    ),
    *(
        SensitivityGroup(
            f"resection{round(100 * rate)}",
            f"Resection {round(100 * rate)} %",
            sensitivity_simulation({"resection_rate": rate}),
            color,
        )
        for rate, color in [
            (0.2, "#fdae6b"),
            (0.4, "#fd8d3c"),
            (0.6, "#e6550d"),
            (0.8, "#a63603"),
        ]
    ),
]


def sensitivity_parameters() -> list[SensitivityParameter]:
    """Get the parameters of the analyses with their bounds and references."""
    # the parameters which are not changed in the sensitivity of the simulation
    # experiments (molecular weight, the degrees of cirrhosis and resection, the
    # injection time), the dosing and the unused volume of the feces; the
    # scaling factors of the functions (`SCALING_FACTORS`) are varied
    parameters = [
        p
        for p in parameters_from_sbml(MODEL_PATH, exclude_na=True, exclude_zero=True)
        if (p.uid in SCALING_FACTORS or not exclude_parameters_icg(p.uid))
        and p.uid not in {"IVDOSE_icg", "Ri_icg", "Vfeces"}
    ]

    # the fitted parameters, their bounds are set relative to the value below
    set_bounds(
        parameters,
        [(fp.pid, np.nan, np.nan, ParameterType.FIT) for fp in fit_parameters],
    )
    # the scaling factors of the hepatic blood flow, the cardiac output and the
    # hepatic uptake
    set_bounds(
        parameters,
        [
            (uid, 1 - BOUNDS_FRACTION, 1 + BOUNDS_FRACTION, ParameterType.SCALING)
            for uid in SCALING_FACTORS
        ],
    )

    # the references of the values from the literature
    references = {
        "HCT": r"\cite{Mondal2025, Fiseha2023}",
        "BW": r"\cite{Ogden2004, Jones2013, Thompson2009, Brown1997}",
        "FQgi": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FQh": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
        "FVar": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVgi": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVli": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
        "FVlu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVve": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "COBW": r"\cite{Cattermole2017, Patel2021, Collis2001}",
    }
    by_uid = {p.uid: p for p in parameters}
    for uid, reference in references.items():
        p = by_uid[uid]
        p.reference = reference
        if p.type == ParameterType.NA:
            p.type = ParameterType.DATA

    # every parameter without bounds varies around its value
    set_relative_bounds(parameters, BOUNDS_FRACTION)
    return parameters


#: the analyses of the script, with their groups and numbers of samples
SENSITIVITY = SensitivityRunner(
    model=MODEL,
    parameters=sensitivity_parameters,
    outputs=OUTPUTS,
    observables=OBSERVABLES,
    analyses={
        "local": LocalAnalysis(groups=SENSITIVITY_GROUPS, delta=0.01),
        "sampling": SamplingAnalysis(groups=SENSITIVITY_GROUPS, n=1000),
        "fast": FastAnalysis(groups=SENSITIVITY_GROUPS[:1], n=1000),
        "sobol": SobolAnalysis(groups=SENSITIVITY_GROUPS[:1], n=4096),
    },
    default=("local",),
)


if __name__ == "__main__":
    SENSITIVITY.main()
