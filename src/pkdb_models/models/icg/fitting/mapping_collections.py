"""Parameter fit problems for indocyanine green.

The data of a fit is selected from the fit mappings of the studies in three
steps, see `sbmlsim.fit.helpers`: the filters on the metadata of the mappings
select the training data of the fit and exclude the rest, the outliers are
named once for all mappings and are the training data which is not usable, and
the validation data is the part of the selected data a fit is evaluated on but
not fitted to. The metadata of a mapping describes its curve, what a fit does
with the curve is decided here.

The training data are the plasma curves and extraction ratios of healthy
subjects without an interacting drug, and the biliary excretion of Meijer1988.
The validation data has four parts, each a different question: the
interactions (`INTERACTION_VALIDATION`, the arms with a coadministered drug
whose effect on the hepatic blood flow the model describes with a factor which
is not fitted), the doses (`DOSE_VALIDATION`, doses which the training data do
not cover), the populations (`POPULATION_VALIDATION`, subjects which differ
from the healthy adults of the training) and the limitations
(`LIMITATION_VALIDATION`, data of a population or a process the model does not
describe). Repeated measurements of the subjects of a training curve are
excluded (`AGGREGATE_MAPPINGS`).
"""

from pkpdsim.fit import mapping_metadata
from pkpdsim.metadata import Coadministration, Health
from sbmlsim.console import console
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.fit import FitMapping, FitMappingCollection, MappingKind
from sbmlsim.fit.helpers import FitMappings, MappingFilter

from pkdb_models.models.icg import DATA_PATHS, ICG_PATH
from pkdb_models.models.icg.experiments.studies import (
    Andersen1999,
    Burns1991,
    Caesar1961,
    Chijiiwa2000,
    Gadano1997,
    Grainger1983,
    Grundmann1992,
    Kamimori2000,
    Klockowski1990,
    Leevy1962,
    Meijer1988,
    Niemann2000,
    Soons1991,
)

#: the keys of fit mappings by study: the keys of the mappings are SIds of their
#: experiment, which two studies can share, so a named mapping is named with
#: its study
type NamedMappings = dict[str, set[str]]

#: the studies with fit mappings
EXPERIMENT_CLASSES: list[type[SimulationExperiment]] = [
    Andersen1999,
    Burns1991,
    Caesar1961,
    Chijiiwa2000,
    Gadano1997,
    Grainger1983,
    Grundmann1992,
    Kamimori2000,
    Klockowski1990,
    Leevy1962,
    Meijer1988,
    Niemann2000,
    Soons1991,
]

#: training data outside of the healthy subjects without a coadministration,
#: admitted by name
TRAINING_MAPPINGS: NamedMappings = {
    # the control of a cocktail of antipyrine (1000 mg) and oxazepam (15 mg),
    # neither interacts with the hepatic uptake of indocyanine green
    "Andersen1999": {"fm_Cve_icg_individual"},
    # patients after a cholecystectomy with a T-drain and a normal liver
    # function: the only dose resolved data of the biliary excretion, which
    # identify the canalicular export and the transit of bile, and the
    # secondary phase of the plasma curves (the reflux)
    "Meijer1988": {
        "fm_Cve_icg_0p5",
        "fm_Cve_icg_1p0",
        "fm_Cve_icg_2p0",
        "fm_Cve_icg_0p5_secondary",
        "fm_Cve_icg_1p0_secondary",
        "fm_Cve_icg_2p0_secondary",
        "fm_LI__ICGLI2BI_0p5",
        "fm_LI__ICGLI2BI_1p0",
        "fm_LI__ICGLI2BI_2p0",
    },
}

#: mappings which are not fitted, the data is not usable
OUTLIER_MAPPINGS: NamedMappings = {
    "Leevy1962": {
        # Table I, the arterial concentration: the removal rate of the table is
        # 1.8 - 2 times the stated infusion rate
        "fm_Car_icg_tab1",
        # Fig 1, the timecourse of one subject with an unknown scale
        "fm_Car_icg_fig1",
        "fm_Chv_icg_fig1",
    },
}

#: the validation of the interactions: the arms with a coadministered drug,
#: simulated with the observed change of the clearance as a factor of the
#: hepatic blood flow, which is not fitted
INTERACTION_VALIDATION: NamedMappings = {
    # isradipine and diltiazem, vasodilation (f_bloodflow 1.19 and 1.14)
    "Klockowski1990": {"fm_Cve_icg_isradipine", "fm_Cve_icg_diltiazem"},
    # propranolol, the reduced hepatic blood flow (f_bloodflow 0.69)
    "Niemann2000": {"fm_Car_icg_individual_propranolol"},
}

#: the validation of the doses: the extraction ratio after boluses of 0.15 -
#: 0.5 mg/kg of two subjects, the time is not reported
DOSE_VALIDATION: NamedMappings = {
    "Leevy1962": {"fm_er_0p15_icg", "fm_er_0p25_icg", "fm_er_0p5_icg"},
}

#: the validation of the populations: subjects which differ from the healthy
#: adults of the training
POPULATION_VALIDATION: NamedMappings = {
    # chronic hepatitis B/C without fibrosis, with pethidine: the extraction
    # ratio and the clearance of the infusion
    "Gadano1997": {"fm_er_hepatitis", "fm_cl_hepatitis"},
    # one subject with a cardiac output of 10.6 l/min, the plasma converted
    # from blood with an assumed hematocrit
    "Niemann2000": {"fm_Car_icg_individual"},
}

#: the validation of the limitations of the model: valid data of a population
#: or a process the model does not describe, which would distort the fitted
#: parameters of the processes it does describe
LIMITATION_VALIDATION: NamedMappings = {
    # the slow phase 75 - 600 min near the limit of quantification
    "Burns1991": {"fm_Cve_icg_individual_tail"},
    # the biliary excretion of livers after a biliary obstruction
    "Chijiiwa2000": {"fm_LI__ICGLI2BI"},
}

#: repeated measurements of the subjects of a training curve, excluded from
#: every fit: they would count the subjects twice
AGGREGATE_MAPPINGS: NamedMappings = {
    # the second identical bolus of the volunteers
    "Grundmann1992": {"fm_Cve_icg2"},
    # the follicular and ovulatory phase of the crossover of the same women
    # (no effect of the phase), the luteal phase is fitted
    "Kamimori2000": {"fm_Cve_icg_Fol", "fm_Cve_icg_Ovu"},
}

#: the validation data by question
VALIDATION_QUESTIONS: dict[str, NamedMappings] = {
    "interactions": INTERACTION_VALIDATION,
    "doses": DOSE_VALIDATION,
    "populations": POPULATION_VALIDATION,
    "limitations": LIMITATION_VALIDATION,
}

#: the validation data
VALIDATION_MAPPINGS: NamedMappings = {
    study: set().union(
        *(mappings.get(study, set()) for mappings in VALIDATION_QUESTIONS.values())
    )
    for mappings in VALIDATION_QUESTIONS.values()
    for study in mappings
}


def study(fit_mapping: FitMapping) -> str:
    """Get the study of a fit mapping, the id of its experiment."""
    return fit_mapping.reference.experiment.sid


def named(
    mappings: NamedMappings, fit_mapping_key: str, fit_mapping: FitMapping
) -> bool:
    """Check whether a fit mapping is one of the named mappings of its study."""
    return fit_mapping_key in mappings.get(study(fit_mapping), set())


def check_named(mappings: NamedMappings, fit_mappings: FitMappings) -> None:
    """Check that the named mappings are fit mappings of their studies.

    Raises:
        ValueError: for a study which is no experiment, or a key which is no fit
            mapping of its study.
    """
    for sid, keys in mappings.items():
        experiment = fit_mappings.experiments.get(sid)
        if experiment is None:
            raise ValueError(f"'{sid}' is no experiment of the fits.")
        if unknown := sorted(keys - set(experiment._fit_mappings)):
            raise ValueError(f"'{unknown}' are no fit mappings of '{sid}'.")


def filter_aggregates(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """No repeated measurements of training subjects, see `AGGREGATE_MAPPINGS`."""
    return not named(AGGREGATE_MAPPINGS, fit_mapping_key, fit_mapping)


def filter_training(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Healthy subjects without a coadministration, or named data."""
    if named(TRAINING_MAPPINGS, fit_mapping_key, fit_mapping) or named(
        VALIDATION_MAPPINGS, fit_mapping_key, fit_mapping
    ):
        return True
    metadata = mapping_metadata(fit_mapping)
    return (
        metadata.health == Health.HEALTHY
        and metadata.coadministration == Coadministration.NONE
    )


def mapping_collections(
    filters: MappingFilter | list[MappingFilter],
    validation: NamedMappings = VALIDATION_MAPPINGS,
) -> dict[str, list[FitMappingCollection]]:
    """Select the data of a fit from the fit mappings of the studies.

    The filters select the training data, the mappings which fail them are
    excluded, and so are the repeated measurements of training subjects
    (`AGGREGATE_MAPPINGS`) and every mapping which is neither data of healthy
    subjects without a coadministered drug nor named (`filter_training`). The
    outliers are `OUTLIER_MAPPINGS`, the validation data is
    `VALIDATION_MAPPINGS`. Neither kind enters the cost; both are simulated and
    evaluated with the training data, so a report says how the fitted
    parameters describe data they were not fitted on.

    Args:
        filters: filters of the training data.
        validation: the keys of the validation data by study.

    Returns:
        The fit mapping collections of all kinds by experiment id.

    Raises:
        ValueError: if an outlier is a mapping of no study or of several.
    """
    fit_mappings = FitMappings(
        experiment_classes=EXPERIMENT_CLASSES,
        base_path=ICG_PATH,
        data_path=DATA_PATHS,
    )
    for mappings in [
        TRAINING_MAPPINGS,
        OUTLIER_MAPPINGS,
        AGGREGATE_MAPPINGS,
        validation,
        *VALIDATION_QUESTIONS.values(),
    ]:
        check_named(mappings, fit_mappings)
    # sbmlsim selects the outliers by their keys alone, which must name one
    # mapping of all studies
    outliers = {key for keys in OUTLIER_MAPPINGS.values() for key in keys}
    for key in outliers:
        studies = [
            sid
            for sid, experiment in fit_mappings.experiments.items()
            if key in experiment._fit_mappings
        ]
        if len(studies) != 1:
            raise ValueError(f"The outlier '{key}' is a mapping of {studies}.")
    selection = fit_mappings.select(
        filters=[
            *(filters if isinstance(filters, list) else [filters]),
            filter_aggregates,
            filter_training,
        ],
        outliers=outliers,
        validation=lambda key, fit_mapping: named(validation, key, fit_mapping),
        print_info=False,
    )
    # the counts per experiment and kind, not the single mappings
    selection.print(detail=False)
    # the validation mappings by the question they answer
    selected = {
        (sid, key)
        for sid, kinds in selection.kinds.items()
        for key, kind in kinds.items()
        if kind is MappingKind.VALIDATION
    }

    def pairs(mappings: NamedMappings) -> set[tuple[str, str]]:
        return {(sid, key) for sid, keys in mappings.items() for key in keys}

    console.print(
        "validation: "
        + ", ".join(
            f"{len(selected & pairs(mappings))} {question}"
            for question, mappings in VALIDATION_QUESTIONS.items()
        )
    )
    return selection.collections


def f_collections_pk() -> dict[str, list[FitMappingCollection]]:
    """The pharmacokinetics of indocyanine green: plasma and bile."""
    return mapping_collections([])


if __name__ == "__main__":
    console.rule(title=f_collections_pk.__name__, align="left", style="white")
    f_collections_pk()
