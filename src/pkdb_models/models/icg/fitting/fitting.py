"""Indocyanine green parameter fitting.

The fit problems of the indocyanine green model: which fit mapping collections
enter a fit and which parameters are adjusted, see `FIT_DEFINITIONS`. The
commands `run_fit`, `run_report` and `run_identifiability` of pkpdsim run them
with the package of the model, a fit is followed by the identifiability
analysis of its parameters (`--no-identifiability` skips it):

    run_fit pkdb_models.models.icg --subset=PK --runs=3 --cores=3 --seed=1234 --method=DE

The identifiability analysis is the profile likelihood of every fitted
parameter around the best parameter set and the Fisher information of the same
set, which end up in the identifiability section of the report of the fit.

The results, figures and the HTML report are written into the fit folder of the
model (`MODEL.fit_dir`), unless `--output_dir` is given. The identifiability of
parameters which were fitted earlier is analysed without fitting again with the
`run_identifiability` command, a report of stored parameters is created with
the `run_report` command:

    run_identifiability pkdb_models.models.icg --subset=PK results/fit/PK_<fit id>/parameters.json
    run_report pkdb_models.models.icg --subset=PK results/fit/PK_<fit id>/parameters.json

The plasma and the bile data are fitted together: with the reflux of
indocyanine green from the liver into the plasma, the parameters of the biliary
excretion act on the secondary phase of the plasma curve.
"""

from pkpdsim.fit import default_fit_settings
from sbmlsim.fit.cli import FitDefinition

from pkdb_models.models.icg import DATA_PATHS, ICG_PATH
from pkdb_models.models.icg.fitting.mapping_collections import f_collections_pk
from pkdb_models.models.icg.fitting.parameters import parameters_pk
from pkdb_models.models.icg.tolerances import absolute_tolerance

#: settings of the fits, stored with the result and read back by the report;
#: the species get a tolerance of their own scale (a dose is about 0.05 mmole,
#: the secondary phase about 1e-4 mM, the extraction ratio of an infusion the
#: difference of two concentrations of about 1e-3 mM), see `tolerances`
FIT_SETTINGS = default_fit_settings(
    absolute_tolerance=absolute_tolerance(1e-6, drug=1e-12),
    relative_tolerance=1e-7,
)

#: the fit problems of the indocyanine green model, the keys are the `--subset`
#: argument
FIT_DEFINITIONS: dict[str, FitDefinition] = {
    # the pharmacokinetics of indocyanine green: plasma and bile
    "PK": FitDefinition(
        mapping_collections=f_collections_pk,
        parameters=parameters_pk,
        base_path=ICG_PATH,
        data_path=DATA_PATHS,
        settings=FIT_SETTINGS,
    ),
}
