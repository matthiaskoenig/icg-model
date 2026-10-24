"""Overview of the parameters and species of the sensitivity analysis.

The reference values of the parameters and of the initial species of the model
which are varied, i.e. which are not excluded by `exclude_parameters_icg` and
are not zero (`sbmlsim.simulation.sampling.parameters_of` and `references`).
"""

from collections.abc import Mapping
from pprint import pprint

from sbmlsim.model import RoadrunnerSBMLModel
from sbmlsim.simulation.sampling import parameters_of, references

from pkdb_models.models.icg import MODEL_PATH
from pkdb_models.models.icg.experiments.base_experiment import exclude_parameters_icg


def magnitudes(values: Mapping[str, object]) -> dict[str, object]:
    """Get the magnitudes of the references, which are quantities or numbers."""
    return {key: getattr(value, "magnitude", value) for key, value in values.items()}


def main() -> None:
    """Print the reference values of the parameters and species."""
    print("Loading model")
    model = RoadrunnerSBMLModel(MODEL_PATH)

    print("Reference dict")
    parameters = parameters_of(model, exclude=exclude_parameters_icg)
    # the species (their amounts) of the targets with the species
    known = set(parameters)
    species = [
        sid
        for sid in parameters_of(model, species=True, exclude=exclude_parameters_icg)
        if sid not in known
    ]
    p_ref = magnitudes(references(model, parameters))
    s_ref = magnitudes(references(model, species))

    pprint(p_ref)
    pprint(s_ref)


if __name__ == "__main__":
    main()
