"""The tolerances of the integrator for the indocyanine green model."""

from sbmlsim.model.tolerances import AbsoluteTolerance


def absolute_tolerance(value: float, drug: float) -> AbsoluteTolerance:
    """Get the absolute tolerance of the integrator for the indocyanine green model.

    The species of indocyanine green have plasma concentrations of about 1e-2 mM
    after a bolus of 0.5 mg/kg, which fall within an hour to the secondary
    phase of about 1e-4 mM and in cirrhosis or after long simulations below;
    the extraction ratio of an infusion is the difference of two concentrations
    of about 1e-3 mM. The doses (the targets of the rate rules) are milligrams.

    Args:
        value: tolerance of the targets of the rate rules, i.e. the doses.
        drug: tolerance of the species, of the amount and of the concentration
            [mmole, mM]; the tolerance of a concentration is multiplied by the
            volume of its compartment.

    Returns:
        The absolute tolerance per state.
    """
    return AbsoluteTolerance(amount=drug, concentration=drug, other=value)


#: the absolute tolerance of the simulation experiments
ABSOLUTE_TOLERANCE = absolute_tolerance(1e-10, drug=1e-13)

#: the relative tolerance of the simulation experiments
RELATIVE_TOLERANCE = 1e-10
