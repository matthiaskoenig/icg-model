"""The scaling factor of the hepatic blood flow (`f_bloodflow`) of a given blood flow.

The hepatic blood flow per body weight of the model is
`COBW * FQh * f_bloodflow`, with the cardiac output per body weight `COBW`
[ml/s/kg] and the fraction of the hepatic blood flow `FQh` [-] of the
whole-body model.
"""

import numpy as np
import pandas as pd

from pkdb_models.models.icg.models.model_body import model_body


def _parameter_value(sid: str) -> float:
    """Get the value of a parameter of the whole-body model.

    Args:
        sid: id of the parameter.

    Raises:
        ValueError: if the parameter has no numeric value.
    """
    pars = {p.sid: p for p in model_body.parameters}
    value = pars[sid].value
    if not isinstance(value, int | float):
        raise ValueError(f"Parameter '{sid}' has no numeric value: '{value}'")
    return value


def f_bloodflow(Qhbw: np.ndarray | pd.Series) -> np.ndarray | pd.Series:
    """Get `f_bloodflow` for a given hepatic blood flow per body weight.

    FIXME: the calculation should be performed with units.

    Args:
        Qhbw: hepatic blood flow per body weight [ml/min/kg].

    Returns:
        The scaling factor of the hepatic blood flow [-].
    """
    # [-] = [ml/min/kg] / [ml/s/kg] * [ml/l] / [s/min] / [ml/l]  =  [ml/l] / [ml/l]
    return (
        Qhbw
        / (_parameter_value("COBW") * 60.0 / 1000.0 * _parameter_value("FQh"))
        / 1000
    )


if __name__ == "__main__":
    qhbws = np.array([12.75, 15, 20, 25])
    print(f_bloodflow(Qhbw=qhbws))

    Qh = 12.75 * 75  # [ml/min]
    print("Qh", Qh)
