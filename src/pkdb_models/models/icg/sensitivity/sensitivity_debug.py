"""Show the blood flows and the volume of the rest of the body of the model.

A check of the reference values of the parameters of the sensitivity analysis.
"""

import roadrunner

from pkdb_models.models.icg import MODEL_PATH


def main() -> None:
    """Print the blood flows and the fractional volume of the rest of the body."""
    r = roadrunner.RoadRunner(str(MODEL_PATH))

    for sid in ["Qre", "Qlu", "Qh", "Qha", "Qpo", "Qgi"]:
        print(f"{sid} = {r[sid]}")

    print("-" * 80)

    for sid in ["FVre"]:
        print(f"{sid} = {r[sid]}")


if __name__ == "__main__":
    main()
