"""Collect the data of the indocyanine green studies from the PK-DB study folders.

The studies are copied from a `pkdb_data` checkout (`PKDB_DATA`, else
`../pkdb_data`, see `pkpdsim.data.studies_dir`) into the data folder of the
model, see `pkpdsim.data.collect_tsv_files`.
Besides the indocyanine green studies, the sampling of virtual individuals uses
the liver volumes and blood flows of Wynne1989 (galactose).
"""

from pathlib import Path

from pkpdsim.data import collect_tsv_files
from pkpdsim.data import studies_dir as pkdb_studies_dir

from pkdb_models.models.icg import DATA_PATH_GALACTOSE, DATA_PATH_ICG

#: the studies of other substances which the icg model uses
GALACTOSE_STUDIES = {"Wynne1989"}


def collect_data(studies_dir: Path | None = None) -> None:
    """Copy the data of the indocyanine green studies into the data folder of the model.

    Args:
        studies_dir: the `studies` folder of pkdb_data, by default the one of
            the checkout `pkpdsim.data.studies_dir` finds.
    """
    if studies_dir is None:
        studies_dir = pkdb_studies_dir()
    collect_tsv_files(source_dir=studies_dir / "icg", target_dir=DATA_PATH_ICG)
    collect_tsv_files(
        source_dir=studies_dir / "galactose",
        target_dir=DATA_PATH_GALACTOSE,
        filter_study=lambda study: study in GALACTOSE_STUDIES,
    )


if __name__ == "__main__":
    collect_data()
