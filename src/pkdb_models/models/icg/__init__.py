"""Indocyanine green: PBPK model of the liver function test with indocyanine green.

The paths of the model and its model object `MODEL`, which the commands
`run_simulation`, `run_fit`, `run_report` and `run_identifiability` of pkpdsim
run by its package `pkdb_models.models.icg`. The paths are those of the model
which is simulated, the data of the studies, the NHANES data of the sampling of
individuals and the results of the simulations, fits, sampling and the
classification; the folders of the results are created by the code which writes
into them.
"""

from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings
from pkpdsim.release.settings import Author, ReleaseSettings

from pkdb_models.models.icg.tolerances import ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE

ICG_PATH = Path(__file__).parent

MODEL_BASE_PATH = ICG_PATH / "models" / "output" / "models"
MODEL_PATH = MODEL_BASE_PATH / "icg_body_flat.xml"

DATA_PATH_BASE = ICG_PATH / "data"
DATA_PATH_ICG = DATA_PATH_BASE / "icg"
DATA_PATH_GALACTOSE = DATA_PATH_BASE / "galactose"
DATA_PATHS = [
    DATA_PATH_ICG,
]

NHANES_DIR = ICG_PATH / "nhanes"
NHANES_INPUT_PATH = NHANES_DIR / "data" / "nhanes.tsv"
NHANES_PROCESSED_PATH = NHANES_DIR / "data" / "nhanes_processed.tsv"

MODEL = PKPDModel(
    name="icg",
    package=__name__,
    title="Indocyanine green PBPK model",
    path=ICG_PATH,
    model_path=MODEL_PATH,
    data_paths=tuple(DATA_PATHS),
    archive_path=MODEL_BASE_PATH.parent / "icg_model.omex",
    factory=".models.factory:create_icg_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
    # the adopted fit, whose values and uncertainty the factory takes
    fits=("PK_20261020_010006__506c2",),
    simulator=SimulatorSettings(
        absolute_tolerance=ABSOLUTE_TOLERANCE, relative_tolerance=RELATIVE_TOLERANCE
    ),
    # the resolution of the PNG figures; at 600 dpi the large figures need up to
    # 12 GB of memory to render
    figures=FigureSettings(formats=("svg",), fig_dpi=300, legend_fontsize=10),
    # the release of the model repository
    release=ReleaseSettings(
        title="Physiologically based pharmacokinetic (PBPK) model of indocyanine green",
        repository="matthiaskoenig/icg-model",
        doi="10.5281/zenodo.5552397",
        authors=(
            Author(
                "Köller, Adrian",
                "0000-0002-4503-5675",
                "Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology",
            ),
            Author(
                "König, Matthias",
                "0000-0003-1725-179X",
                (
                    "Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver; University of Stuttgart, Institute of Structural Mechanics and Dynamics in Aerospace Engineering; University of Southern Denmark, Clinical Pharmacology, Pharmacy and Environmental Medicine, Department of Public Health; University Hospital Schleswig-Holstein, Campus Lübeck, First Department of Medicine, Lübeck, Germany"
                ),
            ),
        ),
        keywords=(
            "model",
            "icg",
            "indocyanine green",
            "liver function test",
            "hepatectomy",
            "pharmacokinetics",
            "PBPK",
            "physiologically based pharmacokinetic model",
            "SBML",
        ),
        docker_image="matthiaskoenig/icg",
        funding=(
            "Adrian Köller and Matthias König were supported by the Federal Ministry of Education and Research (BMBF, Germany) within the research network Systems Medicine of the Liver (LiSyM) by grant number 031L0054. Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B)."
        ),
        description=(
            "This repository provides the indocyanine green (ICG)"
            " physiologically based pharmacokinetic (PBPK) model of the ICG"
            " liver function tests: the distribution of ICG in the plasma, its"
            " hepatic uptake, its reflux into the plasma and its biliary"
            " excretion, with the effects of cirrhosis, hepatectomy, hepatic"
            " blood flow and bilirubin."
        ),
        first_year=2021,
    ),
)
"""The icg model: its commands, see `pkpdsim.model`."""

RESULTS_PATH_PETAB = MODEL.results_dir / "petab"
RESULTS_PATH_CLASSIFICATION = MODEL.results_dir / "classification"
SAMPLES_PATH = MODEL.results_dir / "sampling"
SAMPLES_INDIVIDUALS_NOPK_PATH = SAMPLES_PATH / "individuals_nopk.tsv"
