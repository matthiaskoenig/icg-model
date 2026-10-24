"""The simulation experiments of the indocyanine green model and their groups.

They are the registry of `MODEL`, the experiments its commands simulate.
"""

from sbmlsim.experiment import SimulationExperiment

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.classification import ClassificationData
from pkdb_models.models.icg.experiments.misc import (
    CTPscoremapping,
    DoseDependencyExperiment,
    ICGPKBilirubinScan,
    ICGPKBloodFlowScan,
    ICGPKBodyweightScan,
    ICGPKCirrhosisScan,
    ICGPKcorrelation,
    ICGPKDoseScan,
    ICGPKHepatectomyScan,
    ICGPKLiverVolumeScan,
    ICGPKOatp1b3Scan,
    ReferenceTimecourse,
)
from pkdb_models.models.icg.experiments.studies import (
    Andersen1999,
    Branch1976,
    Burns1991,
    Caesar1961,
    Cherrik1960,
    Chijiiwa2000,
    DOnofrio2014,
    Gadano1997,
    Grainger1983,
    Grundmann1992,
    Huet1983,
    Kamimori2000,
    Kanstrup1987,
    Kawasaki1988,
    Keiding1993,
    Klockowski1990,
    Leevy1962,
    Leevy1967,
    Martin1975,
    Martin1976,
    Meijer1988,
    Niemann2000,
    Sakka2006,
    Soons1991,
)
from pkdb_models.models.icg.experiments.studies.hepatectomy import (
    Kawasaki1991,
    Ohwada2006,
    Thomas2015,
)

EXPERIMENTS: dict[str, list[type[SimulationExperiment]]] = {
    "studies": [
        Andersen1999,
        Branch1976,
        Burns1991,
        Caesar1961,
        Cherrik1960,
        Chijiiwa2000,
        DOnofrio2014,
        Gadano1997,
        Grainger1983,
        Grundmann1992,
        Huet1983,
        Kamimori2000,
        Kanstrup1987,
        Kawasaki1988,
        Keiding1993,
        Klockowski1990,
        Leevy1962,
        Leevy1967,
        Martin1975,
        Martin1976,
        Meijer1988,
        Niemann2000,
        Sakka2006,
        Soons1991,
        # hepatectomy
        Kawasaki1991,
        Ohwada2006,
        Thomas2015,
    ],
    "pharmacokinetics": [
        Branch1976,
        Caesar1961,
        Cherrik1960,
        DOnofrio2014,
        Gadano1997,
        Grundmann1992,
        Kawasaki1988,
        Klockowski1990,
        Leevy1967,
        Martin1975,
        Martin1976,
        Meijer1988,
        Sakka2006,
    ],
    "timecourses": [
        Andersen1999,
        Burns1991,
        Chijiiwa2000,
        Grainger1983,
        Grundmann1992,
        Kamimori2000,
        Keiding1993,
        Klockowski1990,
        Leevy1962,
        Meijer1988,
        Niemann2000,
        Soons1991,
    ],
    "dose": [
        Martin1975,
        Martin1976,
        Meijer1988,
    ],
    "bilirubin": [
        Branch1976,
        Caesar1961,
        DOnofrio2014,
        Kawasaki1988,
        Leevy1962,
    ],
    "bloodflow": [
        Gadano1997,
        Huet1983,
        Kanstrup1987,
        Leevy1962,
    ],
    "cirrhosis": [
        Branch1976,
        Burns1991,
        Caesar1961,
        Cherrik1960,
        Huet1983,
        Kawasaki1988,
        Leevy1967,
    ],
    "extraction_ratio": [
        Caesar1961,
        Gadano1997,
        Grainger1983,
        Keiding1993,
        Leevy1962,
    ],
    "bile": [
        Chijiiwa2000,
        Meijer1988,
    ],
    "hepatectomy": [
        Kawasaki1991,
        Ohwada2006,
        Thomas2015,
    ],
    "misc": [
        DoseDependencyExperiment,
        ICGPKcorrelation,
        ReferenceTimecourse,
    ],
    "scans": [
        CTPscoremapping,
        ICGPKBilirubinScan,
        ICGPKBloodFlowScan,
        ICGPKBodyweightScan,
        ICGPKCirrhosisScan,
        ICGPKDoseScan,
        ICGPKHepatectomyScan,
        ICGPKLiverVolumeScan,
        ICGPKOatp1b3Scan,
    ],
    # the survival after hepatectomy: the classification data of the
    # hepatectomy studies, see `classification`
    "classification": [
        ClassificationData,
    ],
}
EXPERIMENTS["all"] = (
    EXPERIMENTS["studies"]
    + EXPERIMENTS["misc"]
    + EXPERIMENTS["scans"]
    + EXPERIMENTS["classification"]
)


if __name__ == "__main__":
    MODEL.run_simulation_experiments(EXPERIMENTS["all"], MODEL.simulation_dir / "all")
