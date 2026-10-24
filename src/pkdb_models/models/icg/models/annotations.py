"""Annotations of the indocyanine green models: model, compartments and species."""

from pkpdsim import annotations
from sbmlutils.metadata import BQB

model = [
    # taxonomy
    (BQB.HAS_TAXON, "taxonomy/9606"),  # human
    (BQB.HAS_TAXON, "snomedct/337915000"),  # human
    # modelling approach
    (BQB.HAS_PROPERTY, "mamo/MAMO_0000046"),  # ordinary differential equation model
    # biological process explained by model (GO/NCIT)
    (BQB.HAS_PROPERTY, "NCIT:C65913"),  # Indocyanine green
    (BQB.HAS_PROPERTY, "CHEBI:31696"),  # Indocyanine green
    # model relevance to a particular area
    (BQB.HAS_PROPERTY, "OMIT:0009194"),  # Liver function test
    # model reference
    (BQB.IS_DESCRIBED_BY, "doi/10.3389/fphys.2021.730418"),
    (BQB.IS_DESCRIBED_BY, "doi/10.3389/fphys.2021.757293"),
    (BQB.IS_DESCRIBED_BY, "doi/10.5281/zenodo.5552397"),
]

compartments = annotations.compartments

species = {
    "bil": [
        (BQB.IS, "inchikey/BPYKTIZUTYGOLE-IFADSCNNSA-N"),
        (BQB.IS, "chebi/CHEBI:16990"),
        (BQB.IS, "ncit/C305"),
    ],
    "icg": [
        (BQB.IS, "inchikey/MOFVSTNWEDAEEK-UHFFFAOYSA-M"),
        (BQB.IS, "chebi/CHEBI:31696"),
        (BQB.IS, "ncit/C65913"),
    ],
}
