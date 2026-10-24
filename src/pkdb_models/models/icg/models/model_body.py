"""Whole-body PBPK model of indocyanine green, which contains the liver model."""

import logging
from pathlib import Path
from typing import NotRequired, TypedDict

from sbmlutils.factory import (
    PORT_SUFFIX,
    AssignmentRule,
    Compartment,
    ExternalModelDefinition,
    Model,
    ModelUnits,
    Parameter,
    RateRule,
    Reaction,
    ReplacedElement,
    Species,
    Submodel,
    UnitDefinition,
    Units,
    create_model,
)
from sbmlutils.metadata import BQB, SBO

from pkdb_models.models.icg.models import annotations, templates

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------------------------
# Whole body model
# -------------------------------------------------------------------------------------
class U(Units):
    """Units of the whole-body model."""

    mmole = UnitDefinition("mmole", "mmole")
    min = UnitDefinition("min", "min")
    kg = UnitDefinition("kg", "kg")
    m = UnitDefinition("m", "meter")
    m2 = UnitDefinition("m2", "meter^2")
    mM = UnitDefinition("mM", "mmole/liter")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")
    mmole_per_min_l = UnitDefinition("mmole_per_min_l", "mmole/min/l")
    l_per_min = UnitDefinition("l_per_min", "l/min")
    per_min = UnitDefinition("per_min", "1/min")
    s = UnitDefinition("s", "second")
    mg = UnitDefinition("mg", "mg")
    ml = UnitDefinition("ml", "ml")
    hr = UnitDefinition("hr", "hr")
    per_hr = UnitDefinition("per_hr", "1/hr")
    cm = UnitDefinition("cm", "cm")
    mg_per_l = UnitDefinition("mg_per_l", "mg/l")
    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    mg_per_hr = UnitDefinition("mg_per_hr", "mg/hr")
    l_per_hr = UnitDefinition("l_per_hr", "l/hr")
    l_per_kg = UnitDefinition("l_per_kg", "l/kg")
    l_per_ml = UnitDefinition("l_per_ml", "l/ml")
    g_per_mole = UnitDefinition("g_per_mole", "g/mole")
    m2_per_kg = UnitDefinition("m2_per_kg", "meter^2/kg")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")
    mmole_per_hr = UnitDefinition("mmole_per_hr", "mmole/hr")
    mmole_per_min_kg = UnitDefinition("mmole_per_min_kg", "mmole/min/kg")
    mmole_per_hr_ml = UnitDefinition("mmole_per_hr_ml", "mmole/hr/ml")
    ml_per_s = UnitDefinition("ml_per_s", "ml/s")
    ml_per_s_kg = UnitDefinition("ml_per_s_kg", "ml/s/kg")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    mul_per_g = UnitDefinition("mul_per_g", "microliter/g")
    mul_per_min_mg = UnitDefinition("mul_per_min_mg", "microliter/min/mg")
    min_per_hr = UnitDefinition("min_per_hr", "min/hr")
    s_per_min = UnitDefinition("s_per_min", "s/min")
    s_per_hr = UnitDefinition("s_per_hr", "s/hr")


_m = Model(
    "icg_body",
    name="whole-body PBPK model of ICG",
    notes="""
    # Whole-body PBPK model of ICG

    ## Description
    Model for whole-body distribution and elimination of indocyanine green
    encoded in <a href="http://sbml.org">SBML</a> format.

    **Assumptions:**

    - indocyanine green is confined to the plasma (no distribution into the
      tissues); the convective transport is the plasma flow, the blood flow
      times (1 - HCT)
    - indocyanine green is not metabolized; it is taken up by the liver and
      excreted unchanged into bile

    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=ModelUnits(
        time=U.min,
        extent=U.mmole,
        substance=U.mmole,
        length=U.m,
        area=U.m2,
        volume=U.liter,
    ),
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "BTO:0001489"),  # whole body
        (BQB.OCCURS_IN, "FMA:256135"),  # body
        (BQB.HAS_PROPERTY, "ncit/C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
        (BQB.HAS_PROPERTY, "ncit/C79370"),  # Pharmacokinetics: Distribution
        *annotations.model,
    ],
)


class SubstanceBody(TypedDict):
    """Definition of a substance transported in the whole-body model."""

    name: str
    unit: UnitDefinition
    cinit: float  # initial concentration [mmole/l]
    Mr: float  # molecular weight [g/mole]
    IVDOSE: NotRequired[float]  # intravenous dose [mg]
    annotations: list[tuple[BQB, str]]


SUBSTANCES_BODY: dict[str, SubstanceBody] = {
    "icg": {
        "name": "icg",
        "unit": U.mmole,
        # initial concentration
        "cinit": 0.0,  # [mmole/l]
        # Molecular weight
        "Mr": 774.96493,  # [g/mole]
        # doses [mg]
        "IVDOSE": 0,  # flag to add injection & infusion kinetics to the model
        "annotations": annotations.species["icg"],
    },
}

# -----------------------------------------------------------------------------
# Submodels
# -----------------------------------------------------------------------------
COMPARTMENTS_BODY = {
    "bi": "bile",
    "re": "rest",
    "gi": "gastrointestinal tract",
    "li": "liver",
    "lu": "lung",
    "ve": "venous plasma",
    "ar": "arterial plasma",
    "po": "portal vein",
    "hv": "hepatic vein",
}

SUBMODEL_SID_DICT = {  # tissue to submodel mapping
    "li": "LI",  # liver
}

liver_id = "icg_liver"
_m.external_model_definitions = [
    ExternalModelDefinition(
        sid="liver",
        source=f"{liver_id}.xml",
        modelRef=liver_id,
        name="liver model definition",
    ),
]

_m.submodels = [
    Submodel(
        sid=SUBMODEL_SID_DICT["li"],
        modelRef="liver",
        name="liver submodel",
    ),
]

for emd in _m.external_model_definitions:
    logger.info("%s (%s)", emd, Path(emd.source).absolute())


# -------------------------------------------------------------------------------------------------
# Compartments
# -------------------------------------------------------------------------------------------------
_m.compartments = [
    Compartment(
        "Vbi",
        metaId="meta_Vbi",
        value=1,
        unit=U.liter,
        constant=False,
        name="bile",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["bi"],
    ),
    Compartment(
        "Vre",
        metaId="meta_Vre",
        value=1,
        unit=U.liter,
        constant=False,
        name="rest of body",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["re"],
    ),
    Compartment(
        "Vgi",
        metaId="meta_Vgi",
        value=1,
        unit=U.liter,
        constant=False,
        name="gastrointestinal tract",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["gi"],
    ),
    Compartment(
        "Vli",
        metaId="meta_Vli",
        value=1,
        unit=U.liter,
        constant=False,
        name="liver",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
    ),
    Compartment(
        "Vlu",
        metaId="meta_Vlu",
        value=1,
        unit=U.liter,
        constant=False,
        name="lung",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["lu"],
    ),
    Compartment(
        "Vve",
        metaId="meta_Vve",
        value=1,
        unit=U.liter,
        constant=False,
        name="venous plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ve"],
    ),
    Compartment(
        "Var",
        metaId="meta_Var",
        value=1,
        unit=U.liter,
        constant=False,
        name="arterial plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ar"],
    ),
    Compartment(
        "Vpo",
        metaId="meta_Vpo",
        value=1,
        unit=U.liter,
        constant=False,
        name="portal plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["po"],
    ),
    Compartment(
        "Vhv",
        metaId="meta_Vhv",
        value=1,
        unit=U.liter,
        constant=False,
        name="hepatic venous plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["hv"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["feces"],
    ),
]

# create plasma and tissue compartments (for correct blood volume)
for cid in COMPARTMENTS_BODY:
    if cid not in ["ar", "ve", "po", "hv", "bi"]:
        _m.compartments.extend(
            [
                Compartment(
                    f"V{cid}_tissue",
                    value=1,
                    unit=U.liter,
                    constant=False,
                    name=f"{cid} tissue",
                    metaId=f"meta_V{cid}_tissue",
                    port=True,
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                    annotations=annotations.compartments["parenchyma"]
                    + [
                        (BQB.IS_PART_OF, aid[1])
                        for aid in annotations.compartments[cid]
                    ],
                ),
                Compartment(
                    f"V{cid}_plasma",
                    value=1,
                    unit=U.liter,
                    constant=False,
                    name=f"{cid} plasma",
                    metaId=f"meta_V{cid}_plasma",
                    port=True,
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                    annotations=annotations.compartments["plasma"]
                    + [
                        (BQB.IS_PART_OF, aid[1])
                        for aid in annotations.compartments[cid]
                    ],
                ),
            ]
        )

# replace volumes of submodels
_m.replaced_elements = [
    # liver
    ReplacedElement(
        sid="Vli_tissue_RE",
        metaId="Vli_tissue_RE",
        elementRef="Vli_tissue",
        submodelRef=SUBMODEL_SID_DICT["li"],
        portRef=f"Vli{PORT_SUFFIX}",
    ),
    ReplacedElement(
        sid="Vli_plasma_RE",
        metaId="Vli_plasma_RE",
        elementRef="Vli_plasma",
        submodelRef=SUBMODEL_SID_DICT["li"],
        portRef=f"Vext{PORT_SUFFIX}",
    ),
    # bile
    ReplacedElement(
        sid="Vbi_RE",
        metaId="Vbi_RE",
        elementRef="Vbi",
        submodelRef=SUBMODEL_SID_DICT["li"],
        portRef=f"Vbi{PORT_SUFFIX}",
    ),
    # feces
    ReplacedElement(
        sid="Vfeces_RE",
        metaId="Vfeces_RE",
        elementRef="Vfeces",
        submodelRef=SUBMODEL_SID_DICT["li"],
        portRef=f"Vfeces{PORT_SUFFIX}",
    ),
]

# -------------------------------------------------------------------------------------
# Species
# -------------------------------------------------------------------------------------
_m.species = []
for sid, sdict in SUBSTANCES_BODY.items():
    # plasma species
    for cid, cname in COMPARTMENTS_BODY.items():
        if cid == "bi":
            continue
        if cid in ["ve", "ar", "po", "hv"]:
            sid_ex = f"C{cid}_{sid}"
            cid_ex = f"V{cid}"
        else:
            # plasma compartment
            sid_ex = f"C{cid}_plasma_{sid}"
            cid_ex = f"V{cid}_plasma"

        _m.species.append(
            Species(
                sid_ex,
                metaId=f"meta_{sid_ex}",
                initialConcentration=sdict["cinit"],
                compartment=cid_ex,
                substanceUnit=sdict["unit"],
                name=f"{sdict['name']} ({cname})",
                hasOnlySubstanceUnits=False,
                port=True,
                annotations=sdict["annotations"],
                sboTerm=SBO.SIMPLE_CHEMICAL,
            )
        )

    if sid in ["icg"]:
        _m.species.append(
            Species(
                f"Afeces_{sid}",
                metaId=f"meta_Afeces_{sid}",
                initialConcentration=0,
                compartment="Vfeces",
                substanceUnit=sdict["unit"],
                name=f"{sdict['name']} (feces)",
                hasOnlySubstanceUnits=True,
                annotations=sdict["annotations"],
                sboTerm=SBO.SIMPLE_CHEMICAL,
            )
        )

# replace species
replaced_species = {
    "li": ["icg"],  # ['icg', 'bil'],
}
# plasma to external species
for tkey, skey_list in replaced_species.items():
    for skey in skey_list:
        _m.replaced_elements.extend(
            [
                ReplacedElement(
                    sid=f"C{tkey}_plasma_{skey}_RE",
                    metaId=f"C{tkey}_plasma_{skey}_RE",
                    elementRef=f"C{tkey}_plasma_{skey}",
                    submodelRef=SUBMODEL_SID_DICT[tkey],
                    portRef=f"{skey}_ext{PORT_SUFFIX}",
                )
            ]
        )

# replace ICG in feces
_m.replaced_elements.append(
    ReplacedElement(
        sid="Afeces_icg_RE",
        metaId="Afeces_icg_RE",
        elementRef="Afeces_icg",
        submodelRef=SUBMODEL_SID_DICT["li"],
        portRef=f"icg_feces{PORT_SUFFIX}",
    )
)

# -------------------------------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------------------------------
_m.parameters = [
    # --- whole body data ---
    Parameter(
        "BW",
        75,
        U.kg,
        constant=True,
        name="body weight [kg]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS, "efo/0004338"),
            (BQB.IS, "ncit/C25208"),
            (BQB.IS, "cmo/CMO:0000012"),
        ],
    ),
    Parameter(
        "HEIGHT",
        170,
        U.cm,
        constant=True,
        name="height [cm]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS, "ncit/C25347"),
            (BQB.IS, "efo/0004339"),
        ],
    ),
    Parameter(
        "BSA",
        0,
        U.m2,
        constant=False,
        name="body surface area [m^2]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[(BQB.IS, "ncit/C25157"), (BQB.IS, "omit/0003188")],
        notes="""

        """,
    ),
    Parameter(
        "COBW",
        75 / 60,
        U.ml_per_s_kg,
        constant=True,
        name="cardiac output per bodyweight [ml/s/kg]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS, "ncit/C119246"),
            (BQB.IS, "omit/0003637"),
        ],
        notes="""Simone1997 75 [50-100] ml/min/kg""",
    ),
    Parameter(
        "CO",
        75 / 60 * 75,
        U.ml_per_s,
        constant=False,
        name="cardiac output [ml/s]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS, "ncit/C119246"),
            (BQB.IS, "omit/0003637"),
        ],
        notes="""BW*COBW*f_cardiac_output (93.75 ml/s at 75 kg)""",
    ),
    Parameter(
        "QC",
        75 / 60 * 75 / 1000 * 60,
        U.l_per_min,
        constant=False,
        name="cardiac output [L/min]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS, "ncit/C119246"),
            (BQB.IS, "omit/0003637"),
        ],
        notes="""
        CO/1000 ml/l * 60 s/min (5.6 l/min at 75 kg)
        ICRP2002 reference values Cardiac output: 6.5 l/min (male); 5.9 l/min (female)
        """,
    ),
    Parameter(
        "Fblood",
        0.02,
        U.dimensionless,
        constant=False,
        name="blood fraction of organ volume",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        A certain volume of an organ is occupied by large vessels (in the low
        percent range). This parameter describes this fraction.
        """,
    ),
    Parameter(
        "HCT",
        0.45,
        U.dimensionless,
        constant=True,
        metaId="meta_HCT",
        name="hematocrit",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        Normal adult hematocrit 0.40 - 0.50 (men about 0.45, women about 0.41;
        Gadano1997 0.44). Indocyanine green is confined to the plasma: the
        hematocrit sets the plasma volume (43 ml/kg at 0.45) and the plasma
        flows, which are the convective flows of the model.
        """,
        annotations=[
            (BQB.IS, "ncit/C64796"),
            (BQB.IS, "omit/0007571"),
            (BQB.IS, "efo/0004348"),
        ],
    ),
    # --- fractional tissue volumes ---
    Parameter(
        "FVgi",
        0.0297,
        U.l_per_kg,
        constant=True,
        name="fractional tissue volume gastrointestinal tract",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""FVgi = FVgu + FVpa + FVsp = 0.0171 + 0.01 + 0.0026 = 0.0297""",
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["gi"]
        ],
    ),
    Parameter(
        "FVbi",
        0.00071,
        U.l_per_kg,
        constant=True,
        name="fractional tissue volume bile",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["bi"]
        ],
    ),
    Parameter(
        "FVli",
        0.0210,
        U.l_per_kg,
        constant=True,
        name="fractional tissue volume liver",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        1.7-1.8 kg (male); 1.5 kg (female); 75kg; 1.7/75 = 0.0227
        """,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["li"]
        ],
    ),
    Parameter(
        "FVlu",
        0.0076,
        U.l_per_kg,
        constant=True,
        name="fractional tissue volume lung",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["lu"]
        ],
    ),
    # 20% / 64%
    # FVve + FVar = 0.0771
    # 64/84; FVve = 0.0587
    # 20/84; Far  = 0.0184
    Parameter(
        "FVve",
        0.0587,
        U.l_per_kg,
        constant=True,
        name="fractional tissue volume venous",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["ve"]
        ],
    ),
    Parameter(
        "FVar",
        0.0184,
        U.l_per_kg,
        constant=True,
        name="fractional tissue volume arterial",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["ar"]
        ],
    ),
    Parameter(
        "FVpo",
        0.001,
        U.l_per_kg,
        constant=True,
        name="fractional tissue volume portal",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["po"]
        ],
    ),
    Parameter(
        "FVhv",
        0.001,
        U.l_per_kg,
        constant=True,
        name="fractional tissue volume hepatic venous",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["hv"]
        ],
    ),
    Parameter(
        "FVre",
        0,
        U.l_per_kg,
        constant=False,
        name="fractional tissue volume rest of body",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["re"]
        ],
        notes="""
        Fraction is calculated based on other volume fractions.
        """,
    ),
    # fractional tissue blood flows
    # FQgi = FQgu + FQpa + FQsp = 0.146 + 0.017 + 0.017 = 0.18
    # 20 % arteriell; 80 % portal vein (gi) => 0.172
    # FQgi < FQh
    # IRCP2001 reference value arterial: 0.065 (male/female);
    # FQgi = FQpo = FQh - FQha = 0.255 - 0.065 = 0.190
    Parameter(
        "FQgi",
        0.190,
        U.dimensionless,
        constant=True,
        name="gastrointestinal tract fractional tissue blood flow",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["gi"]
        ],
    ),
    # IRCP2001 reference value total liver: 0.255 (male); 0.270 (female)
    Parameter(
        "FQh",
        0.255,
        U.dimensionless,
        constant=True,
        name="hepatic (venous side) fractional tissue blood flow",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["hv"]
        ],
    ),
    Parameter(
        "FQlu",
        1,
        U.dimensionless,
        constant=True,
        name="lung fractional tissue blood flow",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["lu"]
        ],
    ),
    Parameter(
        "FQre",
        0,
        U.dimensionless,
        constant=False,
        name="rest of body fractional tissue blood flow",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, item[1]) for item in annotations.compartments["re"]
        ],
        notes="""
        Calculated based on other tissue blood flow and flow conservation laws.
        """,
    ),
    # --- changes for cirrhotic liver ---
    # the degree of cirrhosis sets the shunts and the tissue loss together
    # (`cirrhosis_map` of the experiments)
    Parameter(
        "f_shunts",
        0,
        U.dimensionless,
        constant=True,
        name="fraction of portal venous blood shunted by the liver",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS, "hp/HP:0002629"),  # Gastrointestinal arteriovenous malformation
        ],
    ),
    Parameter(
        "f_tissue_loss",
        0,
        U.dimensionless,
        constant=True,
        name="fraction of lost parenchymal liver volume",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, "ncit/C2951"),  # Cirrhosis
            (BQB.IS_VERSION_OF, "efo/0001422"),  # cirrhosis of liver
        ],
    ),
    # --- changes for blood flow experiments ---
    Parameter(
        "f_bloodflow",
        1,
        U.dimensionless,
        constant=True,
        name="fractional change of hepatic blood flow",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        Scales the hepatic and the gastrointestinal blood flow. It must be
        below 1/FQh = 3.92, otherwise the hepatic blood flow exceeds the cardiac
        output and the flow of the rest of the body is negative.
        """,
    ),
    Parameter(
        "f_cardiac_output",
        1,
        U.dimensionless,
        constant=True,
        name="fractional change of cardiac output",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
    # --- changes for hepatectomy experiments ---
    Parameter(
        "resection_rate",
        0,
        U.dimensionless,
        constant=True,
        name="resection rate",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS_VERSION_OF, "ncit/C158758"),  # Resection
        ],
    ),
]

# species specific parameters
for sid, sdict in SUBSTANCES_BODY.items():
    _m.parameters.extend(
        [
            # molecular weights
            Parameter(
                f"Mr_{sid}",
                sdict["Mr"],
                U.g_per_mole,
                constant=True,
                name=f"Molecular weight {sid} [g/mole]",
                sboTerm=SBO.MOLECULAR_MASS,
                annotations=[
                    (BQB.IS_VERSION_OF, "opb/OPB_01146"),  # Molecular weight
                    (BQB.HAS_PART, sdict["annotations"][0][1]),  # icg
                ],
            ),
        ]
    )

    if "IVDOSE" in sdict:
        _m.parameters.extend(
            [
                # dosing
                Parameter(
                    f"IVDOSE_{sid}",
                    sdict["IVDOSE"],
                    U.mg,
                    constant=False,
                    name=f"IV bolus dose {sid} [mg]",
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    annotations=[(BQB.IS, "ncit/C38274")],
                ),
                # iv kinetics after application
                Parameter(
                    f"ti_{sid}",
                    5,
                    U.s,
                    constant=True,
                    name=f"injection time {sid} [s]",
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    annotations=[(BQB.IS, "ncit/C69282")],
                ),
                Parameter(
                    f"Ki_{sid}",
                    0.02,
                    U.per_min,
                    constant=False,
                    name=f"Ki [1/min] injection {sid}",
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                ),
                # continuous infusion
                Parameter(
                    f"Ri_{sid}",
                    0,
                    U.mg_per_min,
                    constant=True,
                    name=f"Ri [mg/min] rate of infusion {sid}",
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    annotations=[(BQB.IS, "ncit/C94916")],
                ),
            ]
        )

# -------------------------------------------------------------------------------------------------
# AssignmentRules
# -------------------------------------------------------------------------------------------------
_m.rules = [
    AssignmentRule(
        "FVre",
        "1.0 l_per_kg - (FVbi + FVgi + FVli + FVlu + FVve + FVar + FVpo + FVhv)",
        U.l_per_kg,
        name="rest body volume",
        notes="""
        Remaining volume based on fractional organ volume is used for the rest
        compartment.
        """,
    ),
    AssignmentRule(
        "FQre",
        "1.0 dimensionless - Qh/Qlu",
        U.dimensionless,
        name="rest body perfusion",
        notes="""
        Rest body perfusion .must be calculate based on flux conservation laws.
        ```
        Qlu = Qre + Qh
        Qre = Qlu - Qh
        FQre = (Qlu - Qh) / Qlu = 1 - Qh/Qlu
        ```
        """,
    ),
    AssignmentRule(
        "BSA",
        "0.024265 m2 * power(BW/1 kg, 0.5378) * power(HEIGHT/1 cm, 0.3964)",
        U.m2,
        name="body surface area (BSA)",
        notes="""
        Calculated with Haycock1978 formula.
        """,
    ),
    AssignmentRule(
        "CO",
        "BW*COBW*f_cardiac_output",
        U.ml_per_s,
        name="cardiac output",
    ),
    AssignmentRule(
        "QC", "CO/1000 ml_per_l * 60 s_per_min", U.l_per_min, name="cardiac blood flow"
    ),
    # volumes
    AssignmentRule("Vbi", "BW*FVbi", U.liter, name="volume bile"),
    AssignmentRule("Vgi", "BW*FVgi", U.liter, name="volume gastrointestinal tract"),
    AssignmentRule(
        "Vli",
        "BW*FVli*(1 dimensionless - resection_rate)",
        U.liter,
        name="volume liver",
    ),
    AssignmentRule("Vlu", "BW*FVlu", U.liter, name="volume lung"),
    AssignmentRule("Vre", "BW*FVre", U.liter, name="volume rest body"),
    # venous and arterial plasma volume (corrected for tissue blood volumes)
    AssignmentRule(
        "Vve",
        "(1 dimensionless - HCT) * (BW*FVve - FVve/(FVar+FVve+FVpo+FVhv) * BW * "
        "Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv)))",
        U.liter,
        name="volume venous plasma",
    ),
    AssignmentRule(
        "Var",
        "(1 dimensionless - HCT) * (BW*FVar - FVar/(FVar+FVve+FVpo+FVhv) * BW * "
        "Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv)))",
        U.liter,
        name="volume arterial plasma",
    ),
    AssignmentRule(
        "Vpo",
        "(1 dimensionless - HCT) * (BW*FVpo - FVpo/(FVar+FVve+FVpo+FVhv) * BW * "
        "Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv)))",
        U.liter,
        name="volume hepatic portal plasma",
    ),
    AssignmentRule(
        "Vhv",
        "(1 dimensionless - HCT) * (BW*FVhv - FVhv/(FVar+FVve+FVpo+FVhv) * "
        "BW * Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv)))",
        U.liter,
        name="volume hepatic venous plasma",
    ),
    # blood flows
    AssignmentRule("Qlu", "QC*FQlu", U.l_per_min, name="blood flow lung"),
    AssignmentRule("Qre", "QC*FQre", U.l_per_min, name="blood flow rest of body"),
    AssignmentRule(
        "Qh",
        "QC*FQh*f_bloodflow",
        U.l_per_min,
        name="blood flow hepatic (venous side)",
    ),
    AssignmentRule(
        "Qgi",
        "QC*FQgi*f_bloodflow",
        U.l_per_min,
        name="blood flow gastrointestinal tract",
    ),
    AssignmentRule("Qpo", "Qgi", U.l_per_min, name="blood flow portal vein"),
    AssignmentRule("Qha", "Qh - Qpo", U.l_per_min, name="blood flow hepatic artery"),
]

# plasma flows: indocyanine green is confined to the plasma, so its convective
# transport is the plasma flow, the blood flow times (1 - HCT)
PLASMA_FLOWS = {
    "lu": "lung",
    "re": "rest of body",
    "gi": "gastrointestinal tract",
    "h": "hepatic (venous side)",
    "ha": "hepatic artery",
    "po": "portal vein",
}
for qid, qname in PLASMA_FLOWS.items():
    _m.parameters.append(
        Parameter(
            f"Q{qid}_plasma",
            1,
            U.l_per_min,
            constant=False,
            name=f"plasma flow {qname}",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        )
    )
    _m.rules.append(
        AssignmentRule(
            f"Q{qid}_plasma",
            f"Q{qid} * (1 dimensionless - HCT)",
            U.l_per_min,
            name=f"plasma flow {qname}",
        )
    )

# Volumes for explicit tissue models
for cid, cname in COMPARTMENTS_BODY.items():
    if cid not in ["ve", "ar", "po", "hv", "bi"]:
        _m.rules.append(
            AssignmentRule(
                f"V{cid}_plasma",
                value=f"V{cid} * Fblood * (1 dimensionless - HCT)",
                unit=U.liter,
                name=f"plasma volume of {cname}",
                notes="""plasma volume associated with tissue""",
            ),
        )
        if cid == "li":
            _m.rules.append(
                AssignmentRule(
                    f"V{cid}_tissue",
                    value=f"V{cid}*(1 dimensionless - f_tissue_loss) * (1 dimensionless - Fblood)",
                    unit=U.liter,
                    name=f"tissue volume of {cname}",
                    notes="""Cirrhosis: Adjustment of fractional liver tissue volume.""",
                ),
            )
        else:
            _m.rules.append(
                AssignmentRule(
                    f"V{cid}_tissue",
                    value=f"V{cid}*(1 dimensionless - Fblood)",
                    unit=U.liter,
                    name=f"tissue volume of {cname}",
                ),
            )

for sid, sdict in SUBSTANCES_BODY.items():
    sname = sdict["name"]

    if "IVDOSE" in sdict:
        # injection
        _m.rules.extend(
            [
                AssignmentRule(
                    f"Ki_{sid}",
                    f"0.693 dimensionless/ti_{sid} * 60 s_per_min",
                    U.per_min,
                    name="injection rate IV",
                ),
            ]
        )


_m.rules.extend(
    [
        AssignmentRule(
            "ER_icg",
            "(Car_icg - Chv_icg)/(Car_icg + 1E-12 mM)",
            U.dimensionless,
            name="extraction ratio ICG",
            notes="""handling division by zero with very small offset""",
        ),
        AssignmentRule(
            "CLinfusion_icg",
            "(Ri_icg/Mr_icg) / (Cve_icg + 1E-12 mM)",
            U.l_per_min,
            name="infusion clearance (steady state)",
            notes="""
            Clearance calculated during steady state in continuous infusion.
            handling division by zero with very small offset
            """,
        ),
    ]
)

# -------------------------------------------------------------------------------------
# Reactions
# -------------------------------------------------------------------------------------
_m.reactions = []

for sid, sdict in SUBSTANCES_BODY.items():
    sname = sdict["name"]

    if "IVDOSE" in sdict:
        _m.reactions.extend(
            [
                # --------------------
                # iv application
                # --------------------
                Reaction(
                    sid=f"iv_{sid}",
                    name=f"iv {sname}",
                    formula=(f"Ki_{sid}*IVDOSE_{sid}/Mr_{sid}", U.mmole_per_min),
                    equation=f"-> Cve_{sid}",
                    compartment="Vve",
                    sboTerm=SBO.TRANSPORT_REACTION,
                    annotations=[(BQB.IS, "ncit/C38276")],
                ),
            ]
        )

    # --------------------
    # blood flow model
    # --------------------
    for cid, cname in COMPARTMENTS_BODY.items():
        if cid in ["ve", "ar", "po", "hv"]:
            continue
        # --------------------
        # ve -> lung -> ar
        # --------------------
        if cid == "lu":
            rid_in = f"Flow_ve_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"

            rid_out = f"Flow_{cid}_ar_{sid}"
            name_out = f"outflow {cname} {sname}"

            # distribution in plasma volume
            _m.reactions.extend(
                [
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}_plasma*Cve_{sid}", U.mmole_per_min),
                        equation=f"Cve_{sid} -> C{cid}_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(f"Q{cid}_plasma*C{cid}_plasma_{sid}", U.mmole_per_min),
                        equation=f"C{cid}_plasma_{sid} -> Car_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # ------------------------------------------
        # ar -> arorgan -> organ -> ve
        # ------------------------------------------
        if cid in ["re"]:
            rid_in1 = f"Flow_ar_ar{cid}_{sid}"
            rid_in2 = f"Flow_ar{cid}_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"
            rid_out = f"Flow_{cid}_ve_{sid}"
            name_out = f"outflow {cname} {sname}"

            # only distribution in plasma volume
            _m.reactions.extend(
                [
                    Reaction(
                        sid=rid_in1,
                        name=name_in,
                        formula=(f"Q{cid}_plasma*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> C{cid}_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                        annotations=[(BQB.IS, "ncit/C16353")]
                        + [
                            (BQB.IS_PART_OF, aid[1])
                            for aid in annotations.compartments[cid]
                        ],
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(f"Q{cid}_plasma*C{cid}_plasma_{sid}", U.mmole_per_min),
                        equation=f"C{cid}_plasma_{sid} -> Cve_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                        annotations=[(BQB.IS, "ncit/C16353")]
                        + [
                            (BQB.IS_PART_OF, aid[1])
                            for aid in annotations.compartments[cid]
                        ],
                    ),
                ]
            )
        # -----------------------------
        # ar -> arorgan -> organ -> po
        # -----------------------------
        if cid in ["gi"]:
            rid_in1 = f"Flow_ar_ar{cid}_{sid}"
            rid_in2 = f"Flow_ar{cid}_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"

            rid_out = f"Flow_{cid}_po_{sid}"
            name_out = f"outflow {cname} {sname}"

            _m.reactions.extend(
                [
                    # only distribution in plasma volume
                    Reaction(
                        sid=rid_in1,
                        name=name_in,
                        formula=(f"Q{cid}_plasma*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> C{cid}_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                        annotations=[(BQB.IS, "ncit/C16353")]
                        + [
                            (BQB.IS_PART_OF, aid[1])
                            for aid in annotations.compartments[cid]
                        ],
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(f"Q{cid}_plasma*C{cid}_plasma_{sid}", U.mmole_per_min),
                        equation=f"C{cid}_plasma_{sid} -> Cpo_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                        annotations=[(BQB.IS, "ncit/C16353")]
                        + [
                            (BQB.IS_PART_OF, aid[1])
                            for aid in annotations.compartments[cid]
                        ],
                    ),
                ]
            )
        # --------------------
        # ar -> arli -> li
        # po -> li
        # li -> hv -> ve
        # --------------------
        if cid == "li":
            _m.reactions.extend(
                [
                    Reaction(
                        sid=f"Flow_arli_li_{sid}",
                        name=f"arterial inflow liver {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*Qha_plasma*Car_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Car_{sid} -> Cli_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # shunted arterial flow
                    Reaction(
                        sid=f"Flow_arli_hv_{sid}",
                        name="flow arterial shunts",
                        formula=(f"f_shunts*Qha_plasma*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # (unshunted) portal flow
                    Reaction(
                        sid=f"Flow_po_li_{sid}",
                        name=f"outflow po {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*Qpo_plasma*Cpo_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Cpo_{sid} -> Cli_plasma_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # portal shunts
                    Reaction(
                        sid=f"Flow_po_hv_{sid}",
                        name="flow portal shunts",
                        formula=(f"f_shunts*Qpo_plasma*Cpo_{sid}", U.mmole_per_min),
                        equation=f"Cpo_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=f"Flow_li_hv_{sid}",
                        name=f"outflow liver {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*(Qpo_plasma+Qha_plasma)*Cli_plasma_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Cli_plasma_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=f"Flow_hv_ve_{sid}",
                        name=f"outflow hepatic vein {sname}",
                        formula=(f"Qh_plasma*Chv_{sid}", U.mmole_per_min),
                        equation=f"Chv_{sid} -> Cve_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )

# --------------------------------------------------------------------------------------------------
# RateRules
# --------------------------------------------------------------------------------------------------
_m.rate_rules = []

for sid, sdict in SUBSTANCES_BODY.items():
    if "IVDOSE" in sdict:
        _m.rate_rules.extend(
            [
                # injection of dose
                RateRule(
                    f"IVDOSE_{sid}",
                    f"-iv_{sid}*Mr_{sid} + Ri_{sid}",
                    U.mg_per_min,
                    name="change of injected dose",
                    notes="""
                    The injected dose appears in the blood over time.
                    """,
                ),
            ]
        )

model_body: Model = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.icg import MODEL_BASE_PATH

    results = create_model(
        filepath=MODEL_BASE_PATH / f"{model_body.sid}.xml",
        model=model_body,
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
