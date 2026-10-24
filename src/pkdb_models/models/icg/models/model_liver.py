"""Liver model of indocyanine green: the hepatic uptake, the reflux and the biliary excretion."""

from sbmlutils.factory import (
    Compartment,
    FactoryResult,
    Model,
    ModelUnits,
    Parameter,
    Reaction,
    Species,
    UnitDefinition,
    Units,
    create_model,
)
from sbmlutils.metadata import BQB, SBO

from pkdb_models.models.icg.models import annotations, templates
from pkdb_models.models.icg.models.fits import FITS


class U(Units):
    """Units of the liver model."""

    min = UnitDefinition("min")
    mmole = UnitDefinition("mmole")
    m2 = UnitDefinition("m2", "meter^2")
    mmole_per_min_l = UnitDefinition("mmole_per_min_l", "mmole/min/l")
    mM = UnitDefinition("mM", "mmole/liter")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")
    per_min = UnitDefinition("per_min", "1/min")


_m = Model(
    "icg_liver",
    name="model for icg elimination by the liver",
    notes="""
    # Model for icg elimination by the liver

    ## Description
    Model for hepatic elimination of indocyanine green
    encoded in <a href="http://sbml.org">SBML</a> format.

    **Assumptions:**

    - indocyanine green is not metabolized; it is taken up by the liver
      (OATP1B3) and excreted unchanged into bile (MRP2)
    - a part of the indocyanine green in the liver returns into the sinusoidal
      plasma (reflux), which is the slow secondary phase of the plasma curve
    - bilirubin inhibits the uptake competitively

    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=ModelUnits(
        time=U.min,
        extent=U.mmole,
        substance=U.mmole,
        length=U.meter,
        area=U.m2,
        volume=U.liter,
    ),
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),  # liver
        (BQB.OCCURS_IN, "bto/BTO:0000759"),
        (BQB.OCCURS_IN, "ncit/C12392"),
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
        *annotations.model,
    ],
)

_m.compartments = [
    Compartment(
        "Vext",
        4,
        name="plasma",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True,
    ),
    Compartment(
        "Vli",
        1.5,
        name="liver",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
        port=True,
        notes="""
        BW = 75 [kg]
        FVli = 0.0210 [l/kg]
        """,
    ),
    Compartment(
        "Vbi",
        1.5 * 0.00071 / 0.0210,
        name="bile",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["bi"],
        port=True,
        notes="""
        BW = 75 [kg]
        FVbi = 0.00071 [l/kg]
        """,
    ),
    Compartment(
        "Vfeces",
        1.0,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        annotations=annotations.compartments["feces"],
        port=True,
    ),
    Compartment(
        "Vbil",
        1.0,
        name="bilirubin compartment",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
    ),
]

_m.species = [
    Species(
        "icg_ext",
        initialConcentration=0.0,
        name="icg (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["icg"],
        port=True,
    ),
    Species(
        "icg",
        initialConcentration=0.0,
        name="icg (liver)",
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["icg"],
    ),
    Species(
        "icg_bi",
        initialConcentration=0.0,
        name="icg (bile)",
        compartment="Vbi",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["icg"],
    ),
    Species(
        "icg_feces",
        initialConcentration=0.0,
        name="icg (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["icg"],
        port=True,
    ),
    Species(
        sid="bil_ext",
        initialConcentration=0.01,
        name="bilirubin (extern)",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        compartment="Vbil",
        annotations=annotations.species["bil"],
        boundaryCondition=True,
        notes="""
        Serum bilirubin; 0.01 mmole/l = 0.585 mg/dl, the normal range is
        0.3 - 1.2 mg/dl (5 - 21 µmole/l).
        """,
    ),
]

_m.reactions = [
    Reaction(
        "ICGIM",
        name="ICG import",
        equation="icg_ext -> icg [bil_ext]",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vli",
        pars=[
            Parameter(
                "ICGIM_Vmax",
                **FITS.parameter("LI__ICGIM_Vmax", unit=U.mmole_per_min_l),
                name="Vmax for icg import",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "ICGIM_Km",
                0.025,
                unit=U.mM,
                name="Km icg of icg import",
                sboTerm=SBO.MICHAELIS_CONSTANT_FOR_SUBSTRATE,
                notes="""
                Not identified by the data, fixed: the plasma data identify the
                intrinsic clearance Vmax/Km of the uptake, a fit of both found
                Km of 0.016 - 0.34 mM at the same cost.
                """,
            ),
            Parameter(
                "ICGIM_ki_bil",
                0.02,
                unit=U.mM,
                name="Ki bilirubin of icg import",
                sboTerm=SBO.INHIBITORY_CONSTANT,
                notes="""
                Ki of the competitive inhibition by bilirubin, in the range of
                the serum bilirubin: 0.02 mmole/l = 1.17 mg/dl
                (1 mg/dl = 10 mg/l / 584.66 g/mole = 0.0171 mmole/l).
                """,
            ),
            Parameter(
                "f_oatp1b3",
                1,
                unit=U.dimensionless,
                name="scaling factor protein amount",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""
                Parameter for scaling transport protein amount.
                """,
            ),
        ],
        formula=(
            "f_oatp1b3 * ICGIM_Vmax * Vli * icg_ext/"
            "(ICGIM_Km * (1 dimensionless + bil_ext/ICGIM_ki_bil) + icg_ext)",
            U.mmole_per_min,
        ),
        annotations=[
            (BQB.IS, "uniprot/Q9NPD5"),
            (BQB.IS, "omit/0040643"),
            (BQB.IS_ENCODED_BY, "ncit/C106617"),
        ],
        notes="""
        Solute carrier organic anion transporter family member 1B3 (SLCO1B3).
        Often referred to as OATP1B3.
        Bilirubin effect was modeled as competitive inhibition.
        """,
    ),
    # reflux from the liver into the sinusoidal plasma
    Reaction(
        "ICGEX",
        name="ICG export (plasma)",
        equation="icg -> icg_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vli",
        pars=[
            Parameter(
                "ICGEX_k",
                **FITS.parameter("LI__ICGEX_k", unit=U.per_min),
                name="rate of the icg reflux into plasma",
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
            ),
        ],
        formula=("ICGEX_k * Vli * icg", U.mmole_per_min),
        annotations=[(BQB.IS_VERSION_OF, "uniprot/O15438")],
        notes="""
        Sinusoidal efflux of indocyanine green from the hepatocytes into the
        plasma (the bidirectional transport of the OATPs, MRP3 (ABCC3)), first
        order in the hepatic amount. It is the slow secondary phase of the
        plasma curve after a bolus (Meijer1988: half-life 63 - 89 min).
        """,
    ),
    # transport from the liver into bile
    Reaction(
        "ICGLI2CA",
        name="ICG export (bile)",
        equation="icg -> icg_bi",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vli",
        pars=[
            Parameter(
                "ICGLI2CA_Vmax",
                **FITS.parameter("LI__ICGLI2CA_Vmax", unit=U.mmole_per_min_l),
                name="Vmax bile export icg",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "ICGLI2CA_Km",
                **FITS.parameter("LI__ICGLI2CA_Km", unit=U.mM),
                name="Km bile export icg",
                sboTerm=SBO.MICHAELIS_CONSTANT_FOR_SUBSTRATE,
                notes="""
                The saturation of the export: the maximal biliary excretion is
                about the same at 1 and 2 mg/kg (Meijer1988), and at low
                hepatic concentrations Vmax/Km is the rate of the slow terminal
                phase of the plasma curve.
                """,
            ),
        ],
        formula=("ICGLI2CA_Vmax * Vli * icg/(ICGLI2CA_Km + icg)", U.mmole_per_min),
        annotations=[(BQB.IS, "uniprot/Q92887")],
        notes="""
        Canalicular export of indocyanine green by the multidrug resistance
        associated protein 2, MRP2 (ABCC2), saturable.
        """,
    ),
    Reaction(
        "ICGLI2BI",
        name="ICG bile transport",
        equation="icg_bi -> icg_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vbi",
        pars=[
            Parameter(
                "ICGLI2BI_k",
                **FITS.parameter("LI__ICGLI2BI_k", unit=U.per_min),
                sboTerm=SBO.FORWARD_RATE_CONSTANT,
                name="rate of the icg bile transport",
            ),
        ],
        formula=("ICGLI2BI_k * Vbi * icg_bi", U.mmole_per_min),
        notes="""
        The transit of bile from the canaliculi through the bile ducts, first
        order in the amount in bile, i.e. the excretion of indocyanine green
        with the bile; as simplification directly into the feces.
        """,
    ),
]

model_liver: Model = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.icg import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
