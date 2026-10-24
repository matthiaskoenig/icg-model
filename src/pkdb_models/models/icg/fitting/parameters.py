"""The fit parameters of the indocyanine green model.

Indocyanine green is confined to the plasma, taken up by the liver (OATP1B3,
competitively inhibited by bilirubin), partly returned into the plasma (the
reflux, the slow secondary phase of the plasma curve) and excreted unchanged
into bile (MRP2, saturable), from where it leaves with the bile. The plasma
data (timecourses, extraction ratios) identify the intrinsic clearance of the
uptake (`LI__ICGIM_Vmax` / `LI__ICGIM_Km`), not its saturation: fitted
together, the Km took any value of 0.016 - 0.34 mM at the same cost, so the Km
is fixed and Vmax carries the intrinsic clearance. The secondary phase of the
plasma curves (Meijer1988 30 - 60 min) identifies the reflux and the rate of
the canalicular export at low hepatic concentrations (its Vmax/Km), the bile
data of Meijer1988 (0.5, 1 and 2 mg/kg) the capacity of the export and the
transit of bile. Fixed are the Km of the uptake, the Ki of bilirubin and the
bilirubin baseline (only `Km * (1 + bil/Ki)` is identified by data of healthy
subjects) and the physiology (cardiac output, hepatic blood flow, hematocrit,
liver volume).
"""

from sbmlsim.fit import FitParameter

#: the pharmacokinetics: plasma and bile
parameters_pk: list[FitParameter] = [
    # hepatic uptake (OATP1B3): the intrinsic clearance Vmax/Km * Vli/(1 + bil/Ki)
    # with an extraction ratio of 0.6 - 0.8 (Leevy1962, Caesar1961, Grainger1983)
    # and a plasma disappearance rate of 18 - 22 %/min (Cherrik1960, Martin1976)
    FitParameter(
        pid="LI__ICGIM_Vmax",
        start_value=0.1,
        lower_bound=1e-3,
        upper_bound=10,
        unit="mmole/min/l",
    ),
    # reflux into the plasma: the secondary phase (Meijer1988 t1/2 63 - 89 min,
    # 0.9 - 1.9 µM at 30 - 60 min)
    FitParameter(
        pid="LI__ICGEX_k",
        start_value=0.04,
        lower_bound=1e-4,
        upper_bound=0.5,
        unit="1/min",
    ),
    # canalicular export (MRP2): the saturated biliary excretion, Meijer1988
    # Emax 8.6 / 13.5 / 12.3 µmol per 0.5 h at 0.5 / 1 / 2 mg/kg
    FitParameter(
        pid="LI__ICGLI2CA_Vmax",
        start_value=4e-4,
        lower_bound=1e-5,
        upper_bound=1e-2,
        unit="mmole/min/l",
    ),
    # the saturation of the canalicular export: Emax of the biliary excretion
    # saturates between 1 and 2 mg/kg (Meijer1988), and at low hepatic
    # concentrations Vmax/Km sets the slow terminal phase of the plasma curve
    # (Meijer1988 t1/2 63 - 89 min, Burns1991 75 - 180 min)
    FitParameter(
        pid="LI__ICGLI2CA_Km",
        start_value=0.015,
        lower_bound=1e-3,
        upper_bound=1,
        unit="mM",
    ),
    # transit of bile: tmax of the biliary excretion 1.3 h at 0.5 mg/kg
    # (Meijer1988)
    FitParameter(
        pid="LI__ICGLI2BI_k",
        start_value=0.025,
        lower_bound=1e-3,
        upper_bound=1,
        unit="1/min",
    ),
]
