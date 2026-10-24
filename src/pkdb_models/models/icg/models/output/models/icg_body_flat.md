# whole-body PBPK model of ICG

Model `icg_body`, SBML Level 3 Version 2, read from icg\_body\_flat.xml, written by sbmlode 0.3.0.

Whole-body PBPK model of ICG

Description

Model for whole-body distribution and elimination of indocyanine green encoded in SBML format.

Assumptions:

indocyanine green is confined to the plasma (no distribution into the tissues); the convective transport is the plasma flow, the blood flow times (1 - HCT)

indocyanine green is not metabolized; it is taken up by the liver and excreted unchanged into bile

The content of this model has been carefully created in a manual research effort.

Terms of use

Copyright © 2026 Matthias König.

This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.

Redistribution and use of any part of this model, with or without modification, are permitted provided that the following conditions are met:

Redistributions of this SBML file must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in a different form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

This model is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Units

| Time | Substance | Extent | Volume | Area | Length |
| --- | --- | --- | --- | --- | --- |
| min | mmol | mmol | l | m<sup>2</sup> | m |

## Compartments

| Symbol | Id | Name | Size | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{Vbi}$ | `Vbi` | bile |  | l |  |
| $\mathrm{Vre}$ | `Vre` | rest of body |  | l |  |
| $\mathrm{Vgi}$ | `Vgi` | gastrointestinal tract |  | l |  |
| $\mathrm{Vli}$ | `Vli` | liver |  | l |  |
| $\mathrm{Vlu}$ | `Vlu` | lung |  | l |  |
| $\mathrm{Vve}$ | `Vve` | venous plasma |  | l |  |
| $\mathrm{Var}$ | `Var` | arterial plasma |  | l |  |
| $\mathrm{Vpo}$ | `Vpo` | portal plasma |  | l |  |
| $\mathrm{Vhv}$ | `Vhv` | hepatic venous plasma |  | l |  |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Vre}_{\mathrm{tissue}}$ | `Vre_tissue` | re tissue |  | l |  |
| $\mathrm{Vre}_{\mathrm{plasma}}$ | `Vre_plasma` | re plasma |  | l |  |
| $\mathrm{Vgi}_{\mathrm{tissue}}$ | `Vgi_tissue` | gi tissue |  | l |  |
| $\mathrm{Vgi}_{\mathrm{plasma}}$ | `Vgi_plasma` | gi plasma |  | l |  |
| $\mathrm{Vli}_{\mathrm{tissue}}$ | `Vli_tissue` | li tissue |  | l |  |
| $\mathrm{Vli}_{\mathrm{plasma}}$ | `Vli_plasma` | li plasma |  | l |  |
| $\mathrm{Vlu}_{\mathrm{tissue}}$ | `Vlu_tissue` | lu tissue |  | l |  |
| $\mathrm{Vlu}_{\mathrm{plasma}}$ | `Vlu_plasma` | lu plasma |  | l |  |
| $\mathrm{LI}_{\mathrm{\_Vbil}}$ | `LI__Vbil` | bilirubin compartment | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{Cre}_{\mathrm{plasma\_icg}}$ | `Cre_plasma_icg` | icg (rest) | $\mathrm{Vre}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgi}_{\mathrm{plasma\_icg}}$ | `Cgi_plasma_icg` | icg (gastrointestinal tract) | $\mathrm{Vgi}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{plasma\_icg}}$ | `Cli_plasma_icg` | icg (liver) | $\mathrm{Vli}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{plasma\_icg}}$ | `Clu_plasma_icg` | icg (lung) | $\mathrm{Vlu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{icg}}$ | `Cve_icg` | icg (venous plasma) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{icg}}$ | `Car_icg` | icg (arterial plasma) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{icg}}$ | `Cpo_icg` | icg (portal vein) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{icg}}$ | `Chv_icg` | icg (hepatic vein) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Afeces}_{\mathrm{icg}}$ | `Afeces_icg` | icg (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{LI}_{\mathrm{\_icg}}$ | `LI__icg` | icg (liver) | $\mathrm{Vli}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_icg\_bi}}$ | `LI__icg_bi` | icg (bile) | $\mathrm{Vbi}$ | $0$ | mmol/l | concentration |
| $\mathrm{LI}_{\mathrm{\_bil\_ext}}$ | `LI__bil_ext` | bilirubin (extern) | $\mathrm{LI}_{\mathrm{\_Vbil}}$ | $0.01$ | mmol/l | concentration, boundary |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{BW}$ | `BW` | body weight \[kg\] | $75$ | kg | ✓ |
| $\mathrm{HEIGHT}$ | `HEIGHT` | height \[cm\] | $170$ | cm | ✓ |
| $\mathrm{BSA}$ | `BSA` | body surface area \[m^2\] |  | m<sup>2</sup> |  |
| $\mathrm{COBW}$ | `COBW` | cardiac output per bodyweight \[ml/s/kg\] | $1.25$ | ml/s/kg | ✓ |
| $\mathrm{CO}$ | `CO` | cardiac output \[ml/s\] |  | ml/s |  |
| $\mathrm{QC}$ | `QC` | cardiac output \[L/min\] |  | l/min |  |
| $\mathrm{Fblood}$ | `Fblood` | blood fraction of organ volume | $0.02$ | \- |  |
| $\mathrm{HCT}$ | `HCT` | hematocrit | $0.45$ | \- | ✓ |
| $\mathrm{FVgi}$ | `FVgi` | fractional tissue volume gastrointestinal tract | $0.0297$ | l/kg | ✓ |
| $\mathrm{FVbi}$ | `FVbi` | fractional tissue volume bile | $0.00071$ | l/kg | ✓ |
| $\mathrm{FVli}$ | `FVli` | fractional tissue volume liver | $0.021$ | l/kg | ✓ |
| $\mathrm{FVlu}$ | `FVlu` | fractional tissue volume lung | $0.0076$ | l/kg | ✓ |
| $\mathrm{FVve}$ | `FVve` | fractional tissue volume venous | $0.0587$ | l/kg | ✓ |
| $\mathrm{FVar}$ | `FVar` | fractional tissue volume arterial | $0.0184$ | l/kg | ✓ |
| $\mathrm{FVpo}$ | `FVpo` | fractional tissue volume portal | $0.001$ | l/kg | ✓ |
| $\mathrm{FVhv}$ | `FVhv` | fractional tissue volume hepatic venous | $0.001$ | l/kg | ✓ |
| $\mathrm{FVre}$ | `FVre` | fractional tissue volume rest of body |  | l/kg |  |
| $\mathrm{FQgi}$ | `FQgi` | gastrointestinal tract fractional tissue blood flow | $0.19$ | \- | ✓ |
| $\mathrm{FQh}$ | `FQh` | hepatic (venous side) fractional tissue blood flow | $0.255$ | \- | ✓ |
| $\mathrm{FQlu}$ | `FQlu` | lung fractional tissue blood flow | $1$ | \- | ✓ |
| $\mathrm{FQre}$ | `FQre` | rest of body fractional tissue blood flow |  | \- |  |
| $f_{\mathrm{shunts}}$ | `f_shunts` | fraction of portal venous blood shunted by the liver | $0$ | \- | ✓ |
| $f_{\mathrm{tissue\_loss}}$ | `f_tissue_loss` | fraction of lost parenchymal liver volume | $0$ | \- | ✓ |
| $f_{\mathrm{bloodflow}}$ | `f_bloodflow` | fractional change of hepatic blood flow | $1$ | \- | ✓ |
| $f_{\mathrm{cardiac\_output}}$ | `f_cardiac_output` | fractional change of cardiac output | $1$ | \- | ✓ |
| $\mathrm{resection}_{\mathrm{rate}}$ | `resection_rate` | resection rate | $0$ | \- | ✓ |
| $\mathrm{Mr}_{\mathrm{icg}}$ | `Mr_icg` | Molecular weight icg \[g/mole\] | $774.96493$ | g/mol | ✓ |
| $\mathrm{IVDOSE}_{\mathrm{icg}}$ | `IVDOSE_icg` | IV bolus dose icg \[mg\] | $0$ | mg |  |
| $\mathrm{ti}_{\mathrm{icg}}$ | `ti_icg` | injection time icg \[s\] | $5$ | s | ✓ |
| $\mathrm{Ki}_{\mathrm{icg}}$ | `Ki_icg` | Ki \[1/min\] injection icg |  | 1/min |  |
| $\mathrm{Ri}_{\mathrm{icg}}$ | `Ri_icg` | Ri \[mg/min\] rate of infusion icg | $0$ | mg/min | ✓ |
| $\mathrm{Qlu}_{\mathrm{plasma}}$ | `Qlu_plasma` | plasma flow lung |  | l/min |  |
| $\mathrm{Qre}_{\mathrm{plasma}}$ | `Qre_plasma` | plasma flow rest of body |  | l/min |  |
| $\mathrm{Qgi}_{\mathrm{plasma}}$ | `Qgi_plasma` | plasma flow gastrointestinal tract |  | l/min |  |
| $\mathrm{Qh}_{\mathrm{plasma}}$ | `Qh_plasma` | plasma flow hepatic (venous side) |  | l/min |  |
| $\mathrm{Qha}_{\mathrm{plasma}}$ | `Qha_plasma` | plasma flow hepatic artery |  | l/min |  |
| $\mathrm{Qpo}_{\mathrm{plasma}}$ | `Qpo_plasma` | plasma flow portal vein |  | l/min |  |
| $\mathrm{Qlu}$ | `Qlu` | blood flow lung |  | l/min |  |
| $\mathrm{Qre}$ | `Qre` | blood flow rest of body |  | l/min |  |
| $\mathrm{Qh}$ | `Qh` | blood flow hepatic (venous side) |  | l/min |  |
| $\mathrm{Qgi}$ | `Qgi` | blood flow gastrointestinal tract |  | l/min |  |
| $\mathrm{Qpo}$ | `Qpo` | blood flow portal vein |  | l/min |  |
| $\mathrm{Qha}$ | `Qha` | blood flow hepatic artery |  | l/min |  |
| $\mathrm{ER}_{\mathrm{icg}}$ | `ER_icg` | extraction ratio ICG |  | \- |  |
| $\mathrm{CLinfusion}_{\mathrm{icg}}$ | `CLinfusion_icg` | infusion clearance (steady state) |  | l/min |  |
| $\mathrm{LI}_{\mathrm{\_ICGIM\_Vmax}}$ | `LI__ICGIM_Vmax` | Vmax for icg import | $0.0981895740883823$ | mmol/min/l | ✓ |
| $\mathrm{LI}_{\mathrm{\_ICGIM\_Km}}$ | `LI__ICGIM_Km` | Km icg of icg import | $0.025$ | mmol/l | ✓ |
| $\mathrm{LI}_{\mathrm{\_ICGIM\_ki\_bil}}$ | `LI__ICGIM_ki_bil` | Ki bilirubin of icg import | $0.02$ | mmol/l | ✓ |
| $\mathrm{LI}_{\mathrm{\_f\_oatp1b3}}$ | `LI__f_oatp1b3` | scaling factor protein amount | $1$ | \- | ✓ |
| $\mathrm{LI}_{\mathrm{\_ICGEX\_k}}$ | `LI__ICGEX_k` | rate of the icg reflux into plasma | $0.0375437305083155$ | 1/min | ✓ |
| $\mathrm{LI}_{\mathrm{\_ICGLI2CA\_Vmax}}$ | `LI__ICGLI2CA_Vmax` | Vmax bile export icg | $0.000649489174364763$ | mmol/min/l | ✓ |
| $\mathrm{LI}_{\mathrm{\_ICGLI2CA\_Km}}$ | `LI__ICGLI2CA_Km` | Km bile export icg | $0.00638720325187036$ | mmol/l | ✓ |
| $\mathrm{LI}_{\mathrm{\_ICGLI2BI\_k}}$ | `LI__ICGLI2BI_k` | rate of the icg bile transport | $0.00578221865569546$ | 1/min | ✓ |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{Afeces}_{\mathrm{icg}} &= 0 \cdot \mathrm{Vfeces}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{FVre} &= 1 - \mathrm{FVbi} - \mathrm{FVgi} - \mathrm{FVli} \\
&\quad {} - \mathrm{FVlu} - \mathrm{FVve} - \mathrm{FVar} - \mathrm{FVpo} \\
&\quad {} - \mathrm{FVhv} \\[1ex]
\mathrm{BSA} &= 0.024265 \cdot \mathopen{}\left(\frac{\mathrm{BW}}{1}\right)^{0.5378} \cdot \mathopen{}\left(\frac{\mathrm{HEIGHT}}{1}\right)^{0.3964} \\[1ex]
\mathrm{CO} &= \mathrm{BW} \cdot \mathrm{COBW} \cdot f_{\mathrm{cardiac\_output}} \\[1ex]
\mathrm{QC} &= \frac{\mathrm{CO}}{1000} \cdot 60 \\[1ex]
\mathrm{Vbi} &= \mathrm{BW} \cdot \mathrm{FVbi} \\[1ex]
\mathrm{Vgi} &= \mathrm{BW} \cdot \mathrm{FVgi} \\[1ex]
\mathrm{Vli} &= \mathrm{BW} \cdot \mathrm{FVli} \cdot \mathopen{}\left(1 - \mathrm{resection}_{\mathrm{rate}}\right) \\[1ex]
\mathrm{Vlu} &= \mathrm{BW} \cdot \mathrm{FVlu} \\[1ex]
\mathrm{Vre} &= \mathrm{BW} \cdot \mathrm{FVre} \\[1ex]
\mathrm{Vve} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVve} - \frac{\mathrm{FVve}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Var} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVar} - \frac{\mathrm{FVar}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Vpo} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVpo} - \frac{\mathrm{FVpo}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Vhv} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVhv} - \frac{\mathrm{FVhv}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Qlu} &= \mathrm{QC} \cdot \mathrm{FQlu} \\[1ex]
\mathrm{Qh} &= \mathrm{QC} \cdot \mathrm{FQh} \cdot f_{\mathrm{bloodflow}} \\[1ex]
\mathrm{FQre} &= 1 - \frac{\mathrm{Qh}}{\mathrm{Qlu}} \\[1ex]
\mathrm{Qre} &= \mathrm{QC} \cdot \mathrm{FQre} \\[1ex]
\mathrm{Qgi} &= \mathrm{QC} \cdot \mathrm{FQgi} \cdot f_{\mathrm{bloodflow}} \\[1ex]
\mathrm{Qpo} &= \mathrm{Qgi} \\[1ex]
\mathrm{Qha} &= \mathrm{Qh} - \mathrm{Qpo} \\[1ex]
\mathrm{Qlu}_{\mathrm{plasma}} &= \mathrm{Qlu} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Qre}_{\mathrm{plasma}} &= \mathrm{Qre} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Qgi}_{\mathrm{plasma}} &= \mathrm{Qgi} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Qh}_{\mathrm{plasma}} &= \mathrm{Qh} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Qha}_{\mathrm{plasma}} &= \mathrm{Qha} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Qpo}_{\mathrm{plasma}} &= \mathrm{Qpo} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{plasma}} &= \mathrm{Vre} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{tissue}} &= \mathrm{Vre} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vgi}_{\mathrm{plasma}} &= \mathrm{Vgi} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vgi}_{\mathrm{tissue}} &= \mathrm{Vgi} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{plasma}} &= \mathrm{Vli} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{tissue}} &= \mathrm{Vli} \cdot \mathopen{}\left(1 - f_{\mathrm{tissue\_loss}}\right) \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{plasma}} &= \mathrm{Vlu} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{tissue}} &= \mathrm{Vlu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Ki}_{\mathrm{icg}} &= \frac{0.693}{\mathrm{ti}_{\mathrm{icg}}} \cdot 60 \\[1ex]
\mathrm{ER}_{\mathrm{icg}} &= \frac{\mathrm{Car}_{\mathrm{icg}} - \mathrm{Chv}_{\mathrm{icg}}}{\mathrm{Car}_{\mathrm{icg}} + 10^{-12}} \\[1ex]
\mathrm{CLinfusion}_{\mathrm{icg}} &= \frac{\frac{\mathrm{Ri}_{\mathrm{icg}}}{\mathrm{Mr}_{\mathrm{icg}}}}{\mathrm{Cve}_{\mathrm{icg}} + 10^{-12}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation | Modifiers |
| --- | --- | --- | --- | --- |
| $v_{\mathrm{iv\_icg}}$ | `iv_icg` | iv icg | $\varnothing \longrightarrow \mathrm{Cve}_{\mathrm{icg}}$ |  |
| $v_{\mathrm{Flow\_ar\_arre\_icg}}$ | `Flow_ar_arre_icg` | inflow rest icg | $\mathrm{Car}_{\mathrm{icg}} \longrightarrow \mathrm{Cre}_{\mathrm{plasma\_icg}}$ |  |
| $v_{\mathrm{Flow\_re\_ve\_icg}}$ | `Flow_re_ve_icg` | outflow rest icg | $\mathrm{Cre}_{\mathrm{plasma\_icg}} \longrightarrow \mathrm{Cve}_{\mathrm{icg}}$ |  |
| $v_{\mathrm{Flow\_ar\_argi\_icg}}$ | `Flow_ar_argi_icg` | inflow gastrointestinal tract icg | $\mathrm{Car}_{\mathrm{icg}} \longrightarrow \mathrm{Cgi}_{\mathrm{plasma\_icg}}$ |  |
| $v_{\mathrm{Flow\_gi\_po\_icg}}$ | `Flow_gi_po_icg` | outflow gastrointestinal tract icg | $\mathrm{Cgi}_{\mathrm{plasma\_icg}} \longrightarrow \mathrm{Cpo}_{\mathrm{icg}}$ |  |
| $v_{\mathrm{Flow\_arli\_li\_icg}}$ | `Flow_arli_li_icg` | arterial inflow liver icg | $\mathrm{Car}_{\mathrm{icg}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_icg}}$ |  |
| $v_{\mathrm{Flow\_arli\_hv\_icg}}$ | `Flow_arli_hv_icg` | flow arterial shunts | $\mathrm{Car}_{\mathrm{icg}} \longrightarrow \mathrm{Chv}_{\mathrm{icg}}$ |  |
| $v_{\mathrm{Flow\_po\_li\_icg}}$ | `Flow_po_li_icg` | outflow po icg | $\mathrm{Cpo}_{\mathrm{icg}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_icg}}$ |  |
| $v_{\mathrm{Flow\_po\_hv\_icg}}$ | `Flow_po_hv_icg` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{icg}} \longrightarrow \mathrm{Chv}_{\mathrm{icg}}$ |  |
| $v_{\mathrm{Flow\_li\_hv\_icg}}$ | `Flow_li_hv_icg` | outflow liver icg | $\mathrm{Cli}_{\mathrm{plasma\_icg}} \longrightarrow \mathrm{Chv}_{\mathrm{icg}}$ |  |
| $v_{\mathrm{Flow\_hv\_ve\_icg}}$ | `Flow_hv_ve_icg` | outflow hepatic vein icg | $\mathrm{Chv}_{\mathrm{icg}} \longrightarrow \mathrm{Cve}_{\mathrm{icg}}$ |  |
| $v_{\mathrm{Flow\_ve\_lu\_icg}}$ | `Flow_ve_lu_icg` | inflow lung icg | $\mathrm{Cve}_{\mathrm{icg}} \longrightarrow \mathrm{Clu}_{\mathrm{plasma\_icg}}$ |  |
| $v_{\mathrm{Flow\_lu\_ar\_icg}}$ | `Flow_lu_ar_icg` | outflow lung icg | $\mathrm{Clu}_{\mathrm{plasma\_icg}} \longrightarrow \mathrm{Car}_{\mathrm{icg}}$ |  |
| $v_{\mathrm{LI\_\_ICGIM}}$ | `LI__ICGIM` | ICG import | $\mathrm{Cli}_{\mathrm{plasma\_icg}} \longrightarrow \mathrm{LI}_{\mathrm{\_icg}}$ | $\mathrm{LI}_{\mathrm{\_bil\_ext}}$ |
| $v_{\mathrm{LI\_\_ICGEX}}$ | `LI__ICGEX` | ICG export (plasma) | $\mathrm{LI}_{\mathrm{\_icg}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_icg}}$ |  |
| $v_{\mathrm{LI\_\_ICGLI2CA}}$ | `LI__ICGLI2CA` | ICG export (bile) | $\mathrm{LI}_{\mathrm{\_icg}} \longrightarrow \mathrm{LI}_{\mathrm{\_icg\_bi}}$ |  |
| $v_{\mathrm{LI\_\_ICGLI2BI}}$ | `LI__ICGLI2BI` | ICG bile transport | $\mathrm{LI}_{\mathrm{\_icg\_bi}} \longrightarrow \mathrm{Afeces}_{\mathrm{icg}}$ |  |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{iv\_icg}} &= \frac{\mathrm{Ki}_{\mathrm{icg}} \cdot \mathrm{IVDOSE}_{\mathrm{icg}}}{\mathrm{Mr}_{\mathrm{icg}}} \\[1ex]
v_{\mathrm{Flow\_ar\_arre\_icg}} &= \mathrm{Qre}_{\mathrm{plasma}} \cdot \mathrm{Car}_{\mathrm{icg}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_icg}} &= \mathrm{Qre}_{\mathrm{plasma}} \cdot \mathrm{Cre}_{\mathrm{plasma\_icg}} \\[1ex]
v_{\mathrm{Flow\_ar\_argi\_icg}} &= \mathrm{Qgi}_{\mathrm{plasma}} \cdot \mathrm{Car}_{\mathrm{icg}} \\[1ex]
v_{\mathrm{Flow\_gi\_po\_icg}} &= \mathrm{Qgi}_{\mathrm{plasma}} \cdot \mathrm{Cgi}_{\mathrm{plasma\_icg}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_icg}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha}_{\mathrm{plasma}} \cdot \mathrm{Car}_{\mathrm{icg}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_icg}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha}_{\mathrm{plasma}} \cdot \mathrm{Car}_{\mathrm{icg}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_icg}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo}_{\mathrm{plasma}} \cdot \mathrm{Cpo}_{\mathrm{icg}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_icg}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo}_{\mathrm{plasma}} \cdot \mathrm{Cpo}_{\mathrm{icg}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_icg}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo}_{\mathrm{plasma}} + \mathrm{Qha}_{\mathrm{plasma}}\right) \cdot \mathrm{Cli}_{\mathrm{plasma\_icg}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_icg}} &= \mathrm{Qh}_{\mathrm{plasma}} \cdot \mathrm{Chv}_{\mathrm{icg}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_icg}} &= \mathrm{Qlu}_{\mathrm{plasma}} \cdot \mathrm{Cve}_{\mathrm{icg}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_icg}} &= \mathrm{Qlu}_{\mathrm{plasma}} \cdot \mathrm{Clu}_{\mathrm{plasma\_icg}} \\[1ex]
v_{\mathrm{LI\_\_ICGIM}} &= \frac{\mathrm{LI}_{\mathrm{\_f\_oatp1b3}} \cdot \mathrm{LI}_{\mathrm{\_ICGIM\_Vmax}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{Cli}_{\mathrm{plasma\_icg}}}{\mathrm{LI}_{\mathrm{\_ICGIM\_Km}} \cdot \mathopen{}\left(1 + \frac{\mathrm{LI}_{\mathrm{\_bil\_ext}}}{\mathrm{LI}_{\mathrm{\_ICGIM\_ki\_bil}}}\right) + \mathrm{Cli}_{\mathrm{plasma\_icg}}} \\[1ex]
v_{\mathrm{LI\_\_ICGEX}} &= \mathrm{LI}_{\mathrm{\_ICGEX\_k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_icg}} \\[1ex]
v_{\mathrm{LI\_\_ICGLI2CA}} &= \frac{\mathrm{LI}_{\mathrm{\_ICGLI2CA\_Vmax}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{LI}_{\mathrm{\_icg}}}{\mathrm{LI}_{\mathrm{\_ICGLI2CA\_Km}} + \mathrm{LI}_{\mathrm{\_icg}}} \\[1ex]
v_{\mathrm{LI\_\_ICGLI2BI}} &= \mathrm{LI}_{\mathrm{\_ICGLI2BI\_k}} \cdot \mathrm{Vbi} \cdot \mathrm{LI}_{\mathrm{\_icg\_bi}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{plasma\_icg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_arre\_icg}} - v_{\mathrm{Flow\_re\_ve\_icg}}}{\mathrm{Vre}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgi}_{\mathrm{plasma\_icg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_argi\_icg}} - v_{\mathrm{Flow\_gi\_po\_icg}}}{\mathrm{Vgi}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{plasma\_icg}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Vli}_{\mathrm{plasma}}} \Bigl( v_{\mathrm{Flow\_arli\_li\_icg}} + v_{\mathrm{Flow\_po\_li\_icg}} - v_{\mathrm{Flow\_li\_hv\_icg}} - v_{\mathrm{LI\_\_ICGIM}} \\
&\quad {} + v_{\mathrm{LI\_\_ICGEX}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{plasma\_icg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ve\_lu\_icg}} - v_{\mathrm{Flow\_lu\_ar\_icg}}}{\mathrm{Vlu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{icg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{iv\_icg}} + v_{\mathrm{Flow\_re\_ve\_icg}} + v_{\mathrm{Flow\_hv\_ve\_icg}} - v_{\mathrm{Flow\_ve\_lu\_icg}}}{\mathrm{Vve}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{icg}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_arre\_icg}} - v_{\mathrm{Flow\_ar\_argi\_icg}} - v_{\mathrm{Flow\_arli\_li\_icg}} - v_{\mathrm{Flow\_arli\_hv\_icg}} \\
&\quad {} + v_{\mathrm{Flow\_lu\_ar\_icg}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{icg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_gi\_po\_icg}} - v_{\mathrm{Flow\_po\_li\_icg}} - v_{\mathrm{Flow\_po\_hv\_icg}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{icg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_icg}} + v_{\mathrm{Flow\_po\_hv\_icg}} + v_{\mathrm{Flow\_li\_hv\_icg}} - v_{\mathrm{Flow\_hv\_ve\_icg}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{icg}}}{\mathrm{d} t} &= v_{\mathrm{LI\_\_ICGLI2BI}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_icg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_ICGIM}} - v_{\mathrm{LI\_\_ICGEX}} - v_{\mathrm{LI\_\_ICGLI2CA}}}{\mathrm{Vli}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{LI}_{\mathrm{\_icg\_bi}}}{\mathrm{d} t} &= \frac{v_{\mathrm{LI\_\_ICGLI2CA}} - v_{\mathrm{LI\_\_ICGLI2BI}}}{\mathrm{Vbi}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{IVDOSE}_{\mathrm{icg}}}{\mathrm{d} t} &= -v_{\mathrm{iv\_icg}} \cdot \mathrm{Mr}_{\mathrm{icg}} + \mathrm{Ri}_{\mathrm{icg}} \qquad \text{(rate rule)}
\end{aligned}
$$
