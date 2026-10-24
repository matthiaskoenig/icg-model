# model for icg elimination by the liver

Model `icg_liver`, SBML Level 3 Version 2, read from icg\_liver.xml, written by sbmlode 0.3.0.

Model for icg elimination by the liver

Description

Model for hepatic elimination of indocyanine green encoded in SBML format.

Assumptions:

indocyanine green is not metabolized; it is taken up by the liver (OATP1B3) and excreted unchanged into bile (MRP2)

a part of the indocyanine green in the liver returns into the sinusoidal plasma (reflux), which is the slow secondary phase of the plasma curve

bilirubin inhibits the uptake competitively

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
| $\mathrm{Vext}$ | `Vext` | plasma | $4$ | l | ✓ |
| $\mathrm{Vli}$ | `Vli` | liver | $1.5$ | l | ✓ |
| $\mathrm{Vbi}$ | `Vbi` | bile | $0.0507142857142857$ | l | ✓ |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Vbil}$ | `Vbil` | bilirubin compartment | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{icg}_{\mathrm{ext}}$ | `icg_ext` | icg (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{icg}$ | `icg` | icg (liver) | $\mathrm{Vli}$ | $0$ | mmol/l | concentration |
| $\mathrm{icg}_{\mathrm{bi}}$ | `icg_bi` | icg (bile) | $\mathrm{Vbi}$ | $0$ | mmol/l | concentration |
| $\mathrm{icg}_{\mathrm{feces}}$ | `icg_feces` | icg (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{bil}_{\mathrm{ext}}$ | `bil_ext` | bilirubin (extern) | $\mathrm{Vbil}$ | $0.01$ | mmol/l | concentration, boundary |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{ICGIM}_{\mathrm{Vmax}}$ | `ICGIM_Vmax` | Vmax for icg import | $0.0981895740883823$ | mmol/min/l | ✓ |
| $\mathrm{ICGIM}_{\mathrm{Km}}$ | `ICGIM_Km` | Km icg of icg import | $0.025$ | mmol/l | ✓ |
| $\mathrm{ICGIM}_{\mathrm{ki\_bil}}$ | `ICGIM_ki_bil` | Ki bilirubin of icg import | $0.02$ | mmol/l | ✓ |
| $f_{\mathrm{oatp1b3}}$ | `f_oatp1b3` | scaling factor protein amount | $1$ | \- | ✓ |
| $\mathrm{ICGEX}_{\mathrm{k}}$ | `ICGEX_k` | rate of the icg reflux into plasma | $0.0375437305083155$ | 1/min | ✓ |
| $\mathrm{ICGLI2CA}_{\mathrm{Vmax}}$ | `ICGLI2CA_Vmax` | Vmax bile export icg | $0.000649489174364763$ | mmol/min/l | ✓ |
| $\mathrm{ICGLI2CA}_{\mathrm{Km}}$ | `ICGLI2CA_Km` | Km bile export icg | $0.00638720325187036$ | mmol/l | ✓ |
| $\mathrm{ICGLI2BI}_{\mathrm{k}}$ | `ICGLI2BI_k` | rate of the icg bile transport | $0.00578221865569546$ | 1/min | ✓ |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{icg}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation | Modifiers |
| --- | --- | --- | --- | --- |
| $v_{\mathrm{ICGIM}}$ | `ICGIM` | ICG import | $\mathrm{icg}_{\mathrm{ext}} \longrightarrow \mathrm{icg}$ | $\mathrm{bil}_{\mathrm{ext}}$ |
| $v_{\mathrm{ICGEX}}$ | `ICGEX` | ICG export (plasma) | $\mathrm{icg} \longrightarrow \mathrm{icg}_{\mathrm{ext}}$ |  |
| $v_{\mathrm{ICGLI2CA}}$ | `ICGLI2CA` | ICG export (bile) | $\mathrm{icg} \longrightarrow \mathrm{icg}_{\mathrm{bi}}$ |  |
| $v_{\mathrm{ICGLI2BI}}$ | `ICGLI2BI` | ICG bile transport | $\mathrm{icg}_{\mathrm{bi}} \longrightarrow \mathrm{icg}_{\mathrm{feces}}$ |  |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{ICGIM}} &= \frac{f_{\mathrm{oatp1b3}} \cdot \mathrm{ICGIM}_{\mathrm{Vmax}} \cdot \mathrm{Vli} \cdot \mathrm{icg}_{\mathrm{ext}}}{\mathrm{ICGIM}_{\mathrm{Km}} \cdot \mathopen{}\left(1 + \frac{\mathrm{bil}_{\mathrm{ext}}}{\mathrm{ICGIM}_{\mathrm{ki\_bil}}}\right) + \mathrm{icg}_{\mathrm{ext}}} \\[1ex]
v_{\mathrm{ICGEX}} &= \mathrm{ICGEX}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathrm{icg} \\[1ex]
v_{\mathrm{ICGLI2CA}} &= \frac{\mathrm{ICGLI2CA}_{\mathrm{Vmax}} \cdot \mathrm{Vli} \cdot \mathrm{icg}}{\mathrm{ICGLI2CA}_{\mathrm{Km}} + \mathrm{icg}} \\[1ex]
v_{\mathrm{ICGLI2BI}} &= \mathrm{ICGLI2BI}_{\mathrm{k}} \cdot \mathrm{Vbi} \cdot \mathrm{icg}_{\mathrm{bi}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{icg}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{ICGIM}} + v_{\mathrm{ICGEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{icg}}{\mathrm{d} t} &= \frac{v_{\mathrm{ICGIM}} - v_{\mathrm{ICGEX}} - v_{\mathrm{ICGLI2CA}}}{\mathrm{Vli}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{icg}_{\mathrm{bi}}}{\mathrm{d} t} &= \frac{v_{\mathrm{ICGLI2CA}} - v_{\mathrm{ICGLI2BI}}}{\mathrm{Vbi}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{icg}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{ICGLI2BI}}
\end{aligned}
$$
