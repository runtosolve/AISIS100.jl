# AISIS100.jl

*AISI S100 cold-formed steel design equations, implemented directly in Julia.*

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia](https://img.shields.io/badge/Julia-1.6%2B-blueviolet?logo=julia)](https://julialang.org)
[![AISI S100](https://img.shields.io/badge/AISI-S100--16%20%7C%20S100--24-blue)](https://www.steel.org/steel-markets/specialty-steel/cold-formed-steel/)

---

## The Problem

Cold-formed steel design requires working through dozens of AISI S100 equations — by hand, in spreadsheets, or inside proprietary software that cannot be inspected or extended. This leads to duplicated effort, limited transparency, and increased risk of error.

---

## What It Does

`AISIS100.jl` implements the **AISI S100 North American Specification for the Design of Cold-Formed Steel Structural Members** as callable Julia functions.

- Covers Chapters **D–H (member design)**, **J (connections)**, **K (testing)**, and **Appendix 2 (distortional buckling)**
- Returns both **nominal strength** and **design strength** (ASD, LRFD, LSD)
- Function names map directly to **spec section numbers**

Three code editions are supported:

| Module  | Standard                       |
|---------|--------------------------------|
| `v16`   | AISI S100-16                   |
| `v16S3` | AISI S100-16 with Supplement 3 |
| `v2024` | AISI S100-24                   |

---

## Prerequisites

- [Julia 1.6 or later](https://julialang.org/downloads/)

---

## Installation

Open the Julia REPL, press `]` to enter package mode, then run:

```julia
pkg> add https://github.com/runtosolve/AISIS100.jl.git
```

---

## Quick Start

Function names correspond directly to AISI S100 section numbers. Open the spec, find the equation, call the function.

**Distortional buckling of a stud column — §E4.1 (LRFD)**

```julia
using AISIS100

Pnd, ϕPnd = AISIS100.v16.e41(Py=48.42, Pcrd=13.1, design_code="AISI S100-16 LRFD")

println("Nominal: $Pnd kips")
println("Design:  $ϕPnd kips")
```

```
Nominal: 24.17 kips
Design:  20.54 kips
```

**Lateral-torsional buckling of a beam — §F2.1 (ASD)**

```julia
Mne, eMne = AISIS100.v16.f21(Fcre=300.0, Fy=55.0, Sf=10.0, Sfy=12.6,
                               design_code="AISI S100-16 ASD")

println("Nominal: $Mne kip·in")
println("Allowable: $eMne kip·in")
```

```
Nominal:   550.0 kip·in
Allowable: 329.0 kip·in
```

**Combined axial and biaxial bending — §H1.2.1**

```julia
action_P, action_Mx, action_My, interaction = AISIS100.v16.h121(
    Pbar=10.0, Mxbar=6.0, Mybar=3.0,
    Pa=20.0, Max=40.0, May=35.0)

println("Interaction ratio: $interaction")
```

```
Interaction ratio: 0.71   # ≤ 1.0 → passes
```

---

## Contributing

Found an equation that is missing or incorrectly implemented? Open an issue with the relevant AISI S100 section number and edition. Pull requests are welcome.

---

## License

Released under the [MIT License](LICENSE).

---

## Author

Developed by **Cris Moen** — [runtosolve.com](https://runtosolve.com)
