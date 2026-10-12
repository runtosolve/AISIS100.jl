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

- Covers Chapters **D–H (member design)**, **J (connections)**, **K (testing)**, **Appendix 1 (effective width of elements)**, and **Appendix 2 (elastic buckling, distortional buckling)**
- `v16S3` also implements the **Effective Width Method** for flexure (§F3.1 with Appendix 1, including C-section webs with holes, §1.1.3), the **strength increase from cold work of forming** (§A3.3.2), and the closed-form distortional buckling moment (Appendix 2 §2.3.3.3, §2.3.4.3 for members with web holes)
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

**Effective Width Method, lipped C-section in major-axis bending — §F3.1 with Appendix 1 (ASD)**

```julia
# 362S162-33 stud, Fy = 33 ksi: effective section at f = Fy, then Eq. F3.1-1
sec = AISIS100.v16S3.f311_lipped_cee(H=3.625, B=1.625, D=0.5, t=0.0346, R=0.0765,
                                     f=33.0, E=29500.0, μ=0.3)
Mn, Ma = AISIS100.v16S3.f31(Se=sec.Se, Fn=33.0, Set=sec.Set, Fy=33.0, design_code="ASD")

# same section through a 1.5 in. web punchout (§F3.1.2, Appendix 1 §1.1.3)
net = AISIS100.v16S3.f311_lipped_cee(H=3.625, B=1.625, D=0.5, t=0.0346, R=0.0765,
                                     f=33.0, E=29500.0, μ=0.3, dh=1.5, hole_method=:appendix1_1_2_2)

# cold work of forming, §A3.3.2 (applies only to sections fully effective at Fy)
C = AISIS100.v16S3.a332_C(w=sec.w, t=0.0346, R=0.0765, n_corners=2)
Fya = AISIS100.v16S3.a332(Fyv=33.0, Fuv=45.0, R=0.0765, t=0.0346, C=C).Fya
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
