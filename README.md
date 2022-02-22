S100AISI.jl
==========

*American Iron and Steel Institute (AISI) S100 North American Specification for the Design of Cold-Formed Steel Structural Members*

Use design equations from AISI S100 in your calculations.

Install
-----------------------------

```julia
(v1.7) pkg> add S100AISI
```

(Type `]` to enter package mode.)

Example Usage
-------------

```julia
using S100AISI
```

Calculate the distortional buckling strength of a cold-formed steel stud column using AISI S100-16 Eq. E4.1. 

```julia
 Pnd, ϕPnd = v16.e41(Py=48.42, Pcrd=13.1, design_code = "AISI S100-16 LRFD")
```

