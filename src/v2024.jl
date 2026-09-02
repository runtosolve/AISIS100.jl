module v2024

# Function naming convention
#
# Examples:
#   Eq. E4-4          -> e4__4
#   Eq. E2.1-2         -> e2_1__2      (dot in "E2.1" -> underscore, dash -> double underscore)
#   Eqs. G3-1 to G3-5  -> g3__1to5
#   Eq. 2.3.1.1.4-3    -> appendix2_2_3_1_1_4__3
#
using LinearAlgebra
using Statistics

#General helper (Sections B3.2.1-B3.2.3): converts a nominal strength Rn to an available strength (ASD, LRFD, or LSD)
function calculate_factored_strength(Rn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    if design_code == "ASD"
        eRn  = Rn / Ω
    elseif design_code == "LRFD"
        eRn = Rn * ϕ_LRFD
    elseif design_code == "LSD"
        eRn = Rn * ϕ_LSD
    elseif design_code == "nominal"
        eRn = Rn
    end

    return eRn

end

#---------------------------------------------------------------------------
# CHAPTER E — MEMBERS IN COMPRESSION
#---------------------------------------------------------------------------

#Eqs. E2-1 to E2-4: nominal axial strength for yielding and global (flexural, flexural-torsional, or torsional) buckling (Pne, Fn, λc)
function e2__1to4(Fcre, Fy, Ag, design_code)

    Ω = 1.80
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80

    λc = sqrt(Fy/Fcre)

    if λc <= 1.5
        Fn = 0.658^((λc)^2) * Fy
    else
        Fn = (0.877 /(λc)^2) * Fy
    end

    Pne = Ag * Fn

    ePne = calculate_factored_strength(Pne, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    results = (Rn=Pne, eRn=ePne)

    return results

end


#Eqs. E2.1-1 and E2.1-2: closed-box section reduction to Fcre — length at which local buckling stress equals flexural buckling stress (L0), and the reduction factor (R)
function e2_1__1to2(r, E, Fcr, K, L)

    L0 = π * r * sqrt(E / Fcr)

    R = 0.65 + 0.35 * (K * L / (1.1 * L0))^2

    return L0, R

end

#Eqs. E3.2-1 and E3.2-2: nominal axial strength for local buckling, Direct Strength Method, no holes (Pnℓ, λℓ)
function e3_2__1to2(Pne, Pcrℓ, design_code)

    Ω = 1.80
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80

    λℓ = sqrt(Pne/Pcrℓ)

    Pnℓ = min(1.2*Pne*((1 + 0.10*λℓ^2)/(1 + 0.55*λℓ^2)), Pne)

    ePnℓ = calculate_factored_strength(Pnℓ, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    results = (Rn=Pnℓ, eRn=ePnℓ)

    return results

end


#Eqs. E3.2-3 and E3.2-4: nominal axial strength for local buckling at a net section with holes (Pynet, Pnℓnet)
function e3_2__3to4(Anet, Fy, Pne, Pcrℓnet, design_code)

    Ω = 1.80
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80

    Pynet = Anet * Fy

    λℓ = sqrt(Pne / Pcrℓnet)

    Pnℓnet = min(1.2 * Pynet * ((1 + 0.10 * λℓ^2) / (1 + 0.55 * λℓ^2 * (Pynet / Pne))), Pynet)

    ePnℓnet = calculate_factored_strength(Pnℓnet, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Pynet, Pnℓnet, ePnℓnet

end


#Eqs. E3.3-1, E3.3-2, and E3.3-3: effective area (Ae) for local buckling of closed cylindrical tubes in compression, D/t ≤ 0.441 E/Fy
function e3_3__1to3(D, Fy, t, E, A, Fcre)

    Ao = min((0.037 / ((D * Fy) / (t * E)) + 0.667) * A, A)

    R = min(Fy / (2 * Fcre), 1.0)

    Ae = Ao + R * (A - Ao)

    return Ao, R, Ae

end

#Eqs. E4-1 to E4-3: nominal axial strength for distortional buckling, no holes (Pnd, λd, Py)
function e4__1to3(Py, Pcrd, design_code)

    Ω = 1.80
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80

    λd = sqrt(Py/Pcrd)

    Pnd = min(1.2*Py*((1 + 0.05*λd^2)/(1 + 0.67*λd^2)), Py)

    ePnd = calculate_factored_strength(Pnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    results = (Rn=Pnd, eRn=ePnd)

    return results

end


#Eq. E4-4: nominal axial strength for distortional buckling at a net section with holes (Pndnet)
function e4__4(Py, Anet, Fy, Pcrdnet, design_code)

    Ω = 1.80
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80

    Pynet = Anet * Fy

    λd = sqrt(Py / Pcrdnet)

    Pndnet = min(1.2 * Pynet * ((1 + 0.05 * λd^2) / (1 + 0.67 * λd^2 * (Pynet / Py))), Pynet)

    ePndnet = calculate_factored_strength(Pndnet, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Pynet, Pndnet, ePndnet

end

#---------------------------------------------------------------------------
# CHAPTER F — MEMBERS IN FLEXURE
#---------------------------------------------------------------------------

#Eqs. F2.1-1 to F2.1-5: nominal flexural strength for yielding and global (lateral-torsional) buckling, Effective Width Method (Mne, My, Fn)
function f2_1__1to5(Fcre, Fy, Sfc, Sf, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    My = Sf * Fy

    if Fcre >= 2.78 * Fy
        Fn = Fy
    elseif Fcre > 0.56 * Fy
        Fn = (10 / 9) * Fy * (1 - 10 * Fy / (36 * Fcre))
    else
        Fn = Fcre
    end

    Mne = min(Sfc * Fn, My)

    eMne = calculate_factored_strength(Mne, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return My, Fn, Mne, eMne

end


#Eq. F2.1.1-1: yield strain (εy), used by the F2.1.1 inelastic reserve strength provisions
function f2_1_1__1(Fy, E)

    εy = Fy / E

    return εy

end


#Eqs. F2.1.1-2 to F2.1.1-4: compression strain factor (Cy) for stiffened compression elements without intermediate stiffeners
function f2_1_1__2to4(w, t, Fy, E)

    λ1 = 1.11 / sqrt(Fy / E)
    λ2 = 1.28 / sqrt(Fy / E)

    if w / t <= λ1
        Cy = 3.0
    elseif w / t < λ2
        Cy = 3 - 2 * (w / t - λ1) / (λ2 - λ1)
    else
        Cy = 1.0
    end

    return λ1, λ2, Cy

end


#Eqs. F2.1.1-5 and F2.1.1-6: compression strain factor (Cy) for unstiffened compression elements under a stress gradient (compression at one edge, tension at the other)
function f2_1_1__5to6(λ, ψ)

    λ3 = 0.43
    λ4 = 0.673 * (1 + ψ)

    if λ <= λ3
        Cy = 3.0
    elseif λ < λ4
        Cy = 3 - 2 * (λ - λ3) / (λ4 - λ3)
    else
        Cy = 1.0
    end

    return λ3, λ4, Cy

end


#Eqs. F2.2-1 to F2.2-3: nominal flexural strength for yielding and global (lateral-torsional) buckling, Direct Strength Method (Mne, Mp)
function f2_2__1to3(Mcre, Zf, Fy, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    Mp = Zf * Fy

    if Mcre > 0.5 * Mp
        Mne = Mp * (1 - 0.25 * Mp / Mcre)
    else
        Mne = Mcre
    end

    eMne = calculate_factored_strength(Mne, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mp, Mne, eMne

end


#Eq. F2.3-1: nominal flexural strength for yielding and global buckling of closed cylindrical tubes, D/t ≤ 0.441 E/Fy
function f2_3__1(Sf, Fy, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.95
    ϕ_LSD = 0.90

    Mne = 1.25 * Sf * Fy

    eMne = calculate_factored_strength(Mne, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mne, eMne

end

#Eqs. F3.2-1 and F3.2-2: nominal flexural strength for local buckling, Direct Strength Method, no holes (Mnℓ, λℓ) — ks, αs, βs, My3 are computed per Eqs. F3.2-3 to F3.2-5 and passed in as arguments
function f3_2__1to2(Mne, Mcrℓ, ks, αs, βs, My3, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λℓ=sqrt(Mne/Mcrℓ)


    Mnℓ = min(ks*Mne*((1+0.10*αs*λℓ*λℓ)/(1+0.55*βs*λℓ*λℓ)),My3)


    # eMnℓ = Mnℓ * StrengthFactor

    eMnℓ = calculate_factored_strength(Mnℓ, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnℓ, eMnℓ

end


#Eqs. F3.2-6 to F3.2-9: nominal flexural strength for local buckling at a net section with holes (Mpnet, My3net, Mynet, Mnℓnet)
function f3_2__6to9(Zfnet, Sfnet, Fy, ks, αs, βs, Mne, Mcrℓnet, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    Mpnet = Zfnet * Fy
    Mynet = Sfnet * Fy
    My3net = Mpnet - (Mpnet - Mynet) / 9

    λℓ = sqrt(Mne / Mcrℓnet)

    Mnℓnet = min(Mpnet * ((1 + 0.10 * αs * λℓ^2) / (1 + 0.55 * βs * λℓ^2 * (Mpnet / (ks * Mne)))), My3net)

    eMnℓnet = calculate_factored_strength(Mnℓnet, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mpnet, Mynet, My3net, Mnℓnet, eMnℓnet

end


#Eqs. F3.3-1 to F3.3-4: nominal flexural strength for local buckling of closed cylindrical tubes (Mnℓ, Mnℓe), D/t ≤ 0.441 E/Fy — safety/resistance factors per Section F2.3
function f3_3__1to4(Mne, t, E, I, D, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.95
    ϕ_LSD = 0.90

    Mnℓe = 0.656 * t * E * I / D^2

    if Mnℓe >= 3.67 * Mne
        Mnℓ = Mne
    elseif Mnℓe >= 0.826 * Mne
        Mnℓ = 0.776 * Mne + 0.061 * Mnℓe
    else
        Mnℓ = Mnℓe
    end

    eMnℓ = calculate_factored_strength(Mnℓ, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnℓe, Mnℓ, eMnℓ

end


#Eqs. F4-1 and F4-2: nominal flexural strength for distortional buckling, no holes (Mnd, λd) — αs, βs, My3 are computed per Eq. F4-3 and Section F3.2, passed in as arguments
function f4__1to2(My, Mcrd, Mp, αs, βs, My3, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λd = sqrt(My/Mcrd)

    Mnd = min(Mp*((1 + 0.07*αs*λd^2)/(1 + 0.60*βs*λd^2)), My3)

    eMnd = calculate_factored_strength(Mnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnd, eMnd

end


#Eq. F4-4: nominal flexural strength for distortional buckling at a net section with holes (Mndnet), using Mpnet and My3net from Section F3.2 (Eqs. F3.2-7/F3.2-8)
function f4__4(My, Mcrdnet, Mp, αs, βs, Mpnet, My3net, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λd = sqrt(My / Mcrdnet)

    Mndnet = min(Mpnet * ((1 + 0.07 * αs * λd^2) / (1 + 0.60 * βs * λd^2 * (Mpnet / Mp))), My3net)

    eMndnet = calculate_factored_strength(Mndnet, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mndnet, eMndnet

end

#---------------------------------------------------------------------------
# CHAPTER G — MEMBERS IN SHEAR, WEB CRIPPLING, AND TORSION
#---------------------------------------------------------------------------

#Eqs. G2.1-1 and G2.1-2: nominal shear strength of webs without transverse stiffeners (Vn, λv) — Vy (Eq. G2.1-3) is computed elsewhere and passed in
function g2_1__1to3(Vcr, Vy, design_code)  #no transverse stiffeners

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.75  #Section G2, AISI S100-2024

    λv=sqrt(Vy/Vcr)


    Vn = min(1.2*Vy/(1+0.57*λv*λv),Vy)

    eVn = calculate_factored_strength(Vn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Vn, eVn

end


#Eq. G2.2-1: nominal shear strength of webs with transverse stiffeners (Vn), in lieu of Section G2.1
function g2_2__1(Vcr, Vy, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.75

    λv = sqrt(Vy / Vcr)

    Vn = min(1.2 * Vy * ((1 + 0.05 * λv^2) / (1 + 0.45 * λv^2)), Vy)

    eVn = calculate_factored_strength(Vn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Vn, eVn

end


#Eqs. G2.3-1 to G2.3-3: web elastic critical shear buckling force (Vcr, Aw, Fcr)
function g2_3__1to3(h, t, kv, E, mu)

    Aw = h * t

    Fcr = π^2 * E * kv / (12 * (1 - mu^2) * (h / t)^2)

    Vcr = Aw * Fcr

    return Aw, Fcr, Vcr

end


#Eqs. G2.3-4 and G2.3-5: shear buckling coefficient (kv) for webs with transverse stiffeners meeting Section G4 — for unreinforced webs, kv = 5.34 per Section G2.3(a)
function g2_3__4to5(a, h)

    if a / h <= 1.0
        kv = 4.00 + 5.34 / (a / h)^2
    else
        kv = 5.34 + 4.00 / (a / h)^2
    end

    return kv

end

#Eqs. G3-1 to G3-5
function g3__1to5(d_h, L_h, h, Vy)

    a0 = -0.173 - 0.9252*(L_h/d_h) + 0.0524*(L_h/d_h)^2
    a1 = -3.4095 + 1.9922*(L_h/d_h) - 0.0995*(L_h/d_h)^2
    a2 = 2.684 - 1.084*(L_h/d_h) + 0.0466*(L_h/d_h)^2

    if (d_h/h<=0.10) & (d_h/h > 0.0)

        V_yh = Vy

    elseif (d_h/h>0.10) & (d_h/h <= 0.80)

        V_yh = Vy * (1 + a0*(d_h/h - 0.1) + a1*(d_h/h-0.1)^2 + a2*(d_h/h-0.1)^3)

    end

    return V_yh, a0, a1, a2

end

#Eqs. G3-6 to G3-7
function g3__6to7(L_h, d_h, h, Vcr)

    α_vh = (1-0.4*(L_h-d_h)/h)^2

    V_crh = α_vh * Vcr

    return α_vh, V_crh

end

#Eq. G3-8
function g3__8(a, h, d_h, b_f)

    kv = 4.86 + 6.15 * (h / a) - 3.63 * (d_h / h) - 19.58 * (d_h / a) + 13.88 * (d_h^2 / (a * h)) + 0.57* (b_f / h)

end


#Eqs. G3-9 to G3-10
function g3__9to10(d_h, L_h)

    A_h = π*(d_h/2)^2  #for circular hole only, not slotted hole
    d_h_eq = (0.003*(L_h/d_h) + 0.822) * d_h
    L_h_eq = 0.865*A_h/d_h_eq

    return A_h, d_h_eq, L_h_eq

end

#Eqs. G4.1-1 to G4.1-3: compact transverse web stiffener sizing checks (bst/tst ≤ limit, ρst, Ist1)
function g4_1__1to3(bst, tst, Fyst, E, h, Fy)

    bst_tst_limit = 0.56 * sqrt(E / Fyst)
    bst_tst_ok = (bst / tst) <= bst_tst_limit

    ρst = max(Fy / Fyst, 1.0)

    Ist1 = (h^4 * ρst^1.3 / 40) * (Fy / E)^1.5

    return bst_tst_ok, bst_tst_limit, ρst, Ist1

end


#Eq. G5-1: nominal web crippling strength of webs without holes (Pn)
function g5__1(t, h, Fy, θ, C, C_R, R, C_N, N, C_h, ϕ_w, Ω_w, ϕ_w_LSD, design_code)

    Pn = C * t^2 * Fy * sin(deg2rad(θ)) * (1-C_R*sqrt(R/t)) * (1+C_N*sqrt(N/t)) * (1-C_h*sqrt(h/t))

    ePn = calculate_factored_strength(Pn, Ω_w, ϕ_w, ϕ_w_LSD, design_code)

    return Pn, ePn

end


#Eqs. G5-2 and G5-3: nominal web crippling strength with an overhang, end one-flange loading, C- and Z-sections (Pnc, α)
function g5__2to3(Pn, Lo, h, t)

    α = max(1.34 * (Lo / h)^0.26 / (0.009 * (h / t) + 0.3), 1.0)

    Pnc = α * Pn

    return α, Pnc

end


#Eq. G6-1: web crippling hole reduction factor (Rc), end one-flange reaction, N ≥ 1 in. (25.4 mm)
function g6__1(dh, h, x)

    Rc = min(1.01 - 0.325 * (dh / h) + 0.083 * (x / h), 1.0)

    return Rc

end


#Eq. G6-2: web crippling hole reduction factor (Rc), interior one-flange reaction, N ≥ 3 in. (76.2 mm)
function g6__2(dh, h, x)

    Rc = min(0.90 - 0.047 * (dh / h) + 0.053 * (x / h), 1.0)

    return Rc

end


#Eqs. G7.1-1 to G7.1-3: nominal strength of compact bearing stiffeners (Pn, Ac) — path (a) only; path (b) requires Section E3.1 (Effective Width Method) and is not implemented
function g7_1__1to3(Fwy, t, As, at_end::Bool, design_code)

    Ω = 2.00
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80

    Ac = at_end ? 10 * t^2 + As : 18 * t^2 + As

    Pn = Fwy * Ac

    ePn = calculate_factored_strength(Pn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Ac, Pn, ePn

end


#Eq. G7.2-1: nominal strength of stud and track type bearing stiffeners in C-section flexural members (Pn)
function g7_2__1(Pwc, Pno, design_code)

    Ω = 1.70
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.80

    Pn = max(0.7 * (Pwc + Pno), Pwc)

    ePn = calculate_factored_strength(Pn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Pn, ePn

end

#Eq. G8.1-2: bimoment yield strength (Bn=By), the simple (non-buckling-reserve) path of Section G8.1 — this provision moved from Section H4.1.1 to Section G8.1 in the 2024 edition, function name kept for backward compatibility with dependent packages
function h411(Cw, Fy, Wn, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    #calculate bimoment yield strength
    Bn=Cw.*Fy./Wn

    # #calculate nominal bimoment strength including local buckling
    # #for now use DSM local buckling flexural curve
    # Bn, notused = AISIS10016.f321(By, Bcrℓ, ASDorLRFD)

    eBn = calculate_factored_strength(Bn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Bn, eBn

end


#Eqs. G8.1-3 and G8.1-4: nominal bimoment strength considering inelastic reserve and web buckling (Bn, λw) — alternative (a) to the simple yield-strength path of Eq. G8.1-1/-2
function g8_1__3to4(Bp, By, Bcrw, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λw = sqrt(By / Bcrw)

    Bn = Bp * ((1 + 0.09 * λw^2) / (1 + 0.23 * λw^2))

    eBn = calculate_factored_strength(Bn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Bn, eBn

end


#Eqs. G8.1-5 and G8.1-6: nominal bimoment strength considering inelastic reserve and flange buckling (Bn, λf) — alternative (b) to the simple yield-strength path of Eq. G8.1-1/-2
function g8_1__5to6(Bp, By, Bcrf, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λf = sqrt(By / Bcrf)

    Bn = Bp / (1 + 1.11 * λf^2)

    eBn = calculate_factored_strength(Bn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Bn, eBn

end

#---------------------------------------------------------------------------
# CHAPTER H — MEMBERS UNDER COMBINED FORCES
#---------------------------------------------------------------------------

#Eq. H1.1-3 (a and b): available flexural strength considering tension yielding (Maxt or Mayt, depending on which Sft is supplied)
function h1_1__3(Sft, Fy, design_code)

    Ωb = 1.67
    ϕb_LRFD = 0.90
    ϕb_LSD = 0.90

    Mat = calculate_factored_strength(Sft * Fy, Ωb, ϕb_LRFD, ϕb_LSD, design_code)

    return Mat

end


#Eqs. H1.1-1 and H1.1-2: combined tensile axial load and bending interaction ratios
function h1_1__1to2(T, Mx, My, Maxt, Mayt, Max, May, Ta)

    ratio1 = Mx / Maxt + My / Mayt + T / Ta
    ratio2 = Mx / Max + My / May - T / Ta

    return ratio1, ratio2

end


#Eq. H1.2-1: combined compressive axial load and bending interaction ratio
function h1_2__1(P, Mx, My, Pa, Max, May)

    ratio = P / Pa + Mx / Max + My / May

    return ratio

end


#Eqs. H2-1 and H2-2: combined bending and shear interaction ratios (H2-2 applies only when M/Mao > 0.5 and V/Va > 0.7, for members with shear stiffeners)
function h2__1to2(M, V, Mao, Va)

    ratio1 = (M / Mao)^2 + (V / Va)^2
    ratio2 = 0.6 * (M / Mao) + (V / Va)

    return ratio1, ratio2

end


#Eq. H3-1 (a and b): combined bending and web crippling, shapes with a single unreinforced web
function h3__1(P, M, Pn, Mno, design_code)

    Ω = 1.70
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.75

    LHS = 0.91 * (P / Pn) + (M / Mno)
    RHS = design_code == "ASD" ? 1.33 / Ω : (design_code == "LRFD" ? 1.33 * ϕ_LRFD : 1.33 * ϕ_LSD)

    return LHS, RHS, LHS <= RHS

end


#Eq. H3-2 (a and b): combined bending and web crippling, shapes with multiple unreinforced webs providing a high degree of restraint against web rotation (e.g. back-to-back C-sections)
function h3__2(P, M, Pn, Mno, design_code)

    Ω = 1.70
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.75

    LHS = 0.88 * (P / Pn) + (M / Mno)
    RHS = design_code == "ASD" ? 1.46 / Ω : (design_code == "LRFD" ? 1.46 * ϕ_LRFD : 1.46 * ϕ_LSD)

    return LHS, RHS, LHS <= RHS

end


#Eq. H3-3 (a and b): combined bending and web crippling, two nested Z-shapes
function h3__3(P, M, Pn, Mno, design_code)

    Ω = 1.70
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.80

    LHS = 0.86 * (P / Pn) + (M / Mno)
    RHS = design_code == "ASD" ? 1.65 / Ω : (design_code == "LRFD" ? 1.65 * ϕ_LRFD : 1.65 * ϕ_LSD)

    return LHS, RHS, LHS <= RHS

end


#Not an AISI S100-2024 equation — a supplementary free-flange combined bending/torsion interaction check for Z-purlin standing-seam roof systems, kept because PurlinLine.jl depends on it. See h4__1 (Eq. H4-1) for the actual 2024 Specification combined bending and torsion check.
function h42(Mxbar,Mybar,Bbar,Mybar_freeflange, Maxℓo,Mayℓo,Ba, Mayℓo_freeflange)


    ActionMx = abs.(Mxbar./(Maxℓo))
    ActionMy = abs.(Mybar./(Mayℓo))
    ActionB = abs.(Bbar./(Ba))
    ActionMy_freeflange = abs.(Mybar_freeflange./(Mayℓo_freeflange))

    Interaction = ActionMx .+ ActionMy .+ ActionB .+ ActionMy_freeflange

    return ActionMx, ActionMy, ActionB, ActionMy_freeflange, Interaction

end

#Eq. H4-1, Section H4, AISI S100-2024 (combined bending and torsion interaction)
function h4__1(Mxbar, Mybar, Bbar, Max, May, Ba)

    ActionMx = abs.(Mxbar ./ Max)
    ActionMy = abs.(Mybar ./ May)
    ActionB = abs.(Bbar ./ Ba)

    Interaction = ActionMx .+ ActionMy .+ ActionB

    return ActionMx, ActionMy, ActionB, Interaction

end

#---------------------------------------------------------------------------
# APPENDIX 2 — ELASTIC BUCKLING ANALYSIS OF MEMBERS
#---------------------------------------------------------------------------

#Eq. 2.3.1-1: axial force for flexural buckling about the x-axis (Pex)
function appendix2_2_3_1__1(E, Ix, Kx, Lx)
    Pex = π^2 * E * Ix / (Kx * Lx)^2
    return Pex
end

#Eq. 2.3.1-2: axial force for flexural buckling about the y-axis (Pey)
function appendix2_2_3_1__2(E, Iy, Ky, Ly)
    Pey = π^2 * E * Iy / (Ky * Ly)^2
    return Pey
end

#Eq. 2.3.1-3: axial force for torsional buckling about the shear center (Pt)
function appendix2_2_3_1__3(ro, G, J, E, Cw, Kt, Lt)
    Pt = (G * J + π^2 * E * Cw / (Kt * Lt)^2) / ro^2
    return Pt
end

#Eq. 2.3.1-4: coefficient for flexural-torsional buckling about the x-axis (β)
function appendix2_2_3_1__4(xo, ro, Kt, Lt, Kx, Lx)
    beta = 1 - (xo / ro)^2 * (Kt * Lt / (Kx * Lx))^2
    return beta
end

#Eq. 2.3.1-5: coefficient for flexural-torsional buckling about the y-axis (γ)
function appendix2_2_3_1__5(yo, ro, Kt, Lt, Ky, Ly)
    gamma = 1 - (yo / ro)^2 * (Kt * Lt / (Ky * Ly))^2
    return gamma
end

#Eq. 2.3.1-6: asymmetry property (j), by discretized numerical integration over the section coordinates
function appendix2_2_3_1__6(section_coords, Iy, xo, t, xc, yc)
    integral_x3 = 0.0
    integral_xy2 = 0.0
    for index in 1:(size(section_coords, 1) - 1)
        element_length = norm(section_coords[index + 1, :] - section_coords[index, :])
        x = mean(section_coords[index:index + 1, 1]) - xc
        y = mean(section_coords[index:index + 1, 2]) - yc
        dA = element_length * t
        integral_x3 += x^3 * dA
        integral_xy2 += x * y^2 * dA
    end
    j = (integral_x3 + integral_xy2) / (2 * Iy) - xo
    return j
end

#Helper (not a separately numbered equation): product of inertia (Ixy) by discretized numerical integration over the section coordinates, used to form Ixy,avg (Eq. 2.3.1-12) etc.
function appendix2_Ixy(section_coords, t, xc, yc)
    Ixy = 0.0
    for index in 1:(size(section_coords, 1) - 1)
        element_length = norm(section_coords[index + 1, :] - section_coords[index, :])
        x = mean(section_coords[index:index + 1, 1]) - xc
        y = mean(section_coords[index:index + 1, 2]) - yc
        Ixy += x * y * element_length * t
    end
    return Ixy
end

#Eq. 2.3.1-7: polar radius of gyration about the shear center (ro)
function appendix2_2_3_1__7(Ix, A, Iy, xo, yo)
    ro = sqrt(Ix / A + Iy / A + xo^2 + yo^2)
    return ro
end

#Eq. 2.3.1-8: bending coefficient accounting for moment gradient (Cb)
function appendix2_2_3_1__8(Mmax, MA, MB, MC)
    Cb = 12.5 * Mmax / (2.5 * Mmax + 3 * MA + 4 * MB + 3 * MC)
    return Cb
end

#Eq. 2.3.1-9: weighted average cross-sectional area for members with holes (Aavg)
function appendix2_2_3_1__9(Ag, Anet, Lh, s)
    Aavg = (Ag * (s - Lh) + Anet * Lh) / s
    return Aavg
end

#Eq. 2.3.1-10: weighted average moment of inertia about the x-axis for members with holes (Ix,avg)
function appendix2_2_3_1__10(Ixg, Ixnet, Lh, s)
    Ixavg = (Ixg * (s - Lh) + Ixnet * Lh) / s
    return Ixavg
end

#Eq. 2.3.1-11: weighted average moment of inertia about the y-axis for members with holes (Iy,avg)
function appendix2_2_3_1__11(Iyg, Iynet, Lh, s)
    Iyavg = (Iyg * (s - Lh) + Iynet * Lh) / s
    return Iyavg
end

#Eq. 2.3.1-12: weighted average product of inertia for members with holes (Ixy,avg)
function appendix2_2_3_1__12(Ixyg, Ixynet, Lh, s)
    Ixyavg = (Ixyg * (s - Lh) + Ixynet * Lh) / s
    return Ixyavg
end

#Eq. 2.3.1-13: weighted average Saint-Venant torsion constant for members with holes (Javg)
function appendix2_2_3_1__13(Jg, Jnet, Lh, s)
    Javg = (Jg * (s - Lh) + Jnet * Lh) / s
    return Javg
end

#Eq. 2.3.1-14: weighted average asymmetry property for members with holes (javg)
function appendix2_2_3_1__14(jg, jnet, Lh, s)
    javg = (jg * (s - Lh) + jnet * Lh) / s
    return javg
end

#Eq. 2.3.1-15: weighted average shear center x-coordinate for members with holes (xo,avg)
function appendix2_2_3_1__15(xog, xonet, Lh, s)
    xoavg = (xog * (s - Lh) + xonet * Lh) / s
    return xoavg
end

#Eq. 2.3.1-16: weighted average shear center y-coordinate for members with holes (yo,avg)
function appendix2_2_3_1__16(yog, yonet, Lh, s)
    yoavg = (yog * (s - Lh) + yonet * Lh) / s
    return yoavg
end

#Eq. 2.3.1.1.2-1
function appendix2_2_3_1_1_2__1(Pex, Pt, beta)
    Pcre = ((Pex + Pt) - sqrt((Pex + Pt)^2 - 4 * beta * Pex * Pt)) / (2 * beta)
    return Pcre
end

#Eq. 2.3.1.1.4-2
function appendix2_2_3_1_1_4__2(Pey, Pt, gamma)
    Pcre = ((Pey + Pt) - sqrt((Pey + Pt)^2 - 4 * gamma * Pey * Pt)) / (2 * gamma)
    return Pcre
end

#Eq. 2.3.1.1.4-3
function appendix2_2_3_1_1_4__3(Pcre, Pfx, Pfy, Pt, Pfxy, xo, yo, ro, Kt, Lt, Kf, Lf)
    equation = (Pcre - Pfx) * (Pcre - Pfy) * (Pcre - Pt) -
        Pcre^2 * (Pcre - Pfy) * (xo / ro)^2 * (Kt * Lt / (Kf * Lf))^2 -
        Pcre^2 * (Pcre - Pfx) * (yo / ro)^2 * (Kt * Lt / (Kf * Lf))^2 +
        2 * Pcre^2 * Pfxy * (xo * yo / ro^2) * (Kt * Lt / (Kf * Lf))^2 -
        (Pcre - Pt) * Pfxy^2
    return equation
end

#Eq. 2.3.1.1.4-4
function appendix2_2_3_1_1_4__4(E, Ix, Kf, Lf)
    Pfx = π^2 * E * Ix / (Kf * Lf)^2
    return Pfx
end

#Eq. 2.3.1.1.4-5
function appendix2_2_3_1_1_4__5(E, Iy, Kf, Lf)
    Pfy = π^2 * E * Iy / (Kf * Lf)^2
    return Pfy
end

#Eq. 2.3.1.1.4-6
function appendix2_2_3_1_1_4__6(E, Ixy, Kf, Lf)
    Pfxy = π^2 * E * Ixy / (Kf * Lf)^2
    return Pfxy
end

#Eq. 2.3.1.2.1-1
function appendix2_2_3_1_2_1__1(Cb, ro, Pey, Pt)
    Mcre = Cb * ro * sqrt(Pey * Pt)
    return Mcre
end

#Eq. 2.3.1.2.1-2
function appendix2_2_3_1_2_1__2(Cb, d, E, Iy, Ky, Ly)
    Mcre = Cb * π^2 * E * d * Iy / (2 * (Ky * Ly)^2)
    return Mcre
end

#Eq. 2.3.1.2.2-1
function appendix2_2_3_1_2_2__1(Cb, Pex, Cs, j, ro, Pt)
    Mcre = Cb * Pex * (Cs * j + sqrt(j^2 + ro^2 * Pt / Pex))
    return Mcre
end

#Eq. 2.3.1.2.3-2
function appendix2_2_3_1_2_3__2(Cb, d, E, Iy, Ky, Ly)
    Mcre = Cb * π^2 * E * d * Iy / (4 * (Ky * Ly)^2)
    return Mcre
end

#Eq. 2.3.1.2.4-1
function appendix2_2_3_1_2_4__1(Cb, E, Iy, G, J, Ky, Ly)
    Mcre = Cb * π / (Ky * Ly) * sqrt(E * Iy * G * J)
    return Mcre
end

#Eq. 2.3.1.2.5-1
function appendix2_2_3_1_2_5__1(Cb, Pe, j, ro, Pt)
    Mcre = Cb * Pe * (-j + sqrt(j^2 + ro^2 * Pt / Pe))
    return Mcre
end

#Eq. 2.3.1.2.5-2
function appendix2_2_3_1_2_5__2(E, I1, I2, Ib, K, L)
    Pe = π^2 * E * I1 * I2 / (Ib * (K * L)^2)
    return Pe
end

#Eq. 2.3.1.2.5-4
function appendix2_2_3_1_2_5__4(M1r, M2r)
    Mr = sqrt(M1r^2 + M2r^2)
    return Mr
end

#Eq. 2.3.1.2.5-3
function appendix2_2_3_1_2_5__3(I1, I2, M1r, M2r)
    Mr = appendix2_2_3_1_2_5__4(M1r, M2r)
    Ib = I1 * (M1r / Mr)^2 + I2 * (M2r / Mr)^2
    return Ib
end

#Eq. 2.3.1.2.5-5
function appendix2_2_3_1_2_5__5(xa, ya, xo, yo, M1r, M2r, Mr)
    j = (ya - yo) * M1r / Mr - (xa - xo) * M2r / Mr
    return j
end

#Eq. 2.3.2.1-1
function appendix2_2_3_2_1__1(Ag, Fcrl)
    Pcrl = Ag * Fcrl
    return Pcrl
end

#Eq. 2.3.2.1-2 (this formula is also Eq. 2.3.2.2-2, reused as-is since the two are identical)
function appendix2_2_3_2_1__2(k, E, t, mu, w)
    Fcrl = k * π^2 * E / (12 * (1 - mu^2)) * (t / w)^2
    return Fcrl
end

#Eq. 2.3.2.2-1
function appendix2_2_3_2_2__1(Sfc, Fcrl)
    Mcrl = Sfc * Fcrl
    return Mcrl
end

#Eq. 2.3.3.1-1: elastic distortional buckling force for compression members (Pcrd)
function appendix2_2_3_3_1__1(Ag, Fcrd)
    Pcrd = Ag * Fcrd
    return Pcrd
end

#Eq. 2.3.3.1-2: elastic distortional buckling stress for compression members (Fcrd)
function appendix2_2_3_3_1__2(kphife, kphiwe, kphi, kphifg, kphiwg)
    Fcrd = (kphife + kphiwe + kphi) / (kphifg + kphiwg)
    return Fcrd
end

#Eq. 2.3.3.1-3: elastic rotational stiffness provided by the flange to the flange/web juncture (kφfe)
function appendix2_2_3_3_1__3(E, Ixf, xof, xhf, Cwf, Ixyf, Iyf, Ld, G, Jf)
    kphife = (π / Ld)^4 * (E * Cwf + E * Ixf * (xof - xhf)^2 *
        (1 - Ixyf^2 / (Ixf * Iyf))) + (π / Ld)^2 * G * Jf
    return kphife
end

#Eq. 2.3.3.1-4: elastic rotational stiffness provided by the web to the flange/web juncture (kφwe)
function appendix2_2_3_3_1__4(E, t, mu, ho)
    kphiwe = E * t^3 / (12 * (1 - mu^2)) * (2 / ho)
    return kphiwe
end

#Eq. 2.3.3.1-5: geometric rotational stiffness demanded by the flange from the flange/web juncture (k̃φfg)
function appendix2_2_3_3_1__5(Ld, Ixf, Iyf, Af, xhf, yof, xof, Ixyf)
    kphifg = (π / Ld)^2 * (Ixf + Iyf + Af * (xhf^2 + yof^2 -
        2 * yof * (xof - xhf) * Ixyf / Iyf + (xof - xhf)^2 *
        (Ixyf / Iyf)^2))
    return kphifg
end

#Eq. 2.3.3.1-6: geometric rotational stiffness demanded by the web from the flange/web juncture (k̃φwg)
function appendix2_2_3_3_1__6(Ld, t, ho)
    kphiwg = (π / Ld)^2 * t * ho^3 / 60
    return kphiwg
end

#Eq. 2.3.3.1-7: critical unbraced length for distortional buckling of compression members (Lcrd)
function appendix2_2_3_3_1__7(ho, mu, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
    Lcrd = π * ho * ((6 * (1 - mu^2) / (t^3 * ho^3)) *
        (Cwf + Ixf * (xof - xhf)^2 * (1 - Ixyf^2 / (Ixf * Iyf))))^(1 / 4)
    return Lcrd
end

#Eq. 2.3.3.2-1: elastic distortional buckling moment for flexural members (Mcrd)
function appendix2_2_3_3_2__1(Sfc, Fcrd)
    Mcrd = Sfc * Fcrd
    return Mcrd
end

#Eq. 2.3.3.2-2: elastic distortional buckling stress for flexural members (Fcrd), combining the rotational stiffness terms with the moment-gradient factor β
function appendix2_2_3_3_2__2(kphife, kphiwe, kphi, kphifg, kphiwg, beta)
    Fcrd = beta * (kphife + kphiwe + kphi) / (kphifg + kphiwg)
    return Fcrd
end

#Eq. 2.3.3.2-3: moment-gradient factor (β), 1.0 ≤ β ≤ 1.3
function appendix2_2_3_3_2__3(Ld, Lm, M1, M2)
    beta = clamp(1 + 0.4 * (Ld / Lm)^0.7 * (1 + M1 / M2)^0.7, 1.0, 1.3)
    return beta
end

#Eq. 2.3.3.2-4: critical unbraced length for distortional buckling of flexural members bending about the axis perpendicular to the web (Lcrd)
function appendix2_2_3_3_2__4(ho, mu, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
    Lcrd = π * ho * ((4 * (1 - mu^2) / (t^3 * ho^3)) *
        (Cwf + Ixf * (xof - xhf)^2 * (1 - Ixyf^2 / (Ixf * Iyf))) + 1 / 720)^(1 / 4)
    return Lcrd
end

#Eq. 2.3.3.2-5: elastic rotational stiffness provided by the web to the flange/web juncture (k̄φwe), bending about the axis perpendicular to the web
function appendix2_2_3_3_2__5(E, t, mu, ho, Ld)
    kphiwe = (E * t^3 / (12 * (1 - mu^2))) * (3 / ho) *
        (1 + (2 / 15) * (π * ho / Ld)^2 + (1 / 720) * (π * ho / Ld)^4)
    return kphiwe
end

#Eq. 2.3.3.2-6: geometric rotational stiffness demanded by the flange from the flange/web juncture (k̃φfg), bending about the axis perpendicular to the web
function appendix2_2_3_3_2__6(Ld, Ixf, Iyf, Af, xhf, yof, xof, Ixyf)
    kphifg = (π / Ld)^2 * (Ixf + Iyf + Af * (xhf^2 + yof^2 -
        2 * yof * (xof - xhf) * (Ixyf / Iyf)))
    return kphifg
end

#Eq. 2.3.3.2-7: geometric rotational stiffness demanded by the web from the flange/web juncture (k̃φwg), bending about the axis perpendicular to the web
function appendix2_2_3_3_2__7(Ld, t, ho, ξweb)
    πho_Ld = π * ho / Ld
    kphiwg = (π / Ld)^2 * (t * ho^3 / 240) *
        (1110 + 810 * (1 - ξweb) + 8 * πho_Ld^2 + πho_Ld^4) /
        (420 + 28 * πho_Ld^2 + πho_Ld^4)
    return kphiwg
end

#Eq. 2.3.3.2-8: critical unbraced length for distortional buckling of flexural members bending about the axis parallel to the web, web in tension (Lcrd)
function appendix2_2_3_3_2__8(he, mu, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
    Lcrd = π * he * ((6 * (1 - mu^2) / (t^3 * he^3)) *
        (Cwf + Ixf * (xof - xhf)^2 * (1 - Ixyf^2 / (Ixf * Iyf))) + 1 / 120)^(1 / 4)
    return Lcrd
end

#Eq. 2.3.3.2-9: elastic rotational stiffness provided by the web to the flange/web juncture (k̄φwe), bending about the axis parallel to the web
function appendix2_2_3_3_2__9(E, t, mu, he, Ld)
    kphiwe = (E * t^3 / (12 * (1 - mu^2))) * (2 / he) *
        (1 + (1 / 6) * (π * he / Ld)^2 + (1 / 120) * (π * he / Ld)^4)
    return kphiwe
end

#Eq. 2.3.3.2-10: geometric rotational stiffness demanded by the flange from the flange/web juncture (k̃φfg), bending about the axis parallel to the web
function appendix2_2_3_3_2__10(Ld, Ixf, Iyf, Af, xhf, yof, xof, Ixyf, ψf, ξf)
    kphifg = (π / Ld)^2 * (Ixf + Iyf + Af * (xhf^2 + yof^2 -
        2 * yof * (xof - xhf) * (Ixyf / Iyf))) * ψf + (π / Ld)^2 * Iyf * ξf
    return kphifg
end

#Eq. 2.3.3.2-11: effective web depth for bending about the axis parallel to the web (he) — k̃φwg = 0 for this case, so there is no corresponding function
function appendix2_2_3_3_2__11(Iyf, Af, xhf)
    he = 3.5 * sqrt(Iyf / Af + xhf^2)
    return he
end

#Dispatcher spanning Sections 2.3.1.1.1-2.3.1.1.4 (Fcre, Pcre for compression members) — selects the applicable Pcre formula by section_type
function appendix2_Fcre_compression(section_type, Ag, E, G, Ix, Iy, J, Cw,
    xo, yo, Kx, Lx, Ky, Ly, Kt, Lt)

    if section_type in ("non_symmetric", "nonsymmetric")
        throw(ArgumentError("Non-symmetric compression Fcre calculation is not implemented yet"))
    end

    ro = appendix2_2_3_1__7(Ix, Ag, Iy, xo, yo)
    Pex = appendix2_2_3_1__1(E, Ix, Kx, Lx)
    Pey = appendix2_2_3_1__2(E, Iy, Ky, Ly)
    Pt = appendix2_2_3_1__3(ro, G, J, E, Cw, Kt, Lt)
    beta = appendix2_2_3_1__4(xo, ro, Kt, Lt, Kx, Lx)
    gamma = appendix2_2_3_1__5(yo, ro, Kt, Lt, Ky, Ly)

    if section_type in ("doubly_symmetric", "no_torsional", "flexural")
        Pcre = min(Pex, Pey)
    elseif section_type in ("singly_symmetric", "singly_symmetric_major",
            "singly_symmetric_about_major_axis", "C", "Cee")
        Pcre = min(Pex, Pey, appendix2_2_3_1_1_2__1(Pex, Pt, beta))
    elseif section_type in ("singly_symmetric_minor", "singly_symmetric_about_minor_axis")
        Pcre = min(Pex, Pey, appendix2_2_3_1_1_4__2(Pey, Pt, gamma))
    elseif section_type in ("point_symmetric", "doubly_symmetric_torsional", "Z", "Zee")
        Pcre = min(Pex, Pey, Pt)
    else
        throw(ArgumentError("Unsupported section_type: $(section_type)"))
    end

    Fcre = Pcre / Ag
    return Fcre
end

#Dispatcher spanning Sections 2.3.1.2.1-2.3.1.2.4 (Fcre, Mcre for flexural members) — selects the applicable Mcre formula by section_type
function appendix2_Fcre_flexural(section_type, Sfc, Cb, A, E, G, Ix, Iy, J, Cw,
    xo, yo, Ky, Ly, Kt, Lt)

    if section_type in ("non_symmetric", "nonsymmetric", "singly_symmetric_minor",
            "singly_symmetric_about_minor_axis", "biaxial")
        throw(ArgumentError("Flexural Fcre calculation is not implemented yet for section type: $(section_type)"))
    end

    ro = appendix2_2_3_1__7(Ix, A, Iy, xo, yo)
    Pey = appendix2_2_3_1__2(E, Iy, Ky, Ly)
    Pt = appendix2_2_3_1__3(ro, G, J, E, Cw, Kt, Lt)

    if section_type in ("singly_symmetric", "singly_symmetric_major",
            "singly_symmetric_about_major_axis", "doubly_symmetric", "C", "Cee")
        Mcre = appendix2_2_3_1_2_1__1(Cb, ro, Pey, Pt)
    elseif section_type in ("point_symmetric", "Z", "Zee")
        Mcre = Cb * ro / 2 * sqrt(Pey * Pt)
    elseif section_type == "closed_box"
        Mcre = appendix2_2_3_1_2_4__1(Cb, E, Iy, G, J, Ky, Ly)
    else
        throw(ArgumentError("Unsupported flexural section_type: $(section_type)"))
    end

    Fcre = Mcre / Sfc
    return Fcre
end

end #module
