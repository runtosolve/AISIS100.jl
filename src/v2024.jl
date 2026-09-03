module v2024

export h411, h42

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


function e2(Fcre, Fy, Ag, design_code)

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



function e32(Pne, Pcrℓ, design_code)

    Ω = 1.80
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80

    λℓ = sqrt(Pne/Pcrℓ)

    Pnℓ = min(1.2*Pne*((1 + 0.10*λℓ^2)/(1 + 0.55*λℓ^2)), Pne)

    ePnℓ = calculate_factored_strength(Pnℓ, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    results = (Rn=Pnℓ, eRn=ePnℓ)

    return results

end


function e4(Py, Pcrd, design_code)

    Ω = 1.80
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80

    λd = sqrt(Py/Pcrd)

    Pnd = min(1.2*Py*((1 + 0.05*λd^2)/(1 + 0.67*λd^2)), Py)

    ePnd = calculate_factored_strength(Pnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    results = (Rn=Pnd, eRn=ePnd)

    return results

end


function f32(Mne, Mcrℓ, ks, αs, βs, My3, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λℓ=sqrt(Mne/Mcrℓ)


    Mnℓ = min(ks*Mne*((1+0.10*αs*λℓ*λℓ)/(1+0.55*βs*λℓ*λℓ)),My3)


    # eMnℓ = Mnℓ * StrengthFactor

    eMnℓ = calculate_factored_strength(Mnℓ, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnℓ, eMnℓ

end



function f32_holes(Mne, Mcrℓ, ks, αs, βs, Mpnet, Mynet, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λℓ = sqrt(Mne/Mcrℓ)

    My3net = Mpnet - (Mpnet - Mynet)/9   # Eq. F3.2-8

    Mnℓ = min(Mpnet*((1+0.10*αs*λℓ*λℓ)/(1+0.55*βs*λℓ*λℓ*(Mpnet/(ks*Mne)))), My3net)

    eMnℓ = calculate_factored_strength(Mnℓ, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnℓ, eMnℓ

end


function f4(My, Mcrd, Mp, αs, βs, My3, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λd = sqrt(My/Mcrd)

    Mnd = min(Mp*((1 + 0.07*αs*λd^2)/(1 + 0.60*βs*λd^2)), My3)

    eMnd = calculate_factored_strength(Mnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnd, eMnd

end



function f4_holes(My, Mcrd, Mp, αs, βs, Mpnet, Mynet, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λd = sqrt(My/Mcrd)

    My3net = Mpnet - (Mpnet - Mynet)/9   # Eq. F3.2-8

    Mnd = min(Mpnet*((1 + 0.07*αs*λd^2)/(1 + 0.60*βs*λd^2*(Mpnet/Mp))), My3net)

    eMnd = calculate_factored_strength(Mnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnd, eMnd

end


# Eq. G2.1-3, yield shear force of the cross-section.
function g213(h, t, Fy)
    Vy = 0.6 * h * t * Fy
    return Vy
end

function g2_1__1_2_3(Vcr, Vy, design_code)  #no transverse stiffeners

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.75  #confirmed correct per S100-24 G2 header (φv = 0.75 LSD)

    λv=sqrt(Vy/Vcr)


    Vn = min(1.2*Vy/(1+0.57*λv*λv),Vy)

    eVn = calculate_factored_strength(Vn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Vn, eVn

end


# Eq. G2.3-4, -5: shear buckling coefficient for webs with transverse
# stiffeners satisfying Section G4, a = shear panel length (clear distance
# between stiffeners). Case (a), kv = 5.34 for unreinforced webs -- pass
# a = nothing.
function g234_5(h; a = nothing)
    if isnothing(a)
        return 5.34
    elseif a/h <= 1.0
        return 4.00 + 5.34/(a/h)^2
    else
        return 5.34 + 4.00/(a/h)^2
    end
end

# Eq. G2.3-1, -2, -3: web elastic critical shear buckling force. kv comes
# from g234_5 for webs without holes, or g38 for webs with holes (Section G3).
function g231_2_3(E, μ, kv, h, t)
    Aw  = h * t
    Fcr = (π^2 * E * kv) / (12 * (1 - μ^2) * (h / t)^2)
    Vcr = Aw * Fcr
    return Vcr
end


function g31_5(d_h, L_h, h, Vy)

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

function g36_7(L_h, d_h, h, Vcr)

    α_vh = (1-0.4*(L_h-d_h)/h)^2

    V_crh = α_vh * Vcr

    return α_vh, V_crh

end

function g38(a, h, d_h, b_f)

    kv = 4.86 + 6.15 * (h / a) - 3.63 * (d_h / h) - 19.58 * (d_h / a) + 13.88 * (d_h^2 / (a * h)) + 0.57* (b_f / h)

end


# Eq. G3-9, -10. A_h (area of the circular or slotted hole) is a caller-
# supplied input rather than computed here, since its formula depends on
# hole shape (circular: A_h = π(d_h/2)^2; slotted/oval "stadium" shape:
# A_h = width*(length-width) + π(width/2)^2, etc.) -- this function only
# implements the shape-agnostic dh_eq/Lh_eq equivalence once A_h is known.
function g39_10(A_h, d_h, L_h)

    d_h_eq = (0.003*(L_h/d_h) + 0.822) * d_h
    L_h_eq = 0.865*A_h/d_h_eq

    return d_h_eq, L_h_eq

end





function g51(t, h, Fy, θ, C, C_R, R, C_N, N, C_h, ϕ_w, Ω_w, ϕ_w_LSD, design_code)

    Pn = C * t^2 * Fy * sin(deg2rad(θ)) * (1-C_R*sqrt(R/t)) * (1+C_N*sqrt(N/t)) * (1-C_h*sqrt(h/t))

    ePn = calculate_factored_strength(Pn, Ω_w, ϕ_w, ϕ_w_LSD, design_code)

    return Pn, ePn

end


function h21(Mbar, Vbar, Maℓo, Va)

    interaction = sqrt.((Mbar./Maℓo).^2 .+ (Vbar./Va).^2)

    return interaction

end


function h31(Pbar, Mbar, Pn, Mnℓo, design_code)

    actionP = 0.91 * (Pbar/Pn)
    actionM = Mbar / Mnℓo
    interaction = actionP + actionM

    Ω = 1.70
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.75

    if design_code == "ASD"
        DC = interaction / (1.33/Ω)
    elseif design_code == "LRFD"
        DC = interaction / (1.33*ϕ_LRFD)
    elseif design_code == "LSD"
        DC = interaction / (1.33*ϕ_LSD)
    elseif design_code == "nominal"
        DC = interaction / 1.33
    end

    return actionP, actionM, interaction, DC

end

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


function h42(Mxbar,Mybar,Bbar,Mybar_freeflange, Maxℓo,Mayℓo,Ba, Mayℓo_freeflange)


    ActionMx = abs.(Mxbar./(Maxℓo))
    ActionMy = abs.(Mybar./(Mayℓo))
    ActionB = abs.(Bbar./(Ba))
    ActionMy_freeflange = abs.(Mybar_freeflange./(Mayℓo_freeflange))

    Interaction = ActionMx .+ ActionMy .+ ActionB .+ ActionMy_freeflange

    return ActionMx, ActionMy, ActionB, ActionMy_freeflange, Interaction

end



function l21(Md, M, Ig)

    Ieff = Ig * (Md/M)

    if Ieff > Ig

        Ieff = Ig
    end

    return Ieff

end



end #module
