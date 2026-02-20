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


function g39_10(d_h, L_h)  

    A_h = π*(d_h/2)^2  #for circular hole only, not slotted hole
    d_h_eq = (0.003*(L_h/d_h) + 0.822) * d_h
    L_h_eq = 0.865*A_h/d_h_eq

    return A_h, d_h_eq, L_h_eq

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



end #module
