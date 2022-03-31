module v24

export h411, h42

function calculate_factored_strength(Rn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    if design_code == "AISI S100-16 ASD"
        eRn  = Rn / Ω
    elseif design_code == "AISI S100-16 LRFD"
        eRn = Rn * ϕ_LRFD
    elseif design_code == "AISI S100-16 LSD"
        eRn = Rn * ϕ_LSD
    elseif design_code == "nominal"
        eRn = Rn
    end

    return eRn

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
