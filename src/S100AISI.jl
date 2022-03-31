module S100AISI


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

export v16
include("v16.jl")
using .v16

export v24
include("v24.jl")
using .v24

end # module
