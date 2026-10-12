module v16S3

using CSV, DataFrames, Unitful, Statistics, LinearAlgebra


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


function d21(;Ag, Fy, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    Tn = Ag * Fy

    eTn = calculate_factored_strength(Tn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Tn, eTn

end


function d31(;An, Fu, design_code)

    Ω = 2.00
    ϕ_LRFD = 0.75
    ϕ_LSD = 0.75

    Tn = An * Fu

    eTn = calculate_factored_strength(Tn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Tn, eTn

end

function e2(;Fcre, Fy, Ag, design_code)

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

    return Pne, ePne

end

function e321(;Pne, Pcrℓ, design_code)


    λℓ = sqrt(Pne/Pcrℓ)

    if λℓ <= 0.776
        Pnℓ = Pne
    elseif λℓ > 0.776
        Pnℓ = (1 - 0.15 * (Pcrℓ/Pne)^0.4) * (Pcrℓ/Pne)^0.4 * Pne
    end

    Ω = 1.80
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80

    if design_code == "AISI S100-16 ASD"
        ePnℓ  = Pnℓ / Ω
     elseif design_code == "AISI S100-16 LRFD"
        ePnℓ = Pnℓ * ϕ_LRFD
     elseif design_code == "AISI S100-16 LSD"
        ePnℓ = Pnℓ * ϕ_LSD
     elseif design_code == "nominal"
         ePnℓ = Pnℓ
    end

    return Pnℓ, ePnℓ

end


function e41(;Py, Pcrd, design_code)

    Ω = 1.80
    ϕ_LRFD = 0.85
    ϕ_LSD = 0.80
  
    λc = sqrt(Py/Pcrd)

    if λc <= 0.561
        Pnd = Py
    else
        Pnd = (1 - 0.25 * (Pcrd/Py)^0.6) * (Pcrd/Py)^0.6 * Py
    end

    ePnd = calculate_factored_strength(Pnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Pnd, ePnd

end


function f2_1__3_4_5(Fcre, Fy)

    if Fcre>=2.78*Fy
        Fn = Fy
    elseif (Fcre<(2.78*Fy)) & (Fcre>(0.56*Fy))
        Fn = 10/9 * Fy * (1 - (10*Fy)/(36*Fcre))
    elseif Fcre <= 0.56*Fy
        Fn = Fcre
    end

    return Fn

end


function f2_1__1_2(;Fcre, Fy, Sf, Sfy, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90
  
    # if Fcre>=2.78*Fy
    #     Fn = Fy
    # elseif (Fcre<2.78*Fy) & (Fcre>0.56*Fy)
    #     Fn = 10/9 * Fy * (1 - (10*Fy)/(36*Fcre))
    # elseif Fcre <= 0.56*Fy
    #     Fn = Fcre
    # end

    Fn = f2_1__3_4_5(Fcre, Fy)

    Mne = minimum([Sf*Fn, Sfy*Fy])

    eMne = calculate_factored_strength(Mne, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mne, eMne
 
end



function f321(Mne, Mcrℓ, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λℓ=sqrt(Mne/Mcrℓ)

    if λℓ <= 0.776
        Mnℓ=Mne
    else
        Mnℓ=(1-0.15*(Mcrℓ/Mne)^0.4)*(Mcrℓ/Mne)^0.4*Mne
    end

    # eMnℓ = Mnℓ * StrengthFactor

    eMnℓ = calculate_factored_strength(Mnℓ, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnℓ, eMnℓ

end

function f322(Mne, Mcrℓ, My_net, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    #Mcrℓ includes the influence of holes here.
    Mnℓ, eMnℓ = f321(Mne, Mcrℓ, design_code)

    if Mnℓ > My_net

        Mnℓ = My_net 
        eMnℓ = calculate_factored_strength(Mnℓ, Ω, ϕ_LRFD, ϕ_LSD, design_code)
        
    end

    # eMnℓ = calculate_factored_strength(Mnℓ, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnℓ, eMnℓ

end


#inelastic reserve flexural local buckling
function f323(My, Mcrl, Sc, St, Z, Fy, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    lambda_l = sqrt(My/Mcrl)

    Cyl = sqrt(0.776/lambda_l)

    if Cyl > 3

        Cyl = 3

    end

    Mp = Fy * Z

    Myc = Fy * Sc

    Cyt = 3
    Myt3 = My + (1 - 1 / Cyt^2) * (Mp - My)

    if Sc <= St  #first yield in compression
        Mnl = My + (1 - 1/Cyl^2) * (Mp-My)
    elseif Sc > St  #first yield in tension
        Mnl = Myc + (1-1/Cyl^2) * (Mp-My)

        if Mnl > Myt3
            Mnl = Myt3
        end

    end

    eMnl = calculate_factored_strength(Mnl, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return lambda_l, Cyl, Mp, Myc, Myt3, Mnl, eMnl

end

f323(;My, Mcrl, Sc, St, Z, Fy, design_code) = f323(My, Mcrl, Sc, St, Z, Fy, design_code)


function f411(My, Mcrd, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λd=sqrt(My/Mcrd)

    if λd <= 0.673
        Mnd=My
    else
        Mnd=(1-0.22*(Mcrd/My)^0.5)*(Mcrd/My)^0.5*My   #Eq. F4.1-2 (exponent 0.5; was 0.4 before September 2026)
    end

    eMnd = calculate_factored_strength(Mnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mnd, eMnd

end


function f421_6(My, Mcrd, My_net, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    λd = sqrt(My/Mcrd)
    λd1 = 0.673*(My_net/My)^3
    λd2 = 0.673*(1.7*(My/My_net)^(2.7)-0.7)

    if λd <= λd1
        Mnd = My_net
        eMnd = calculate_factored_strength(Mnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)
    
    elseif (λd > λd1) & (λd <= λd2)

        Md2 = (1-0.22*(1/λd2)) * (1/λd2) * My
        Mnd = My_net - ((My_net - Md2) / (λd2 - λd1)) * (λd - λd1)

        Mnd_nh, eMnd_nh = f411(My, Mcrd, design_code)

        if Mnd > Mnd_nh

            Mnd = Mnd_nh

        end

        eMnd = calculate_factored_strength(Mnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    else

        Mnd, eMnd = f411(My, Mcrd, design_code)

    end

    return Mnd, eMnd

end



#inelastic reserve flexural distortional buckling
function f43(My, Mcrd, Sc, St, Z, Fy, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    lambda_d = sqrt(My/Mcrd)
    Cyd = sqrt(0.673/lambda_d)

    if Cyd > 3
        Cyd = 3
    end

    Mp = Fy * Z

    Myc = Fy * Sc

    Cyt = 3
    Myt3 = My + (1-1/Cyt^2)*(Mp - My)

    if Sc <= St  #first yield in compression
        Mnd = My + (1 - 1/Cyd^2)*(Mp-My)
    elseif Sc > St  #first yield in tension
        Mnd = Myc + (1-1/Cyd^2)*(Mp-My)
        if Mnd > Myt3
            Mnd = Myt3
        end
    end

    eMnd = calculate_factored_strength(Mnd, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return lambda_d, Cyd, Mp, Myc, Myt3, Mnd, eMnd

end


# function g21_3(Vcr, Vy, design_code)

#     Ω = 1.60
#     ϕ_LRFD = 0.95
#     ϕ_LSD = 0.80

#     # Aw, Vy = g215_6(h, t, Fy)
#     λv=sqrt.(Vy/Vcr)

#     if λv <= 0.815
#         Vn = Vy
#     elseif (λv>0.815) & (λv<=1.227)
#         Vn = 0.815 *sqrt(Vcr*Vy)
#     elseif λv > 1.227
#         Vn = Vcr
#     end

#     eVn = calculate_factored_strength(Vn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

#     return Vn, eVn

# end

function g2_1__1_2_3(Vcr, Vy, design_code)  #no transverse stiffeners

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.75  #check this?  seems low compared to f3 S100-16

    λv=sqrt(Vy/Vcr)

    if λv <= 0.587
        Vn=Vy
    else
        Vn=(1-0.25*(Vcr/Vy)^0.65)*(Vcr/Vy)^0.65*Vy
    end

    eVn = calculate_factored_strength(Vn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Vn, eVn
    
end




function g215_6(h, t, Fy)

    Aw = h * t
    Vy=0.6 * Aw * Fy

    return Aw, Vy

end


function g231(h, t, Fcr)

    Aw = h.*t
    Vcr = Aw.*Fcr

    return Vcr

end

function g232(E, μ, kv, h, t)

    Fcr = (π^2 .*E.*kv)./(12 .*(1-μ.^2).*(h./t).^2)

end

function g233(a, h)

    if a./h <= 1.0
        kv = 4.00 .+ 5.34./(a./h).^2
    elseif a./h > 1.0
        kv = 5.34 .+ 4.00./(a./h).^2
    end

    return kv

end



function g51(t, h, Fy, θ, C, C_R, R, C_N, N, C_h, ϕ_w, Ω_w, ϕ_w_LSD, design_code)

    Pn = C * t^2 * Fy * sin(deg2rad(θ)) * (1-C_R*sqrt(R/t)) * (1+C_N*sqrt(N/t)) * (1-C_h*sqrt(h/t))

    ePn = calculate_factored_strength(Pn, Ω_w, ϕ_w, ϕ_w_LSD, design_code)

    return Pn, ePn

end


function table_g52()

    filename = string(@__DIR__, "/assets/AISI_S100_16_Table_G5_2.csv")
    data = CSV.File(filename)

    table = DataFrame(support_condition = Vector(data.support_condition), flange_condition = Vector(data.flange_condition), load_case = Vector(data.load_case), load_location=Vector(data.load_location), C=Vector(data.C), C_R = Vector(data.C_R), C_N = Vector(data.C_N), C_h = Vector(data.C_h), ASD = Vector(data.ASD), LRFD = Vector(data.LRFD), LSD = Vector(data.LSD), limits = Vector(data.limits))

    return table

end


function table_g53()

    filename = string(@__DIR__, "/assets/AISI_S100_16_Table_G5_3.csv")
    data = CSV.File(filename)

    table = DataFrame(support_condition = Vector(data.support_condition), flange_condition = Vector(data.flange_condition), load_case = Vector(data.load_case), load_location=Vector(data.load_location), C=Vector(data.C), C_R = Vector(data.C_R), C_N = Vector(data.C_N), C_h = Vector(data.C_h), ASD = Vector(data.ASD), LRFD = Vector(data.LRFD), LSD = Vector(data.LSD), limits = Vector(data.limits))

    return table

end


function h21(Mbar, Vbar, Maℓo, Va)

    Interaction = sqrt.((Mbar./Maℓo).^2 .+ (Vbar./Va).^2)

end

function h121(Pbar, Mxbar, Mybar, Pa, Max, May)

    action_P = abs.(Pbar./Pa)
    action_Mx = abs.(Mxbar./Max)
    action_My = abs.(Mybar./May)

    interaction = action_P .+ action_Mx .+ action_My

    return action_P, action_Mx, action_My, interaction

end

h121(;Pbar, Mxbar, Mybar, Pa, Max, May) = h121(Pbar, Mxbar, Mybar, Pa, Max, May)



function h31(Pbar, Mbar, Pn, Mnℓo, design_code)

    actionP = 0.91 * (Pbar/Pn)
    actionM = Mbar / Mnℓo
    interaction = actionP + actionM

    Ω = 1.70
    ϕ_LRFD = 0.95
    ϕ_LSD = 0.75

    if design_code == "AISI S100-16 ASD"
        DC = interaction / (1.33/Ω)
    elseif design_code == "AISI S100-16 LRFD"
        DC = interaction / (1.33*ϕ_LRFD)
    elseif design_code == "AISI S100-16 LSD"
        DC = interaction / (1.33*ϕ_LSD)
    elseif design_code == "nominal"
        DC = interaction / 1.33
    end

    return actionP, actionM, interaction, DC

end

function h33(Pbar, Mbar, Pn, Mnℓo, design_code)

    actionP = 0.86 * (Pbar/Pn)
    actionM = Mbar / Mnℓo
    interaction = actionP + actionM

    Ω = 1.70
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.80

    if design_code == "AISI S100-16 ASD"
        DC = interaction / (1.65/Ω)
    elseif design_code == "AISI S100-16 LRFD"
        DC = interaction / (1.65*ϕ_LRFD)
    elseif design_code == "AISI S100-16 LSD"
        DC = interaction / (1.65*ϕ_LSD)
    end

    return actionP, actionM, interaction, DC

end

function j26(L, t1, Fu1, t2, Fu2, tw, Fxx, loading_direction, design_code, thickness_units)

    t = minimum([t1, t2])

    if loading_direction == "longitudinal"

        if L/t < 25.0

            Pnv1 = (1 - 0.01*L/t1) * L * t1 * Fu1

            Pnv2 = (1 - 0.01*L/t2) * L * t2 * Fu2

            Pnv = minimum([Pnv1, Pnv2])

            Ω = 2.55
            ϕ_LRFD = 0.60
            ϕ_LSD = 0.50

        elseif L/t >= 25.0

            Pnv1 = 0.75 * t1 * L * Fu1
            Pnv2 = 0.75 * t2 * L * Fu2

            Pnv = minimum([Pnv1, Pnv2])

            Ω = 3.05
            ϕ_LRFD = 0.50
            ϕ_LSD = 0.40

        end

    elseif loading_direction == "transverse"

        Pnv1 = t1 * L * Fu1
        Pnv2 = t2 * L * Fu2

        Pnv = minimum([Pnv1, Pnv2])

        Ω = 2.35
        ϕ_LRFD = 0.65
        ϕ_LSD = 0.60

    end
        
    if ((thickness_units == "inches") & (t > 0.10)) | ((thickness_units == "mm") &  (t > 2.54)) 

        Pn = 0.75 * tw * L * Fxx

        Ω = 2.55
        ϕ_LRFD = 0.60
        ϕ_LSD = 0.50

        if Pn < Pnv

            Pnv = Pn

        end

    else

        Pn = NaN

    end

    # if design_code == "AISI S100-16 ASD"
    #     ePnv  = Pnv / Ω
    # elseif design_code == "AISI S100-16 LRFD"
    #     ePnv = Pnv * ϕ_LRFD
    # elseif design_code == "AISI S100-16 LSD"
    #     ePnv = Pnv * ϕ_LSD
    # elseif design_code == "nominal"
    #     ePnv = Pnv
    # end

    ePnv = calculate_factored_strength(Pnv, Ω, ϕ_LRFD, ϕ_LSD, design_code)


    return Pnv1, Pnv2, Pn, Pnv, ePnv

end


function tablej3311(;d, t, hole_shape)

    if (t >= 0.024u"inch") & (t <= 0.1875u"inch")

        if hole_shape == "standard hole"
    
            if d/t < 10.0

                C = 3.0

            elseif (d/t >= 10.0) & (d/t <= 22.0)

                C = 4.0 - 0.1 * (d/t)

            else d/t > 22

                C = 1.8

            end

        elseif hole_shape == "oversized hole"

            if d/t < 7.0

                C = 3.0

            elseif (d/t >= 7.0) & (d/t <= 18.0)

                C = 1.0 + 14.0 * (d/t)

            else d/t > 18

                C = 1.8

            end

        end

    end

    return C

end

function tablej3312(;connection_type, hole_shape, washers)

    #Only consider standard holes for now.

    if (connection_type == "single shear") & (hole_shape=="standard hole") & (washers == "yes")

        mf = 1.0

    elseif (connection_type == "single shear") & (hole_shape=="standard hole") & (washers == "no")

        mf = 0.75

    end

    return mf

end




function j3311(;C, mf, d, t, Fu, design_code)

    Pnb = C * mf * d * t * Fu 

    Ω = 2.50
    ϕ_LRFD = 0.60
    ϕ_LSD = 0.50

    if design_code == "AISI S100-16 ASD"
        ePnb  = Pnb / Ω
     elseif design_code == "AISI S100-16 LRFD"
        ePnb = Pnb * ϕ_LRFD
     elseif design_code == "AISI S100-16 LSD"
        ePnb = Pnb * ϕ_LSD
     elseif design_code == "nominal"
         ePnb = Pnb
    end

    return Pnb, ePnb 

end


function j611(Anv, Fu, design_code, connection_type)

    Pnv = 0.6 * Fu * Anv

    if connection_type == "bolts"
        Ω = 2.22
        ϕ_LRFD = 0.65
        ϕ_LSD = 0.75

        if design_code == "AISI S100-16 ASD"
            ePnv  = Pnv / Ω
         elseif design_code == "AISI S100-16 LRFD"
            ePnv = Pnv * ϕ_LRFD
         elseif design_code == "AISI S100-16 LSD"
            ePnv = Pnv * ϕ_LSD
         elseif design_code == "nominal"
             ePnv = Pnv
        end

    end

    return Pnv, ePnv

end

function k2112(Cϕ, Mm, Fm, Pm, βo, VM, VF, CP, VP, VQ)

    ϕ = Cϕ*(Mm*Fm*Pm)exp(-βo*sqrt(VM^2+VF^2+CP*VP^2+VQ^2))

    return ϕ

end

function k2114(n)

    m = n-1

    if n>=4
        CP = ((1+1/n)*m)/(m-2)
    elseif n==3
        CP = 5.7
    end
    
    return CP

end


function k2117(n, Rt, Rn)

    Cc = (n*sum(Rt.*Rn) - sum(Rt)*sum(Rn))/((sqrt(n*sum(Rt.^2)-sum(Rt)^2))*sqrt(n*sum(Rn.^2)-sum(Rn)^2))

    return Cc

end



function l21(Md, M, Ig)

    Ieff = Ig * (Md/M)

    if Ieff > Ig

        Ieff = Ig
    end

    return Ieff

end


function appAj341(;Ab, Fn, design_code)

    Ω = 2.00  
    ϕ = 0.75

    Pn = Ab * Fn

    if design_code == "AISI S100-16 ASD"
       ePn  = Pn / Ω
    elseif design_code == "AISI S100-16 LRFD"
        ePn = Pn * ϕ
    elseif design_code == "nominal"
        ePn = Pn
    end

    return Pn, ePn

end



function table_2_3_3__1(CorZ,t,b,d,θ)

    CorZ = convert(Int8, CorZ)

    Af=(b+d)*t
    Jf=1/3*b*t^3+1/3*d*t^3
    Cwf=0.0

    θ = deg2rad(θ)  #convert degrees to radians

    if CorZ==0

        Ixf=t*(t^2*b^2+4*b*d^3+t^2*b*d+d^4)/(12*(b+d))   #Table 2.3.1.3-1 (duplicated t^2*b*d term removed September 2026)
        Iyf=t*(b^4+4*d*b^3)/(12*(b+d))
        Ixyf=(t*b^2*d^2)/(4*(b+d))
        xof=b^2/(2*(b+d))
        xhf=-(b^2+2*d*b)/(2*(b+d))
        yhf=-d^2/(2*(b+d))
        yof=yhf

    else

        Ixf=t*(t^2*b^2+4*b*d^3-4*b*d^3*cos(θ)^2+t^2*b*d+d^4-d^4*cos(θ)^2)/(12*(b+d))
        Iyf=t*(b^4+4*d*b^3+6*d^2*b^2*cos(θ)+4*d^3*b*cos(θ)^2+d^4*cos(θ)^2)/(12*(b+d))
        Ixyf=t*b*d^2*sin(θ)*(b+d*cos(θ))/(4*(b+d))
        xof=(b^2-d^2*cos(θ))/(2*(b+d))
        xhf=-(b^2+2*d*b+d^2*cos(θ))/(2*(b+d))
        yhf=-d^2*sin(θ)/(2*(b+d))
        yof=yhf
    end

    return Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,xhf,yhf,yof

end


function appendix2_2_3_1__1(E, Ix, Kx, Lx)

    Pex = π^2 * E * Ix / (Kx * Lx)^2

end

function appendix2_2_3_1__2(E, Iy, Ky, Ly)

    Pey = π^2 * E * Iy / (Ky * Ly)^2

end

function appendix2_2_3_1__3(ro, G, J, E, Cw, Kt, Lt)

    Pt = (1 / ro^2) * (G * J + π^2 * E * Cw / (Kt * Lt)^2)

end

function appendix2_2_3_1__4(xo, ro, Kt, Lt, Kx, Lx)

   β = 1 - (xo / ro)^2 * ((Kt * Lt) / (Kx * Lx))^2

end

function appendix2_2_3_1__6(section_coords, Iy, xo, t, xc, yc)

    num_elements = size(section_coords)[1] - 1

    integral_x3 = 0.0
    integral_xy2 = 0.0
    for i =1:num_elements

        element_length = norm(section_coords[i+1, :] - section_coords[i, :])
        x_element = mean([section_coords[i, 1], section_coords[i+1, 1]]) - xc
        y_element = mean([section_coords[i, 2], section_coords[i+1, 2]]) - yc
        dA = element_length * t
        integral_x3 += x_element^3 * dA
        integral_xy2 += x_element * y_element^2 * dA 

    end

    j = 1 / (2 * Iy) * (integral_x3 + integral_xy2) - xo

    return j 

end

function appendix2_2_2_3_1__7(Ix, A, Iy, xo, yo)

    ro = sqrt(Ix/A + Iy/A + xo^2 + yo^2)

end

function appendix2_2_3_3_1__1(Ag, Fcrd)

    Pcrd = Ag * Fcrd

end

function appendix2_2_3_3_1__2(kϕfe, kϕwe, kϕ, kϕfg, kϕwg)

    Fcrd = (kϕfe +kϕwe +kϕ) / (kϕfg + kϕwg)

end

function appendix2_2_3_3_1__3(E, Ixf, xof, xhf, Cwf, Ixyf, Iyf, Ld, G, Jf)

    kϕfe = (π / Ld)^4 * (E * Cwf + E * Ixf * (xof - xhf)^2 * (1 - Ixyf^2 / (Ixf * Iyf))) + (π / Ld)^2 * G * Jf

end

function appendix2_2_3_3_1__4(E, t, μ, ho)

    kϕwe = E * t^3 / (12 * (1 - μ^2)) * (2 / ho)

end

function appendix2_2_3_3_1__5(Ld, Ixf, Iyf, Af, xhf, yof, xof, Ixyf)

    kϕfg = (π / Ld)^2 * (Ixf + Iyf + Af * (xhf^2 + yof^2 - 2 * yof * (xof - xhf) * (Ixyf / Iyf)+ (xof - xhf)^2 * (Ixyf / Iyf)^2))

end

function appendix2_2_3_3_1__6(Ld, t, ho)

    kϕwg = (π / Ld)^2 * t * ho^3 / 60

end


function appendix2_2_3_3_1__7(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)

    Lcrd = π * ho * (((6 * (1 - μ^2)) / (t^3 * ho^3)) * (Cwf + Ixf * (xof - xhf)^2 * (1 - Ixyf^2 / (Ixf * Iyf)))) ^ (1/4)

    # Lcrd = ((6 * π^4 * ho * (1 - μ^2))/t^3) * (Ixf * (xof - xhf)^2 + Cwf - Ixyf^2/Iyf*(xof - xhf)^2)^ (1/4)

end

function appendix2_2_3_3_2__4(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)

    Lcrd = π * ho * (((4 * (1 - μ^2)) / (t^3 * ho^3)) * (Cwf + Ixf * (xof - xhf)^2 * (1 - Ixyf^2 / (Ixf * Iyf))) + 1/720) ^ (1/4)

end

function appendix2_2_3_1_2_1__1(Cb, ro, Pey, Pt)

    Mcre = Cb * ro * sqrt(Pey * Pt)

end

function appendix2_2_3_1_2_2__1(Cb, Pex, Cs, j, ro, Pt)

    Mcre = Cb * Pex * (Cs * abs(j) + sqrt(j^2 + ro^2 * Pt/Pex))

end


function appendix2_3311(CorZ, t, ho, b, d, θc, E, μ, G, f1, f2, M1, M2, CurvatureSign, Lm, kϕ, Sf)

    θc=deg2rad(θc)
    bc=b-t/2-t/2*tan(θc/2)
    dc=d-t/2*tan(θc/2)

    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,xhf,yhf,yof=table2331(CorZ,t,bc,dc,θc)

    Lcrd, Ld = app23317(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf, Lm)

    kϕfe=app23133(E,Ixf,xof,hxf,Cwf,Ixyf,Iyf,L,G,Jf)

    kϕwe=app23335(E,t,μ,ho,L)

    kϕfg=app23315(L,Af,xof,hxf,Ixyf,Iyf,yof,Ixf)

    kϕwg=app23316(f1,f2,ho,t,L)

    Fcrd=app23311(kϕfe, kϕwe, kϕ, kϕfg, kϕwg)

    Mcrd=Sf*Fcrd

    return Mcrd

end





######

function app23331(CorZ, t, ho, b, d, θc, E, μ, G, f1, f2, M1, M2, CurvatureSign, Lm, kϕ, Sf)

    θc=deg2rad(θc)
    bc=b-t/2-t/2*tan(θc/2)
    dc=d-t/2*tan(θc/2)

    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof=table23131(CorZ,t,bc,dc,θc)

    L = app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

    β=app23333(L, Lm, M1, M2, CurvatureSign)

    kϕfe=app23133(E,Ixf,xof,hxf,Cwf,Ixyf,Iyf,L,G,Jf)

    kϕwe=app23335(E,t,μ,ho,L)

    kϕfg=app23135(L,Af,xof,hxf,Ixyf,Iyf,yof,Ixf)

    kϕwg=app23336(f1,f2,ho,t,L)

    Fcrd=app23332(β, kϕfe, kϕwe, kϕ, kϕfg, kϕwg)

    Mcrd=Sf*Fcrd

    return Mcrd

end


function app23333(L, Lm, M1, M2)

    β=1+0.4(L/Lm)^0.7*(1+M1/M2)^0.7

    if β>=1.3
        β=1.3
    end

    if β<=1.0
        β=1.0
    end

    return β

end







function app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

    Lcrd=(4*π^4*ho*(1-μ^2)/t^3*(Ixf*(xof-hxf)^2+Cwf-Ixyf^2/Iyf*(xof-hxf)^2)+π^4*ho^4/720)^(1/4)

    L=minimum([Lcrd, Lm])

    return Lcrd, L

end

function app23133(E,Ixf,xof,hxf,Cwf,Ixyf,Iyf,L,G,Jf)

    kϕfe=(π/L)^4*(E*Ixf*(xof-hxf)^2+E*Cwf-E*Ixyf^2/Iyf*(xof-hxf)^2)+(π/L)^2*G*Jf

end

function app23335(E,t,μ,ho,L)

    kϕwe=E*t^3/(12*(1-μ^2))*(3/ho+(π/L)^2*19/60*ho+(π/L)^4*ho^3/240)

end

function app23135(L,Af,xof,hxf,Ixyf,Iyf,yof,Ixf)

    kϕfg=(π/L)^2*(Af*((xof-hxf)^2*(Ixyf/Iyf)^2-2*yof*(xof-hxf)*(Ixyf/Iyf)+hxf^2+yof^2)+Ixf+Iyf)

end

function app23336(f1,f2,ho,t,L)

    ξweb=(f1-f2)/f1

    kϕwg=ho*t*π^2/13440*((((45360*(1-ξweb)+62160)*(L/ho)^2)+448*π^2+(ho/L)^2*(53+3*(1-ξweb))*π^4)/(π^4+28*π^2*(L/ho)^2+420*(L/ho)^4))

    return kϕwg

end

function app23332(β, kϕfe, kϕwe, kϕ, kϕfg, kϕwg)

    Fcrd=β*(kϕfe+kϕwe+kϕ)/(kϕfg+kϕwg)

end

function app2C2262(t, Lh, Lcrd)

    tr = t * (1 - Lh/Lcrd)^(1/3)

end


# ===========================================================================================
# AISI S100-16 Section A3.3.2 -- Strength Increase From Cold Work of Forming
# ===========================================================================================

"""
    a332(; Fyv, Fuv, R, t, C, θ=90.0)

AISI S100-16 Section A3.3.2, Eqs. A3.3.2-1 to A3.3.2-4. Average yield stress `Fya` of the full
unreduced section (compression members) or of the full flange (flexural members) from the
virgin yield stress `Fyv`, virgin tensile strength `Fuv`, inside bend radius `R`, thickness `t`,
and `C`, the ratio of total corner area to full area (of the section, or of the controlling
flange for flexural members -- see `a332_C`). `Fyf` is taken as the virgin yield stress
(no flat tests). Returns a NamedTuple with `Fya`, `Fyc`, `Bc`, `m`, `C`, `Fyf` and
`eq_A332_2_applies`, which is `false` when Eq. A3.3.2-2 is outside its limits
(Fuv/Fyv ≥ 1.2, R/t ≤ 7, included angle `θ` ≤ 120°).

Section A3.3.2 permits the increase only for sections not subject to strength reduction from
local or distortional buckling at stress level Fy (for beams, Mn = Mne from Section F3 and
Mnd = My from Section F4); see `a332_applicable_flexure`.
"""
function a332(;Fyv, Fuv, R, t, C, θ=90.0)

    Bc = 3.69 * (Fuv/Fyv) - 0.819 * (Fuv/Fyv)^2 - 1.79      #Eq. A3.3.2-3
    m = 0.192 * (Fuv/Fyv) - 0.068                            #Eq. A3.3.2-4
    Fyc = Bc * Fyv / (R/t)^m                                 #Eq. A3.3.2-2
    Fyf = Fyv
    Fya = min(C * Fyc + (1 - C) * Fyf, Fuv)                  #Eq. A3.3.2-1

    eq_A332_2_applies = (Fuv/Fyv >= 1.2) & (R/t <= 7) & (θ <= 120.0)

    return (Fya=Fya, Fyc=Fyc, Bc=Bc, m=m, C=C, Fyf=Fyf, eq_A332_2_applies=eq_A332_2_applies)

end

"""
    a332_corner_area(R, t; θ=90.0)

Cross-sectional area of one bent corner of included angle `θ` (degrees), inside radius `R`,
thickness `t`: arc length at the centerline radius (R + t/2) times `t`.
"""
a332_corner_area(R, t; θ=90.0) = deg2rad(θ) * (R + t/2) * t

"""
    a332_C(; w, t, R, n_corners=2, d=0.0, θ=90.0)

Ratio `C` of total corner cross-sectional area to full cross-sectional area (Section A3.3.2)
for a flat of width `w` plus `n_corners` corners, and optionally the adjacent flat `d` (e.g. a
stiffening lip). For the controlling flange of a lipped C-section flexural member use
`w` = flange flat width, `n_corners` = 2 (web/flange and flange/lip corners), and
`d` = 0 to count only the flange flat as the "full flange" (the convention that reproduces
the SFIA Technical Guide tables) or `d` = lip flat width to include the lip.
"""
function a332_C(;w, t, R, n_corners=2, d=0.0, θ=90.0)

    Ac = n_corners * a332_corner_area(R, t; θ=θ)
    C = Ac / ((w + d) * t + Ac)

    return C

end

"""
    a332_applicable_flexure(; Se, Sf, Mnd, My, rtol=1e-6)

Section A3.3.2 applicability for a flexural member: the cold work increase is permitted only
when the section is not subject to strength reduction from local or distortional buckling at
stress level Fy, i.e. the effective section modulus at Fy equals the gross section modulus
(Mn = Mne = Sf Fy from Section F3.1) and Mnd = My from Section F4.
"""
function a332_applicable_flexure(;Se, Sf, Mnd, My, rtol=1e-6)

    local_ok = Se >= Sf * (1 - rtol)
    distortional_ok = Mnd >= My * (1 - rtol)

    return local_ok & distortional_ok

end


# ===========================================================================================
# AISI S100-16 Appendix 1 -- Effective Width of Elements
# ===========================================================================================

"""
    appendix1_1_1(; w, t, f, k, E, μ)

Appendix 1 Section 1.1(a), uniformly compressed stiffened element, Eqs. 1.1-1 to 1.1-4.
`w` flat width, `t` thickness, `f` compressive stress in the element, `k` plate buckling
coefficient (4.0 for a stiffened element supported by a web on each edge). Returns
`(b, ρ, λ, Fcr, k)`. For serviceability (Section 1.1(b), Procedure I) call with `f = fd`.
"""
function appendix1_1_1(;w, t, f, k, E, μ)

    Fcr = k * π^2 * E / (12 * (1 - μ^2)) * (t/w)^2     #Eq. 1.1-4
    λ = sqrt(max(f, 0.0) / Fcr)                         #Eq. 1.1-3

    if λ <= 0.673
        ρ = 1.0
    else
        ρ = (1 - 0.22/λ) / λ                            #Eq. 1.1-2
    end

    b = ρ * w                                           #Eq. 1.1-1

    return (b=b, ρ=ρ, λ=λ, Fcr=Fcr, k=k)

end

"""
    appendix1_1_1_2(; w, t, f1, f2, ho, bo, E, μ)

Appendix 1 Section 1.1.2(a), webs and other stiffened elements under stress gradient,
Eqs. 1.1.2-1 to 1.1.2-10. Stresses are computed on the effective section, compression
positive: `f1` at the more compressed edge (f1 ≥ f2), `f2` at the other edge, tension negative.
`ho` out-to-out web depth, `bo` out-to-out compression flange width. Returns
`(b1, b2, be, k, ψ, ρ, λ, Fcr)` with `b1` adjacent to the f1 edge and `b2` adjacent to the
f2 edge (webs, Figure 1.1.2-1(a)) or to the neutral axis side. For webs (f2 in tension), Section
1.1.2 also requires b1 + b2 not to exceed the compression portion of the web, which the caller
checks because it depends on the effective neutral axis (see `f311_lipped_cee`).
"""
function appendix1_1_1_2(;w, t, f1, f2, ho, bo, E, μ)

    ψ = abs(f2 / f1)                                        #Eq. 1.1.2-1

    if f2 < 0.0   #(1) webs, f1 compression and f2 tension

        k = 4 + 2 * (1 + ψ)^3 + 2 * (1 + ψ)                 #Eq. 1.1.2-2
        be = appendix1_1_1(w=w, t=t, f=f1, k=k, E=E, μ=μ)

        if ho/bo <= 4
            b1 = be.b / (3 + ψ)                             #Eq. 1.1.2-3
            b2 = ψ > 0.236 ? be.b / 2 : be.b - b1           #Eqs. 1.1.2-4, 1.1.2-5
        else
            b1 = be.b / (3 + ψ)                             #Eq. 1.1.2-6
            b2 = be.b / (1 + ψ) - b1                        #Eq. 1.1.2-7
        end

    else   #(2) other stiffened elements, f1 and f2 both compression

        k = 4 + 2 * (1 - ψ)^3 + 2 * (1 - ψ)                 #Eq. 1.1.2-8
        be = appendix1_1_1(w=w, t=t, f=f1, k=k, E=E, μ=μ)
        b1 = be.b / (3 - ψ)                                 #Eq. 1.1.2-9
        b2 = be.b - b1                                      #Eq. 1.1.2-10

    end

    return (b1=b1, b2=b2, be=be.b, k=k, ψ=ψ, ρ=be.ρ, λ=be.λ, Fcr=be.Fcr)

end

"""
    appendix1_1_1_3(; w, t, f1, f2, ho, bo, dh, wc, E, μ, method=:appendix1_1_2_1)

Appendix 1 Section 1.1.3(a), C-section web with a hole under stress gradient (strength
determination). `w` = h, the flat web depth, `dh` the hole depth (hole centered at mid-depth),
`wc` the compression portion of the web flat measured from the compression flange corner to
the neutral axis, `f1`/`f2` the compression/tension stresses at the ends of the flat computed
on the effective section, `ho`/`bo` as in Section 1.1.2.

- dh/h < 0.38: `b1`, `b2` from Section 1.1.2(a) assuming no hole exists in the web
  (`hole_governs = false`).
- dh/h ≥ 0.38: the compression portion of the web between the flange corner and the hole edge
  (`wu` = min(wc, h/2 - dh/2)) is an unstiffened element (Section 1.2.1, k = 0.43) with f = f1;
  `b1` is its effective width measured from the flange corner, `b2` = 0, and the ineffective
  portion is adjacent to the hole (`hole_governs = true`). `method = :appendix1_1_2_2` instead
  treats that element as an unstiffened element under stress gradient (Section 1.2.2,
  Eq. 1.2.2-2, both edges in compression, stress decreasing toward the hole edge), a less
  conservative alternative that is not the wording of Section 1.1.3.

Warns when the geometry is outside the Section 1.1.3 limits (dh/h ≤ 0.7, h/t ≤ 200). The
serviceability determination of Section 1.1.3(b) uses Section 1.1.2 assuming no hole.
Returns `(b1, b2, k, ρ, λ, ψ, hole_governs, wu)`.
"""
function appendix1_1_1_3(;w, t, f1, f2, ho, bo, dh, wc, E, μ, method=:appendix1_1_2_1)

    dh / w > 0.7 && @warn "Section 1.1.3 applies for dh/h ≤ 0.7; dh/h = $(dh / w)"
    w / t > 200 && @warn "Section 1.1.3 applies for h/t ≤ 200; h/t = $(w / t)"

    if dh / w < 0.38

        e = appendix1_1_1_2(w=w, t=t, f1=f1, f2=f2, ho=ho, bo=bo, E=E, μ=μ)
        return (b1=e.b1, b2=e.b2, k=e.k, ρ=e.ρ, λ=e.λ, ψ=e.ψ, hole_governs=false, wu=NaN)

    else

        wu = min(wc, w / 2 - dh / 2)          #unstiffened element between the flange corner and the hole
        if method == :appendix1_1_2_1
            e = appendix1_1_2_1(w=wu, t=t, f=f1, E=E, μ=μ)                                #Section 1.1.3(a) as written, k = 0.43
        elseif method == :appendix1_1_2_2
            f2u = f1 * (wc - wu) / wc                                                   #compressive stress at the hole edge
            e = appendix1_1_2_2(w=wu, t=t, f1=f1, f2=f2u, E=E, μ=μ, f1_at_supported_edge=true)   #alternative: stress gradient, Eq. 1.2.2-2
        else
            error("method must be :appendix1_1_2_1 or :appendix1_1_2_2")
        end
        return (b1=e.b, b2=0.0, k=e.k, ρ=e.ρ, λ=e.λ, ψ=abs(f2 / f1), hole_governs=true, wu=wu)

    end

end

"""
    appendix1_1_2_1(; w, t, f, E, μ)

Appendix 1 Section 1.2.1(a), uniformly compressed unstiffened element: Section 1.1(a) with
k = 0.43. Returns `(b, ρ, λ, Fcr, k)`.
"""
appendix1_1_2_1(;w, t, f, E, μ) = appendix1_1_1(w=w, t=t, f=f, k=0.43, E=E, μ=μ)

"""
    appendix1_1_2_2(; w, t, f1, f2, E, μ, f1_at_supported_edge)

Appendix 1 Section 1.2.2(a), unstiffened elements and edge stiffeners with stress gradient,
Eqs. 1.2.2-1 to 1.2.2-7. `f1` is the maximum compressive stress on the element (positive),
`f2` the stress at the other longitudinal edge (compression positive, tension negative), and
`f1_at_supported_edge` says which edge carries f1. The four cases of the section are selected
from the sign of `f2` and the location of f1:

- both compression, f1 at the supported edge (stress decreases toward the unsupported edge,
  Figure 1.2.2-1(a)): k = 0.578/(ψ + 0.34), Eq. 1.2.2-2
- both compression, f1 at the unsupported edge (Figure 1.2.2-1(b)): Eq. 1.2.2-3
- f2 tension, unsupported edge in compression (Figure 1.2.2-2(a)): Eqs. 1.2.2-4, 1.2.2-5
- f2 tension, supported edge in compression (Figure 1.2.2-2(b)): Eqs. 1.2.2-6, 1.2.2-7

The effective width `b` is measured from the supported edge. Returns `(b, ρ, λ, k, ψ, Fcr)`.
"""
function appendix1_1_2_2(;w, t, f1, f2, E, μ, f1_at_supported_edge=true)

    ψ = abs(f2 / f1)                                        #Eq. 1.2.2-1

    if f2 >= 0.0   #(1) both edges in compression

        if f1_at_supported_edge
            k = 0.578 / (ψ + 0.34)                          #Eq. 1.2.2-2
        else
            k = 0.57 - 0.21 * ψ + 0.07 * ψ^2                #Eq. 1.2.2-3
        end

        e = appendix1_1_1(w=w, t=t, f=f1, k=k, E=E, μ=μ)
        ρ = e.ρ
        λ = e.λ
        Fcr = e.Fcr

    else   #(2) f1 compression, f2 tension

        if !f1_at_supported_edge   #(i) unsupported edge in compression

            k = 0.57 + 0.21 * ψ + 0.07 * ψ^2                #Eq. 1.2.2-5
            e = appendix1_1_1(w=w, t=t, f=f1, k=k, E=E, μ=μ)
            λ = e.λ
            Fcr = e.Fcr
            if λ <= 0.673 * (1 + ψ)
                ρ = 1.0
            else
                ρ = (1 + ψ) * (1 - 0.22 * (1 + ψ) / λ) / λ  #Eq. 1.2.2-4
            end

        else   #(ii) supported edge in compression

            if ψ < 1.0
                k = 1.70 + 5 * ψ + 17.1 * ψ^2               #Eq. 1.2.2-7
                e = appendix1_1_1(w=w, t=t, f=f1, k=k, E=E, μ=μ)
                λ = e.λ
                Fcr = e.Fcr
                if λ <= 0.673
                    ρ = 1.0
                else
                    ρ = (1 - ψ) * (1 - 0.22 / λ) / λ + ψ    #Eq. 1.2.2-6
                end
            else
                k = NaN
                λ = NaN
                Fcr = NaN
                ρ = 1.0
            end

        end

    end

    ρ = min(ρ, 1.0)
    b = ρ * w

    return (b=b, ρ=ρ, λ=λ, k=k, ψ=ψ, Fcr=Fcr)

end

"""
    appendix1_1_3(; w, t, f, D, d, θ, E, μ, ds_prime)

Appendix 1 Section 1.3(a), uniformly compressed element with a simple lip edge stiffener,
Eqs. 1.3-1 to 1.3-11 and Table 1.3-1. `w` flange flat width, `t` thickness, `f` stress in the
compression flange, `D` overall lip dimension and `d` lip flat width (Figure 1.3-1), `θ` lip
angle in degrees (140° ≥ θ ≥ 40°), `ds_prime` = d′s, the effective width of the lip from
Section 1.2.1 or 1.2.2. Returns `(b, b1, b2, ds, k, ρ, λ, S, Ia, Is, RI, n)`: `b1` is the
portion of the effective flange width adjacent to the web and `b2` the portion adjacent to
the lip; `ds` is the reduced effective width of the stiffener used for section properties.
"""
function appendix1_1_3(;w, t, f, D, d, θ, E, μ, ds_prime)

    (θ < 40.0 || θ > 140.0) && @warn "Section 1.3, Table 1.3-1 applies for 140° ≥ θ ≥ 40°; θ = $θ"
    D/w > 0.8 && @warn "Section 1.3, Table 1.3-1 applies for D/w ≤ 0.8; D/w = $(D/w)"

    S = 1.28 * sqrt(E / f)                                          #Eq. 1.3-7
    Is = d^3 * t * sind(θ)^2 / 12                                   #Eq. 1.3-10

    if w/t <= 0.328 * S

        Ia = 0.0                                                    #no edge stiffener needed
        RI = 1.0
        n = NaN
        k = 4.0
        ρ = 1.0
        λ = NaN
        b = w                                                       #Eq. 1.3-1
        b1 = w / 2                                                  #Eq. 1.3-2
        b2 = w / 2
        ds = ds_prime                                               #Eq. 1.3-3

    else

        Ia = min(399 * t^4 * ((w/t) / S - 0.328)^3, t^4 * (115 * (w/t) / S + 5))   #Eq. 1.3-8
        RI = min(Is / Ia, 1.0)                                      #Eq. 1.3-9
        n = max(0.582 - (w/t) / (4 * S), 1/3)                       #Eq. 1.3-11

        if D/w <= 0.25                                              #Table 1.3-1
            k = min(3.57 * RI^n + 0.43, 4.0)
        else
            k = min((4.82 - 5 * D/w) * RI^n + 0.43, 4.0)
        end

        e = appendix1_1_1(w=w, t=t, f=f, k=k, E=E, μ=μ)
        b = e.b
        ρ = e.ρ
        λ = e.λ
        b1 = (b / 2) * RI                                           #Eq. 1.3-4
        b2 = b - b1                                                 #Eq. 1.3-5
        ds = ds_prime * RI                                          #Eq. 1.3-6

    end

    return (b=b, b1=b1, b2=b2, ds=ds, k=k, ρ=ρ, λ=λ, S=S, Ia=Ia, Is=Is, RI=RI, n=n)

end


# ===========================================================================================
# AISI S100-16 Section F3.1 -- Effective Width Method, flexural strength for local buckling
# ===========================================================================================

"""
    f31(; Se, Fn, Set, Fy, design_code)

AISI S100-16 Section F3.1, Eq. F3.1-1: Mn = Se Fn ≤ Set Fy, where `Se` is the effective section
modulus at extreme fiber compressive stress `Fn` (Section F3.1.1, Appendix 1) and `Set` the
effective section modulus at extreme fiber tension stress Fy. Returns `(Mn, eMn)` with
Ω = 1.67, ϕ = 0.90 (LRFD), ϕ = 0.90 (LSD).
"""
function f31(;Se, Fn, Set, Fy, design_code)

    Ω = 1.67
    ϕ_LRFD = 0.90
    ϕ_LSD = 0.90

    Mn = min(Se * Fn, Set * Fy)                                     #Eq. F3.1-1

    eMn = calculate_factored_strength(Mn, Ω, ϕ_LRFD, ϕ_LSD, design_code)

    return Mn, eMn

end

"""
    f311_lipped_cee(; H, B, D, t, R, f, E, μ, θ=90.0, lip_method=:appendix1_1_2_2, dh=0.0, hole_method=:appendix1_1_2_1, tol=1e-8, maxiter=100)

AISI S100-16 Section F3.1.1 (members without holes): effective section properties of a lipped
C-section bending about its major axis (strong axis, compression flange at the top) at extreme
fiber compressive stress `f` (= Fn for strength, = fd for serviceability Procedure I), with the
effective width of every element from Appendix 1. `H` out-to-out web depth, `B` out-to-out
flange width, `D` overall lip length, `t` design thickness, `R` inside bend radius, `θ` lip
angle (degrees, 90 = perpendicular lip).

Geometry is the centerline "linear method" of the AISI Cold-Formed Steel Design Manual: flats
between the corner arcs (w = B - 2(R + t), h = H - 2(R + t), d = D - (R + t)) plus four 90°
corner arcs of centerline radius r = R + t/2, all of thickness `t`. Distances `y` are measured
from the extreme compression fiber (outside of the top flange).

Effective widths, iterating on the effective neutral axis location `ycg`:

- compression flange: Section 1.3 with `f` (the flange stress is taken conservatively as the
  extreme fiber stress) and the lip effective width d′s from `lip_method`:
  `:appendix1_1_2_2` (edge stiffener under stress gradient, Section 1.2.2, with f1 at the
  flange/lip corner and f2 at the lip tip computed on the effective section) or
  `:appendix1_1_2_1` (uniform compression at the f1 stress, k = 0.43);
- web: Section 1.1.2 with f1 and f2 at the ends of the web flat computed on the effective
  section, b1 + b2 limited to the compression portion of the web flat; with a web hole of depth
  `dh` > 0 centered at mid-depth (Section F3.1.2), Section 1.1.3 instead (`hole_method`, see
  `appendix1_1_1_3`), and the hole strip is removed from the section (net section through the
  hole, the SFIA Technical Guide convention for Sxe, Ma);
- tension flange, tension lip and the four corners: fully effective.

Returns a NamedTuple with `Se` (effective section modulus, extreme compression fiber),
`Set` (effective section modulus, extreme tension fiber), `Ie`, `Ae`, `ycg`, the gross
properties of the same centerline model (`A`, `Ix`, `Sx`), `fully_effective`, the
`flange`, `lip`, `web` element results, and `iterations`.
"""
function f311_lipped_cee(;H, B, D, t, R, f, E, μ, θ=90.0, lip_method=:appendix1_1_2_2, dh=0.0, hole_method=:appendix1_1_2_1, tol=1e-8, maxiter=100)

    r = R + t/2                    #centerline corner radius
    w = B - 2 * (R + t)            #flange flat
    h = H - 2 * (R + t)            #web flat
    d = D - (R + t)                #lip flat
    c = R + t                      #distance from an outside face to the start of the adjacent flat

    #90° corner arc: length, centroid offset from the arc center along each axis, own inertia
    Lc = π * r / 2
    ec = 2 * r / π
    Ic = r^3 * (π/4 - 2/π)

    #elements that never change: (length, y of centroid, own I about a horizontal centroidal axis)
    fixed = [
        (Lc, c - ec, Ic),                #top web/flange corner
        (Lc, c - ec, Ic),                #top flange/lip corner
        (Lc, H - c + ec, Ic),            #bottom web/flange corner
        (Lc, H - c + ec, Ic),            #bottom flange/lip corner
        (w, H - t/2, 0.0),               #tension flange flat
        (d, H - c - d/2, d^3/12),        #tension lip flat
    ]

    gross = [fixed; (w, t/2, 0.0); (d, c + d/2, d^3/12); (h, H/2, h^3/12)]
    A_g = t * sum(e[1] for e in gross)
    y_g = sum(e[1] * e[2] for e in gross) / sum(e[1] for e in gross)
    Ix_g = t * sum(e[1] * (e[2] - y_g)^2 + e[3] for e in gross)
    Sx_g = Ix_g / y_g

    ycg = y_g
    iterations = 0
    flange = lip = web = nothing
    elements = gross

    while true

        iterations += 1

        #compression lip: stresses at the flange/lip corner (f1) and at the lip tip (f2), effective section
        f1_lip = f * (ycg - c) / ycg
        f2_lip = f * (ycg - c - d) / ycg
        if lip_method == :appendix1_1_2_2
            lip = appendix1_1_2_2(w=d, t=t, f1=f1_lip, f2=f2_lip, E=E, μ=μ, f1_at_supported_edge=true)
        elseif lip_method == :appendix1_1_2_1
            lip = appendix1_1_2_1(w=d, t=t, f=f1_lip, E=E, μ=μ)
        else
            error("lip_method must be :appendix1_1_2_2 or :appendix1_1_2_1")
        end

        #compression flange with the lip edge stiffener
        flange = appendix1_1_3(w=w, t=t, f=f, D=D, d=d, θ=θ, E=E, μ=μ, ds_prime=lip.b)

        #web under stress gradient: f1 compression at the top of the flat, f2 tension at the bottom
        f1_web = f * (ycg - c) / ycg
        f2_web = -f * (H - ycg - c) / ycg
        wc = ycg - c                                     #compression portion of the web flat
        if dh > 0.0
            w_ = appendix1_1_1_3(w=h, t=t, f1=f1_web, f2=f2_web, ho=H, bo=B, dh=dh, wc=wc, E=E, μ=μ, method=hole_method)
            hole_governs = w_.hole_governs
        else
            e_ = appendix1_1_1_2(w=h, t=t, f1=f1_web, f2=f2_web, ho=H, bo=B, E=E, μ=μ)
            w_ = (b1=e_.b1, b2=e_.b2, k=e_.k, ρ=e_.ρ, λ=e_.λ, ψ=e_.ψ, hole_governs=false, wu=NaN)
            hole_governs = false
        end
        if hole_governs
            web_fully_effective = w_.b1 >= w_.wu
            ineffective_width = web_fully_effective ? 0.0 : w_.wu - w_.b1
        else
            web_fully_effective = (w_.b1 + w_.b2) >= wc
            ineffective_width = web_fully_effective ? 0.0 : wc - w_.b1 - w_.b2
        end
        web = (b1=w_.b1, b2=w_.b2, k=w_.k, ψ=w_.ψ, ρ=w_.ρ, λ=w_.λ,
               f1=f1_web, f2=f2_web, compression_portion=wc, fully_effective=web_fully_effective,
               ineffective_width=ineffective_width, dh=dh, hole_governs=hole_governs, wu=w_.wu)

        #effective elements
        elements = copy(fixed)
        push!(elements, (flange.b, t/2, 0.0))                              #compression flange flat, b = b1 + b2
        push!(elements, (flange.ds, c + flange.ds/2, flange.ds^3/12))      #compression lip, ds from the corner
        if dh > 0.0
            #net section through the hole: web flat above and below the hole
            y_hole_top = H/2 - dh/2
            y_hole_bot = H/2 + dh/2
            if hole_governs
                L1 = min(w_.b1, y_hole_top - c)                            #effective, from the flange corner toward the hole
                push!(elements, (L1, c + L1/2, L1^3/12))
                if ycg < y_hole_top                                        #neutral axis above the hole: tension part above the hole
                    L1t = y_hole_top - ycg
                    push!(elements, (L1t, ycg + L1t/2, L1t^3/12))
                end
            else
                Lw = y_hole_top - c
                if web_fully_effective || (w_.b1 + w_.b2) >= min(wc, Lw)
                    push!(elements, (Lw, c + Lw/2, Lw^3/12))
                else
                    L1 = w_.b1
                    push!(elements, (L1, c + L1/2, L1^3/12))
                    y_top = ycg - w_.b2
                    L2 = y_hole_top - y_top
                    push!(elements, (L2, y_top + L2/2, L2^3/12))
                end
            end
            L3 = (H - c) - y_hole_bot                                      #below the hole, fully effective
            push!(elements, (L3, y_hole_bot + L3/2, L3^3/12))
        elseif web_fully_effective
            push!(elements, (h, H/2, h^3/12))
        else
            L1 = w_.b1                                                     #adjacent to the compression flange
            push!(elements, (L1, c + L1/2, L1^3/12))
            y_top = ycg - w_.b2                                            #b2 adjacent to the neutral axis
            L2 = (H - c) - y_top                                           #b2 plus the tension portion
            push!(elements, (L2, y_top + L2/2, L2^3/12))
        end

        ycg_new = sum(e[1] * e[2] for e in elements) / sum(e[1] for e in elements)

        if abs(ycg_new - ycg) < tol || iterations >= maxiter
            ycg = ycg_new
            break
        end
        ycg = ycg_new

    end

    Ae = t * sum(e[1] for e in elements)
    Ie = t * sum(e[1] * (e[2] - ycg)^2 + e[3] for e in elements)
    Se = Ie / ycg
    Set = Ie / (H - ycg)

    fully_effective = (flange.ρ == 1.0) & (lip.ρ == 1.0) & web.fully_effective

    return (Se=Se, Set=Set, Ie=Ie, Ae=Ae, ycg=ycg, A=A_g, Ix=Ix_g, Sx=Sx_g, ycg_gross=y_g,
            fully_effective=fully_effective, flange=flange, lip=lip, web=web,
            w=w, h=h, d=d, r=r, f=f, dh=dh, lip_method=lip_method, iterations=iterations)

end

"""
    f311_lipped_cee_serviceability(; M, Fy, kwargs...)

Appendix 1 Section 1.1(b) Procedure I serviceability determination for a lipped C-section
in major-axis bending: the effective moment of inertia `Ie` at the service moment `M`, i.e.
at the extreme fiber stress `fd` for which `Se(fd) fd = M` (the SFIA Technical Guide
convention, where M is the allowable moment and the section is the unpunched section, i.e.
`dh` = 0; Section 1.1.3(b) likewise assumes no hole for serviceability). `fd` is found by bisection on (0, Fy] (`appendix1_1_1b_procedure_I`); if
Se(Fy) Fy ≤ M then fd = Fy. Remaining keyword arguments are those of `f311_lipped_cee`.
Returns `(Ie, Se, fd, section)` where `section` is the `f311_lipped_cee` result at fd.
"""
f311_lipped_cee_serviceability(;M, Fy, tol=1e-6, kwargs...) = appendix1_1_1b_procedure_I(f311_lipped_cee; M=M, Fy=Fy, tol=tol, kwargs...)

"""
    f311_unlipped_cee(; H, B, t, R, f, E, μ, tol=1e-8, maxiter=100)

AISI S100-16 Section F3.1.1 (members without holes): effective section properties of an
unlipped C-section (track) bending about its major axis (compression flange at the top) at
extreme fiber compressive stress `f` (= Fn for strength, = fd for serviceability Procedure I).
`H` out-to-out web depth, `B` out-to-out flange width, `t` design thickness, `R` inside bend
radius.

Geometry is the centerline linear model: flats w = B - (R + t) (one bend, free edge),
h = H - 2(R + t), plus two 90° corner arcs of centerline radius r = R + t/2. Distances `y` are
measured from the extreme compression fiber.

Effective widths, iterating on the effective neutral axis location `ycg`:

- compression flange: uniformly compressed unstiffened element, Section 1.2.1 (k = 0.43), at
  the extreme fiber stress `f` (conservative; the ineffective portion is at the free edge,
  which does not change the major-axis properties);
- web: Section 1.1.2 with f1 and f2 at the ends of the web flat computed on the effective
  section, b1 + b2 limited to the compression portion of the web flat;
- tension flange and the two corners: fully effective.

Returns a NamedTuple with `Se`, `Set`, `Ie`, `Ae`, `ycg`, the gross properties of the same
model (`A`, `Ix`, `Sx`), `fully_effective`, the `flange` and `web` element results, and
`iterations`.
"""
function f311_unlipped_cee(;H, B, t, R, f, E, μ, tol=1e-8, maxiter=100)

    r = R + t/2                    #centerline corner radius
    w = B - (R + t)                #flange flat, one bend and a free edge
    h = H - 2 * (R + t)            #web flat
    c = R + t                      #distance from an outside face to the start of the adjacent flat

    Lc = π * r / 2
    ec = 2 * r / π
    Ic = r^3 * (π/4 - 2/π)

    fixed = [
        (Lc, c - ec, Ic),                #top corner
        (Lc, H - c + ec, Ic),            #bottom corner
        (w, H - t/2, 0.0),               #tension flange flat
    ]

    gross = [fixed; (w, t/2, 0.0); (h, H/2, h^3/12)]
    A_g = t * sum(e[1] for e in gross)
    y_g = sum(e[1] * e[2] for e in gross) / sum(e[1] for e in gross)
    Ix_g = t * sum(e[1] * (e[2] - y_g)^2 + e[3] for e in gross)
    Sx_g = Ix_g / max(y_g, H - y_g)

    ycg = y_g
    iterations = 0
    flange = web = nothing
    elements = gross

    while true

        iterations += 1

        #compression flange, unstiffened under uniform compression
        flange = appendix1_1_2_1(w=w, t=t, f=f, E=E, μ=μ)

        #web under stress gradient
        f1_web = f * (ycg - c) / ycg
        f2_web = -f * (H - ycg - c) / ycg
        wc = ycg - c
        e_ = appendix1_1_1_2(w=h, t=t, f1=f1_web, f2=f2_web, ho=H, bo=B, E=E, μ=μ)
        web_fully_effective = (e_.b1 + e_.b2) >= wc
        ineffective_width = web_fully_effective ? 0.0 : wc - e_.b1 - e_.b2
        web = (b1=e_.b1, b2=e_.b2, k=e_.k, ψ=e_.ψ, ρ=e_.ρ, λ=e_.λ,
               f1=f1_web, f2=f2_web, compression_portion=wc, fully_effective=web_fully_effective,
               ineffective_width=ineffective_width)

        elements = copy(fixed)
        push!(elements, (flange.b, t/2, 0.0))                              #compression flange flat
        if web_fully_effective
            push!(elements, (h, H/2, h^3/12))
        else
            L1 = e_.b1                                                     #adjacent to the compression flange
            push!(elements, (L1, c + L1/2, L1^3/12))
            y_top = ycg - e_.b2                                            #b2 adjacent to the neutral axis
            L2 = (H - c) - y_top                                           #b2 plus the tension portion
            push!(elements, (L2, y_top + L2/2, L2^3/12))
        end

        ycg_new = sum(e[1] * e[2] for e in elements) / sum(e[1] for e in elements)

        if abs(ycg_new - ycg) < tol || iterations >= maxiter
            ycg = ycg_new
            break
        end
        ycg = ycg_new

    end

    Ae = t * sum(e[1] for e in elements)
    Ie = t * sum(e[1] * (e[2] - ycg)^2 + e[3] for e in elements)
    Se = Ie / ycg
    Set = Ie / (H - ycg)

    fully_effective = (flange.ρ == 1.0) & web.fully_effective

    return (Se=Se, Set=Set, Ie=Ie, Ae=Ae, ycg=ycg, A=A_g, Ix=Ix_g, Sx=Sx_g, ycg_gross=y_g,
            fully_effective=fully_effective, flange=flange, web=web,
            w=w, h=h, r=r, f=f, iterations=iterations)

end

"""
    appendix1_1_1b_procedure_I(section_fn; M, Fy, tol=1e-6, kwargs...)

Appendix 1 Section 1.1(b) Procedure I serviceability determination: the effective moment of
inertia `Ie` at the service moment `M`, i.e. at the extreme fiber stress `fd` for which
`Se(fd) fd = M`, where `section_fn(; f, kwargs...)` returns the effective section at stress
`f` (`f311_lipped_cee`, `f311_unlipped_cee`). `fd` is found by bisection on (0, Fy]; if
Se(Fy) Fy ≤ M then fd = Fy. Returns `(Ie, Se, fd, section)` where `section` is the
`section_fn` result at fd.
"""
function appendix1_1_1b_procedure_I(section_fn; M, Fy, tol=1e-6, kwargs...)

    g(fd) = section_fn(;f=fd, kwargs...).Se * fd - M

    if g(Fy) <= 0.0
        fd = Fy
    else
        lo = 1e-6 * Fy
        hi = Fy
        while (hi - lo) > tol * Fy
            mid = (lo + hi) / 2
            if g(mid) > 0.0
                hi = mid
            else
                lo = mid
            end
        end
        fd = (lo + hi) / 2
    end

    section = section_fn(;f=fd, kwargs...)

    return (Ie=section.Ie, Se=section.Se, fd=fd, section=section)

end

"""
    f311_unlipped_cee_serviceability(; M, Fy, kwargs...)

Appendix 1 Section 1.1(b) Procedure I serviceability determination for an unlipped C-section
(track) in major-axis bending at the service moment `M` (see `appendix1_1_1b_procedure_I`;
the SFIA Technical Guide convention takes M as the allowable moment). Remaining keyword
arguments are those of `f311_unlipped_cee`. Returns `(Ie, Se, fd, section)`.
"""
f311_unlipped_cee_serviceability(;M, Fy, tol=1e-6, kwargs...) = appendix1_1_1b_procedure_I(f311_unlipped_cee; M=M, Fy=Fy, tol=tol, kwargs...)



# ===========================================================================================
# AISI S100-16 Appendix 2, Section 2.3.3.3 -- Distortional buckling of flexural members
# ===========================================================================================

"""
    appendix2_2_3_3_3(; t, ho, b, d, θ, E, μ, Sf, f1, f2, kϕ=0.0, Lm=Inf, M1=nothing, M2=nothing, CorZ=0, t_web=t)

AISI S100-16 Appendix 2 Section 2.3.3.3, elastic distortional buckling moment of a C- or
Z-section with a single web and a simple lip edge stiffener, Eqs. 2.3.3.3-1 to 2.3.3.3-6 with
the flange properties of Table 2.3.1.3-1 (`table_2_3_3__1`) and the flange stiffnesses of
Eqs. 2.3.1.3-3 and 2.3.1.3-5. `b` and `d` are centerline flange width and lip length
(b = B - t, d = D - t/2 for a 90° lip), `ho` out-to-out web depth, `θ` lip angle in
degrees, `Sf` gross section modulus to the extreme compression fiber, `f1`/`f2` the web end
stresses on the gross section (pure symmetric bending: f1 = -f2, ξweb = 2), `kϕ` restraint
stiffness (0 = unrestrained), `Lm` distance between discrete distortional restraints, and
`M1`/`M2` the end moments for the moment gradient factor β (β = 1.0 when omitted). `t_web`
is the web thickness used in Eqs. 2.3.3.3-5 and 2.3.3.3-6 (the modified thickness tr of
Section 2.3.4.3 for members with web holes, see `appendix2_2_3_4_3`).
Returns `(Mcrd, Fcrd, Lcrd, L, β, kϕfe, kϕwe, kϕfg, kϕwg, flange, t_web)`.
"""
function appendix2_2_3_3_3(;t, ho, b, d, θ, E, μ, Sf, f1, f2, kϕ=0.0, Lm=Inf, M1=nothing, M2=nothing, CorZ=0, t_web=t)

    G = E / (2 * (1 + μ))

    Af, Jf, Ixf, Iyf, Ixyf, Cwf, xof, hxf, hyf, yof = table_2_3_3__1(CorZ, t, b, d, θ)   #Table 2.3.1.3-1

    Lcrd = appendix2_2_3_3_2__4(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf)                #Eq. 2.3.3.3-4
    L = min(Lcrd, Lm)

    if M1 === nothing
        β = 1.0
    else
        β = min(max(1 + 0.4 * (L / Lm)^0.7 * (1 + M1 / M2)^0.7, 1.0), 1.3)               #Eq. 2.3.3.3-3
    end

    kϕfe = app23133(E, Ixf, xof, hxf, Cwf, Ixyf, Iyf, L, G, Jf)                          #Eq. 2.3.1.3-3
    kϕwe = app23335(E, t_web, μ, ho, L)                                                  #Eq. 2.3.3.3-5
    kϕfg = app23135(L, Af, xof, hxf, Ixyf, Iyf, yof, Ixf)                                #Eq. 2.3.1.3-5
    kϕwg = app23336(f1, f2, ho, t_web, L)                                                #Eq. 2.3.3.3-6

    Fcrd = app23332(β, kϕfe, kϕwe, kϕ, kϕfg, kϕwg)                                      #Eq. 2.3.3.3-2
    Mcrd = Sf * Fcrd                                                                     #Eq. 2.3.3.3-1

    flange = (Af=Af, Jf=Jf, Ixf=Ixf, Iyf=Iyf, Ixyf=Ixyf, Cwf=Cwf, xof=xof, hxf=hxf, hyf=hyf, yof=yof)

    return (Mcrd=Mcrd, Fcrd=Fcrd, Lcrd=Lcrd, L=L, β=β, kϕfe=kϕfe, kϕwe=kϕwe, kϕfg=kϕfg, kϕwg=kϕwg, flange=flange, t_web=t_web)

end

"""
    appendix2_2_3_4_3(; Lh, t, kwargs...)

AISI S100-16 Appendix 2 Section 2.3.4.3, distortional buckling moment of a flexural member with
a web hole of length `Lh`: Section 2.3.3.3 (`appendix2_2_3_3_3`, same keyword arguments) with
the web thickness in Eqs. 2.3.3.3-5 and 2.3.3.3-6 replaced by tr = t (1 - Lh/Lcrd)^(1/3)
(Eq. 2.3.4.3-1), where Lcrd is the gross-section half-wavelength of Eq. 2.3.3.3-4. Returns the
`appendix2_2_3_3_3` NamedTuple plus `tr` and `Lh`.
"""
function appendix2_2_3_4_3(;Lh, t, kwargs...)

    gross = appendix2_2_3_3_3(;t=t, kwargs...)
    tr = app2C2262(t, Lh, gross.Lcrd)                                                     #Eq. 2.3.4.3-1
    net = appendix2_2_3_3_3(;t=t, t_web=tr, kwargs...)

    return merge(net, (tr=tr, Lh=Lh, Mcrd_gross=gross.Mcrd))

end


end #module
