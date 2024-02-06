module v16S3

using CSV, DataFrames, Unitful, Statistics, LinearAlgebra


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
        Mnd=(1-0.22*(Mcrd/My)^0.4)*(Mcrd/My)^0.4*My
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

function j25(;L, t1, Fu1, t2, Fu2, tw, Fxx, loading_direction, design_code)

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

    if (t > 0.10u"inch") | (t > 2.54u"mm")

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

    if design_code == "AISI S100-16 ASD"
        ePnv  = Pnv / Ω
    elseif design_code == "AISI S100-16 LRFD"
        ePnv = Pnv * ϕ_LRFD
    elseif design_code == "AISI S100-16 LSD"
        ePnv = Pnv * ϕ_LSD
    elseif design_code == "nominal"
        ePnv = Pnv
    end


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

        Ixf=t*(t^2*b^2+4*b*d^3+t^2*b*d+t^2*b*d+d^4)/(12*(b+d))
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

end #module
