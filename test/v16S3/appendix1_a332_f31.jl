# AISI S100-16 Appendix 1 (effective width), Section A3.3.2 (cold work), Section F3.1 (EWM
# flexural strength), Section F4.1 and Appendix 2 Sections 2.3.3.3 / 2.3.4.3, checked against
# hand calculations and the SFIA Technical Guide (FINAL 2026-01-07) section property tables.
using AISIS100, Test
const S = AISIS100.v16S3
E = 29500.0
μ = 0.30

@testset "Appendix 1 Section 1.1, stiffened element" begin
    # w/t = 60, k = 4, f = 33 ksi: Fcr = 4 π² E /(12(1-μ²)) (1/60)², λ = sqrt(f/Fcr)
    e = S.appendix1_1_1(w=60.0, t=1.0, f=33.0, k=4.0, E=E, μ=μ)
    Fcr = 4 * π^2 * E / (12 * (1 - μ^2)) / 60^2
    @test e.Fcr ≈ Fcr
    @test e.λ ≈ sqrt(33.0 / Fcr)
    @test e.ρ ≈ (1 - 0.22 / e.λ) / e.λ
    @test e.b ≈ e.ρ * 60.0
    e2 = S.appendix1_1_1(w=20.0, t=1.0, f=33.0, k=4.0, E=E, μ=μ)
    @test e2.ρ == 1.0 && e2.b == 20.0        # λ ≤ 0.673, fully effective
end

@testset "Appendix 1 Section 1.1.2, web under stress gradient" begin
    w_ = S.appendix1_1_1_2(w=3.0, t=0.0346, f1=30.0, f2=-30.0, ho=3.625, bo=1.625, E=E, μ=μ)
    @test w_.ψ ≈ 1.0
    @test w_.k ≈ 4 + 2 * 2^3 + 2 * 2                    # 24 for ψ = 1
    @test w_.b1 ≈ w_.be / 4                             # ho/bo ≤ 4: b1 = be/(3 + ψ)
    @test w_.b2 ≈ w_.be / 2                             # ψ > 0.236
end

@testset "Section A3.3.2, cold work of forming" begin
    cw = S.a332(Fyv=50.0, Fuv=65.0, R=0.107, t=0.0713, C=0.3)
    @test cw.Bc ≈ 3.69 * 1.3 - 0.819 * 1.3^2 - 1.79
    @test cw.m ≈ 0.192 * 1.3 - 0.068
    @test cw.Fyc ≈ cw.Bc * 50.0 / (0.107 / 0.0713)^cw.m
    @test cw.Fya ≈ min(0.3 * cw.Fyc + 0.7 * 50.0, 65.0)
    @test cw.eq_A332_2_applies
    @test !S.a332(Fyv=50.0, Fuv=55.0, R=0.107, t=0.0713, C=0.3).eq_A332_2_applies   # Fu/Fy < 1.2
    @test S.a332_corner_area(0.107, 0.0713) ≈ π/2 * (0.107 + 0.0713/2) * 0.0713
end

@testset "Section F4.1, Eq. F4.1-2 exponent 0.5" begin
    My = 10.0; Mcrd = 15.0
    Mnd, = S.f411(My, Mcrd, "nominal")
    @test Mnd ≈ (1 - 0.22 * (Mcrd / My)^0.5) * (Mcrd / My)^0.5 * My
    @test S.f411(10.0, 40.0, "nominal")[1] == 10.0       # λd ≤ 0.673
end

@testset "362S162-33 vs SFIA Technical Guide" begin
    H = 3.625; B = 1.625; D = 0.5; t = 0.0346; R = 0.0765; Fy = 33.0
    g = S.f311_lipped_cee(H=H, B=B, D=D, t=t, R=R, f=Fy, E=E, μ=μ)
    @test isapprox(g.A, 0.262; atol=0.001)               # gross A
    @test isapprox(g.Ix, 0.551; atol=0.001)              # gross Ix
    @test isapprox(g.Sx, 0.304; atol=0.001)              # gross Sx
    @test !g.fully_effective && g.flange.ρ < 1.0 && g.web.fully_effective
    # SFIA Sxe, Ma, Mad include the standard 1.5 x 4 in. punchout (Appendix 1 Section 1.1.3 with the
    # unstiffened web element treated with the Section 1.2.2 stress gradient; Appendix 2 Section 2.3.4.3)
    n = S.f311_lipped_cee(H=H, B=B, D=D, t=t, R=R, f=Fy, E=E, μ=μ, dh=1.5, hole_method=:appendix1_1_2_2)
    @test n.web.hole_governs
    @test isapprox(n.Se, 0.275; atol=0.001)              # SFIA Sxe
    Mn, Ma = S.f31(Se=n.Se, Fn=Fy, Set=n.Set, Fy=Fy, design_code="ASD")
    @test isapprox(Ma, 5.44; atol=0.02)                  # SFIA Ma
    dist = S.appendix2_2_3_4_3(Lh=4.0, t=t, ho=H, b=B - t, d=D - t/2, θ=90.0, E=E, μ=μ, Sf=g.Sx, f1=Fy, f2=-Fy)
    @test dist.tr ≈ t * (1 - 4.0 / dist.Lcrd)^(1/3)
    Mnd, Mad = S.f411(g.Sx * Fy, dist.Mcrd, "ASD")
    @test isapprox(Mad, 5.20; atol=0.02)                 # SFIA Mad
    serv = S.f311_lipped_cee_serviceability(M=Ma, Fy=Fy, H=H, B=B, D=D, t=t, R=R, E=E, μ=μ)
    @test isapprox(serv.Ie, 0.551; atol=0.002)           # SFIA Ixe (unpunched, Procedure I at Ma)
    # no-hole strength for comparison: larger than the punched value
    Mn_g, Ma_g = S.f31(Se=g.Se, Fn=Fy, Set=g.Set, Fy=Fy, design_code="ASD")
    @test Ma_g > Ma
end

@testset "250S137-68 (50 ksi) vs SFIA: fully effective, cold work applied" begin
    H = 2.5; B = 1.375; D = 0.375; t = 0.0713; R = 0.107; Fy = 50.0; Fu = 65.0
    g = S.f311_lipped_cee(H=H, B=B, D=D, t=t, R=R, f=Fy, E=E, μ=μ)
    @test g.flange.ρ == 1.0 && g.lip.ρ == 1.0 && g.web.fully_effective
    C = S.a332_C(w=g.w, t=t, R=R, n_corners=2)          # corners / (flange flat + corners)
    cw = S.a332(Fyv=Fy, Fuv=Fu, R=R, t=t, C=C)
    @test isapprox(cw.Fya, 57.8; atol=0.1)               # back-calculated from the SFIA table: 10.62 x 1.67 / 0.307
    n = S.f311_lipped_cee(H=H, B=B, D=D, t=t, R=R, f=cw.Fya, E=E, μ=μ, dh=0.75)
    @test !n.web.hole_governs                            # dh/h < 0.38
    @test isapprox(n.Se, 0.307; atol=0.0015)             # 0.3060 at Fya (SFIA Sxe is tabulated at Fy)
    Mn, Ma = S.f31(Se=n.Se, Fn=cw.Fya, Set=n.Set, Fy=cw.Fya, design_code="ASD")
    @test isapprox(Ma, 10.62; rtol=0.005)
end

@testset "Track (unlipped C) vs SFIA Technical Guide: Sxe at Fy, Ma, Ixe by Procedure I at Ma" begin
    # SFIA track out-to-out depth H = nominal inside depth + 2t + R (Section Properties Table note 8)
    for (H_in, B, td, t, R, Fy, Sxe_s, Ma_s, Ixe_s) in ((3.625, 1.25, 33, 0.0346, 0.0765, 33.0, 0.174, 3.44, 0.385),
                                                       (6.0, 2.0, 54, 0.0566, 0.0849, 50.0, 0.718, 21.51, 2.643),
                                                       (12.0, 2.0, 97, 0.1017, 0.1526, 33.0, 4.303, 85.02, 29.808))
        H = H_in + 2 * t + R
        g = S.f311_unlipped_cee(H=H, B=B, t=t, R=R, f=Fy, E=E, μ=μ)
        @test isapprox(g.Se, Sxe_s; rtol=0.006)
        Mn, Ma = S.f31(Se=g.Se, Fn=Fy, Set=g.Set, Fy=Fy, design_code="ASD")
        @test isapprox(Ma, Ma_s; rtol=0.006)
        serv = S.f311_unlipped_cee_serviceability(M=Ma, Fy=Fy, H=H, B=B, t=t, R=R, E=E, μ=μ)
        @test isapprox(serv.Ie, Ixe_s; rtol=0.006)
        @test serv.fd < Fy && isapprox(serv.section.Se * serv.fd, Ma; rtol=1e-4)
    end
    # fully effective section: Ie = Ix and fd = Fy when Se(Fy) Fy <= M
    H = 1.625 + 2 * 0.1017 + 0.1526                     # 162T125-97: flange w/t = 9.8, λ < 0.673 at 33 ksi
    g = S.f311_unlipped_cee(H=H, B=1.25, t=0.1017, R=0.1526, f=33.0, E=E, μ=μ)
    @test g.fully_effective && g.Ie ≈ g.Ix
    serv = S.f311_unlipped_cee_serviceability(M=g.Se * 33.0, Fy=33.0, H=H, B=1.25, t=0.1017, R=0.1526, E=E, μ=μ)
    @test serv.fd == 33.0 && serv.Ie ≈ g.Ix
end
