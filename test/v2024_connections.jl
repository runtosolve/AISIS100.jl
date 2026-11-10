# AISI S100-2024 Chapters D and J functions added in v0.4.2, checked by hand calculation
using AISIS100, Test
const inch = AISIS100.v16.Unitful.uparse("inch")

v = AISIS100.v2024

@testset "S100-2024 D2, D3" begin
    @test v.d2__1(1.0, 50.0, "ASD").eRn ≈ 50.0 / 1.67
    @test v.d3__1(1.0, 62.0, "LRFD").eRn ≈ 0.75 * 62.0
end

@testset "S100-2024 J2.6 fillet weld" begin
    # L/t = 2.5/0.105 < 25, t > 0.10 in.: Pnv1 = (1 - 0.01*2.5/0.105)*2.5*0.105*62 = 12.4, Pn = 0.75*0.088375*2.5*60 = 9.942
    r = v.j2_6__1to7(2.5, 0.105, 62.0, 1.0, 52.0, 0.707 * 0.125, 60.0, "longitudinal", "ASD")
    @test r.Pnv1 ≈ 12.4
    @test r.Rn ≈ 0.75 * 0.707 * 0.125 * 2.5 * 60.0
    @test r.eRn ≈ r.Rn / 2.55
    # L/t ≥ 25 and t > 0.10: Pnv carries Ω = 3.05, Pn Ω = 2.55 -> lesser available strength
    r = v.j2_6__1to7(4.0, 0.105, 62.0, 0.105, 62.0, 0.2, 70.0, "longitudinal", "ASD")
    @test r.eRn ≈ min(0.75 * 0.105 * 4.0 * 62.0 / 3.05, 0.75 * 0.2 * 4.0 * 70.0 / 2.55)
end

@testset "S100-2024 J3.3.1 bearing, Tables J3.3.1-1 and -2" begin
    @test v.table_j3_3_1__1(0.375, 0.061, "standard hole") == 3.0
    @test v.table_j3_3_1__1(0.5, 0.030, "standard hole") ≈ 4.0 - 0.1 * 0.5 / 0.030
    @test v.table_j3_3_1__1(0.75, 0.030, "standard hole") == 1.8
    @test v.table_j3_3_1__1(0.375, 0.040, "oversized hole") ≈ 1.0 + 14.0 / (0.375 / 0.040)
    @test v.table_j3_3_1__2("single shear", "standard hole", "no") == 0.75
    @test v.table_j3_3_1__2("inside sheet of double shear", "standard hole", "no") == 1.33
    @test v.j3_3_1__1(3.0, 0.75, 0.375, 0.061, 62.0, "ASD").eRn ≈ 3.0 * 0.75 * 0.375 * 0.061 * 62.0 / 2.50
end

@testset "S100-2024 App. A J3.4 and J6.1" begin
    Ab = π * 0.375^2 / 4
    @test v.appendixA_j3_4__1(Ab, 48.0, "ASD").eRn ≈ Ab * 48.0 / 2.0
    @test v.j6_1__2(1, 0.061, 0.59) ≈ 0.07198
    @test v.j6_1__1(62.0, 0.07198, "bolts", "ASD").eRn ≈ 0.6 * 62.0 * 0.07198 / 2.00     # S100-16: Ω = 2.22
    @test v.j6_1__1(62.0, 0.07198, "bolts", "LRFD").eRn ≈ 0.75 * 0.6 * 62.0 * 0.07198    # S100-16: ϕ = 0.65
    @test v.j6_1__1(62.0, 0.07198, "welds", "ASD").eRn ≈ 0.6 * 62.0 * 0.07198 / 2.50
end

@testset "S100-16 Table J3.3.1-1 oversized holes" begin
    @test AISIS100.v16.tablej3311(d=0.375inch, t=0.040inch, hole_shape="oversized hole") ≈ 1.0 + 14.0 / (0.375 / 0.040)
end
