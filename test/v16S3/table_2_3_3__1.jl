using AISIS100

#D100-17E AISI Cold-Formed Steel Design Manual Volume 1
# pg II-184
# 800S200-54

CorZ = 0
t = 0.0566
d = 0.597
b = 1.943
θ = 90.0
Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,xhf,yhf,yof = AISIS100.v16S3.table_2_3_3__1(CorZ,t,b,d,θ)


isapprox(Af, 0.144, rtol = 0.03)
isapprox(Ixf, 0.00334, rtol = 0.03)
isapprox(Iyf, 0.0590, rtol = 0.03)
isapprox(Ixyf, 0.00750, rtol = 0.03)
isapprox(xof, 0.743, rtol = 0.03)
isapprox(yof, -0.0702, rtol = 0.03)
isapprox(xhf, -1.20, rtol = 0.03)
isapprox(Jf, 0.000154, rtol = 0.03)
isapprox(Cwf, 0.0, rtol = 0.03)