using AISIS100

#D100-17E AISI Cold-Formed Steel Design Manual Volume 1
# pg II-184
#pg II-15
# 800S200-54

# CorZ = 0
# t = 0.0566
# d = 0.625 - t/2
# b = 2.000 - t
# θ = 90.0

ho = 8.0
μ = 0.30

CorZ = 0
t = 0.0566
b = 2.00 - t
d = 0.625 - t/2
θ = 90.0
Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,xhf,yhf,yof = AISIS100.v16S3.table_2_3_3__1(CorZ,t,b,d,θ)

Lcrd = AISIS100.v16S3.appendix2_2_3_3_2__4(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)






ho = 8.0
μ = 0.30
t = 0.058
Ixf = 0.006997521335553253
xof = 1.171860817763336
xhf = -1.770139182236664
yhf = -0.07636081776333603
yof = -0.07636081776333603
Cwf = 0.0
Lcrd = AISIS100.v16S3.appendix2_2_3_3_1__7(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)
