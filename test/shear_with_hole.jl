using S100AISI


d_h = 30.0
L_h = 30.0
h = 108.1
Vy = 34.275
a = h  #supposed to be twice the distance to end of the section
b_f = 50.0
E = 200000.0
μ = 0.30
t = 1.95

A_h, d_h_eq, L_h_eq = S100AISI.v24.g39_10(d_h, L_h)  

V_yh, a0, a1, a2 = S100AISI.v24.g31_5(d_h_eq, L_h_eq, h, Vy)

kv = S100AISI.v24.g38(a, h, d_h_eq, b_f)

Fcr = S100AISI.v16.g232(E, μ, kv, h, t)

A_w = h * t

Vcr = Fcr * A_w 


α_vh, V_crh = S100AISI.v24.g36_7(L_h_eq, d_h_eq, h, Vcr)


Vn, eVn = S100AISI.v24.g211_212(V_crh, V_yh, "nominal")