using AISIS100

My = 66.3
Mcrd = 47.1
My_net = 60.7
design_code = "AISI S100-16 ASD"


λd = sqrt(My/Mcrd)
λd1 = 0.673*(My_net/My)^3
λd2 = 0.673*(1.7*(My/My_net)^(2.7)-0.7)

Mnd, eMnd = v16.f421_6(My, Mcrd, My_net, design_code)