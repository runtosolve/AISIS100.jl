using AISIS100

M_test_t25 = [18053212, 18068901, 18348153]  #N-mm, t=2.5mm tests
M_test_t30 = [24899589, 24899589, 24973628]  #N-mm, t=3.0mm tests
M_test_t40 = [33497101, 34984155, 34604900]  #$N-mm, t=4.0mm tests

M_FEA_t25 = 17400542 #N-mm, t=2.5mm 
M_FEA_t30 = 23508658 #N-mm, t=3mm 
M_FEA_t40 = 32959452 #$N-mm, t=4.0mm 

Rt = [M_test_t25; M_test_t30; M_test_t40]

Rn = [ones(Float64, 3)* M_FEA_t25; ones(Float64, 3)* M_FEA_t30; ones(Float64, 3)* M_FEA_t40;]

n = length(Rt)

Cc = AISIS100.v16.k2117(n, Rt, Rn)