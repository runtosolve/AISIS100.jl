module AISIS100TablesExt

using AISIS100, CSV, DataFrames

function AISIS100.read_g5_table(filename)

    data = CSV.File(filename)

    table = DataFrame(support_condition = Vector(data.support_condition), flange_condition = Vector(data.flange_condition), load_case = Vector(data.load_case), load_location=Vector(data.load_location), C=Vector(data.C), C_R = Vector(data.C_R), C_N = Vector(data.C_N), C_h = Vector(data.C_h), ASD = Vector(data.ASD), LRFD = Vector(data.LRFD), LSD = Vector(data.LSD), limits = Vector(data.limits))

    return table

end

end # module
