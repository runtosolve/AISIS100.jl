module AISIS100


#Hooks for optional dependencies. They have no methods here; the package extensions add them when the
#optional package is loaded, so `using AISIS100` alone does not load CSV, DataFrames or Unitful.

"""
    inches(x), millimeters(x)

Length `x` as a Unitful quantity. Defined by the Unitful extension (ext/AISIS100UnitfulExt.jl); the
functions that use them (v16.j25, v16.tablej3311, v16S3.tablej3311) take Unitful thicknesses, so callers
already have Unitful loaded.
"""
function inches end
function millimeters end

"""
    read_g5_table(filename)

Read an AISI S100-16 Table G5 CSV file into a DataFrame. Defined by the CSV + DataFrames extension
(ext/AISIS100TablesExt.jl); `table_g52()` and `table_g53()` need `using CSV, DataFrames`.
"""
function read_g5_table end



export v16
include("v16.jl")
using .v16

export v16S3
include("v16S3.jl")
using .v16S3

export v2024
include("v2024.jl")
using .v2024

end # module
