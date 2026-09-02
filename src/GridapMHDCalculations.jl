module GridapMHDCalculations

using Revise

include("../../master/GridapMHD.jl/src/GridapMHD.jl")       #TDB: Add a global path or something along those lines
using .GridapMHD: main, default_solver_params, uses_petsc     #.GridapMHD = GridapMHDCalculations.GridapMHD This is necessary to get the exports from GridapMHD

using PartitionedArrays
using SparseArrays
using SparseMatricesCSR

using Gridap
using Gridap.Helpers, Gridap.Algebra, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField, Gridap.ODEs
using Gridap.Fields

using GridapDistributed
using GridapGmsh
using GridapPETSc
using GridapP4est
using GridapSolvers

using FileIO
using BSON

include("Data.jl")
export output_info

include("models/models.jl")

include("SteadyState.jl")
export SteadyState

include("MagneticFields.jl")
include("PostProcess.jl")
export exec_post_process

include("utils.jl")

#Tests
include("../test/serial/serialTest.jl")

end # module GridapMHDCalculations
