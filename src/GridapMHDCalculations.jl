module GridapMHDCalculations


include("../../master/GridapMHD.jl/src/GridapMHD.jl")
using .GridapMHD: main, default_solver_params, uses_petsc     #.GridapMHD = GridapMHDCalculations.GridapMHD

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



include("models/models.jl")

include("SteadyState.jl")
export SteadyState

include("MagneticFields.jl")
include("PostProcess.jl")
include("utils.jl")

end # module GridapMHDCalculations
