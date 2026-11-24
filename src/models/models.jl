module models

using Gridap
using GridapSolvers
using GridapDistributed
using PartitionedArrays
using GridapMHDCalculations: mounted_models
using GridapMHDCalculations
using Gmsh: gmsh
using GridapGmsh

include("BC.jl")
export BC_tags, BC_values, BC

include("maps.jl")
export map_Roberts

include("gmsh_models.jl")
export gmsh_model

include("channel_models.jl")
include("expansion_models.jl")

#include("deprecated/channel_model.jl")  #Deprecated

end 