module models

using Gridap
using GridapSolvers
using GridapDistributed
using PartitionedArrays
using GridapMHDCalculations: mounted_models

include("BC.jl")
export BC_tags, BC_values, BC

include("maps.jl")
export map_Roberts


include("channel_models.jl")

#include("deprecated/channel_model.jl")  #Deprecated

end 