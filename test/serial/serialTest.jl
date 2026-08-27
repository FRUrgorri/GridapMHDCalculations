module serialTest

using GridapMHDCalculations
using GridapMHDCalculations.models: map_Roberts, channel_model
using Gridap
using SparseArrays

include("channel.jl")
include("channel_multigrid.jl")

end 