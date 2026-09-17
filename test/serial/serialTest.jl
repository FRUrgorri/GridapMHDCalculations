module serialTest

using GridapMHDCalculations
using GridapMHDCalculations: kp_shercliff_cartesian, u_parabolic, gradp_check, noSlip_check, writeFields_vtk 
using GridapMHDCalculations.models
using GridapMHDCalculations.models: insulated_channel, channel_geom, channel_mesh
using Gridap
using SparseArrays

include("channel.jl")
include("channel_blocks.jl")
include("channel_multigrid.jl")

end 
