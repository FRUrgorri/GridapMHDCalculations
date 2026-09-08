module serialTest

using GridapMHDCalculations
using GridapMHDCalculations: kp_shercliff_cartesian, u_parabolic, pp_gradp_check, pp_Noslip_check, post_process_basic 
using GridapMHDCalculations.models: map_Roberts, channel_model
using Gridap
using SparseArrays

include("mpi_channel.jl")
include("mpi_channel_multigrid.jl")

end