module mpiTest

using GridapMHDCalculations
using GridapMHDCalculations: kp_shercliff_cartesian, u_parabolic, gradp_check, noSlip_check, writeFields_vtk 
using GridapMHDCalculations.models
using GridapMHDCalculations.models: insulated_channel, channel_geom, channel_mesh
using Gridap, GridapPETSc
using SparseArrays, SparseMatricesCSR
using MPI

#Add new methods to the equivalent test in serial
include("mpi_channel.jl")
include("mpi_channel_multigrid.jl")

end
