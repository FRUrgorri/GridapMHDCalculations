using GridapMHDCalculations
using Gridap
using GridapPETSc
using SparseArrays, SparseMatricesCSR

#Communication variables

#@assert isfile("../Transfer.jl")  #Investigating why this assertion fails
include("../Transfer.jl")


#Saving Path
_path = "../results/channel_Ha10Re1_H1H1_direct"
#@assert isdir(_path)

#Inputs of the calculation
b = 1.0         #Channel aspect ratio
L = 6.0         #Channel lenght ratio
Ha = 10         #Hartmann number
Re = 1          #Reynolds number

#Mesh inputs
nX = 16
nY = 16
nZ = 30

#Define the boundary fields
U_inlet((x,y,z))=VectorValue(0.0,0.0,GridapMHDCalculations.u_parabolic(b)(x,y))
B((x,y,z))=VectorValue(0.0,1.0,0.0)

#Define the Gridap models
map = GridapMHDCalculations.models.map_Roberts(b,Ha)
Model_ins = GridapMHDCalculations.models.channel_model(
                (nX,nY,nZ);
                b = b,
                L = L,
                mesh_map = map,
                cw = 0
                )
Model_cond = GridapMHDCalculations.models.channel_model(
                (nX,nY,nZ);
                b = b,
                L = L,
                mesh_map = map,
                cw = 200
                )

#Define solvers (direct solvers MUMPS)
solver_direct = Dict(
    :solver => :petsc,
    :matrix_type    => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type    => Vector{PetscScalar},
    :petsc_options  => "-snes_monitor -ksp_error_if_not_converged true \\
                        -ksp_converged_reason -ksp_type preonly -pc_type lu \\
                        -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 0 \\
                        -mat_mumps_icntl_28 1 -mat_mumps_icntl_29 2 -mat_mumps_icntl_4 3 \\
                        -mat_mumps_cntl_1 0.001 -mat_mumps_icntl_14 30 ",
    :niter          => 100,
    :rtol           => 1e-5,
    :initial_values => Dict(
      :u => U_inlet,
      :j => VectorValue(0.0,0.0,0.0),
      :p => 0.0,
      :φ => 0.0,
    ),
)

solver_blocks = Dict(
    :solver => :h1h1blocks,
    :matrix_type    => SparseMatrixCSR{Float64,Int},
    :vector_type    => Vector{Float64},
    :block_solvers  => [:petsc_mumps, :petsc_cg_jacobi, :petsc_gmres_amg],
    :petsc_options  => "-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason",
    )

#Define the FE spaces
FE_spaces = Dict(
    :order_u => 2,
    :order_j => 1,
    :fluid_disc => :Qk_dPkm1,
    :current_disc => :H1,
    )

println("Starting steady state driver")

SteadyState(;
  title = "channel_ins_direct",
  path = _path,
  backend = :mpi,
  np = _np,
  modelGen = Model_ins,
  Ha = Ha,
  N = Ha^2/Re,
  Bfield = B,
  u_inlet = U_inlet,
  source = VectorValue(0.0,0.0,0.0),
  mesh2vtk = false,
  solver = solver_direct,
  convection = :newton,
  fespaces = FE_spaces,
  post_process = GridapMHDCalculations.post_process_basic,
  order_pp = 2,
)

SteadyState(;
  title = "channel_cond_direct",
  path = _path,
  backend = :mpi,
  np = _np,
  modelGen = Model_cond,
  Ha = Ha,
  N = Ha^2/Re,
  Bfield = B,
  u_inlet = U_inlet,
  source = VectorValue(0.0,0.0,0.0),
  mesh2vtk = false,
  solver = solver_direct,
  convection = :newton,
  fespaces = FE_spaces,
  post_process = GridapMHDCalculations.post_process_basic,
  order_pp = 2,
)

SteadyState(;
  title = "channel_ins_blocks",
  path = _path,
  backend = :mpi,
  np = _np,
  modelGen = Model_ins,
  Ha = Ha,
  N = Ha^2/Re,
  Bfield = B,
  u_inlet = U_inlet,
  source = VectorValue(0.0,0.0,0.0),
  mesh2vtk = false,
  solver = solver_blocks,
  convection = :newton,
  fespaces = FE_spaces,
  post_process = GridapMHDCalculations.post_process_basic,
  order_pp = 2,
)

SteadyState(;
  title = "channel_cond_blocks",
  path = _path,
  backend = :mpi,
  np = _np,
  modelGen = Model_cond,
  Ha = Ha,
  N = Ha^2/Re,
  Bfield = B,
  u_inlet = U_inlet,
  source = VectorValue(0.0,0.0,0.0),
  mesh2vtk = false,
  solver = solver_blocks,
  convection = :newton,
  fespaces = FE_spaces,
  post_process = GridapMHDCalculations.post_process_basic,
  order_pp = 2,
)



