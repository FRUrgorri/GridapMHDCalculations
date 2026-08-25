using GridapMHDCalculations
using Gridap
using SparseArrays

#Inputs of the calculation

b = 1.0         #Channel aspect ratio
L = 4.0         #Channel lenght ratio
Ha = 10         #Hartmann number
Re = 1          #Reynolds number

#Mesh inputs (Coarser level)

nX = 2
nY = 2
nZ = 4


#Define the boundary fields
U_inlet((x,y,z))=VectorValue(0.0,0.0,GridapMHDCalculations.u_parabolic(b)(x,y))
B((x,y,z))=VectorValue(0.0,1.0,0.0)

#Define the Gridap model 

#map = GridapMHDCalculations.models.map_Roberts(b,Ha)
Model = GridapMHDCalculations.models.channel_model(
                (nX,nY,nZ),        # Number of cells 
                2;                 # Number of multigrid levels
                nrefs = 2,         # Refinement factor per level
                b = b,
                L = L,
#                mesh_map = map
                )

#Define the FE formulation
FE_spaces = Dict(:order_u => 1, :fluid_disc => :RT,
                 :order_j => 1, :current_disc => :H1
                 )

#Define multigrid solver 
solver_multigrid = Dict(
    :solver => :h1h1blocks,
    :niter => 1,
    :matrix_type    => SparseMatrixCSC{Float64,Int},
    :vector_type    => Vector{Float64},
    :block_solvers  => [:gmg, :petsc_cg_jacobi, :petsc_gmres_amg],
    :petsc_options  => "-ksp_monitor -ksp_error_if_not_converged true -ksp_converged_reason",
    :initial_values => Dict(
      :u => U_inlet,
      :j => VectorValue(0.0,0.0,0.0),
      :p => 0.0,
      :φ => 0.0,
    ),
)


#Call the steady state driver

SteadyState(;
  title = "channel_multigrid_test",
  path = "./results_test",
#  backend = :sequential,
#  np = (2,2,2),
  modelGen = Model,
  Ha = Ha,
  N = Ha^2/Re,
  Bfield = B,
  u_inlet = U_inlet,
  source = VectorValue(0.0,0.0,0.0),
  mesh2vtk = false,
  solver = solver_multigrid,
  convection = :none,
  fespaces = FE_spaces,
  post_process = GridapMHDCalculations.post_process_basic,
  order_pp = max(FE_spaces[:order_u],FE_spaces[:order_j]),
#  solve = false,
)