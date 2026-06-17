using GridapMHD

using SparseMatricesCSR
using GridapPETSc
using Gridap

#Inputs of the calculation

b = 1.0         #Channel aspect ratio
L = 4.0         #Channel lenght ratio
Ha = 10         #Hartmann number
Re = 1          #Reynolds number

#Mesh inputs

nX = 6
nY = 6
nZ = 12


#Define the boundary fields
U_inlet((x,y,z))=VectorValue(0.0,0.0,GridapMHD.u_parabolic(b)(x,y))
B((x,y,z))=VectorValue(0.0,1.0,0.0)

#Define the Gridap model 

map = GridapMHD.Meshers.map_Roberts(b,Ha)
Model = GridapMHD.Meshers.channel_model(
                (nX,nY,nZ);
                b = b,
                L = L,
                mesh_map = map
                )

#Define solver (direct solver MUMPS in H1Hdiv formulation)

solver = Dict(
    :solver => :petsc,
    :matrix_type    => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type    => Vector{PetscScalar},
    :petsc_options  => "-snes_monitor -ksp_error_if_not_converged true 
                        -ksp_converged_reason -ksp_type preonly -pc_type lu 
                        -pc_factor_mat_solver_type mumps -mat_mumps_icntl_28 1 
                        -mat_mumps_icntl_29 2 -mat_mumps_icntl_4 3 
                        -mat_mumps_cntl_1 0.001",
    :niter          => 100,
    :rtol           => 1e-5,
    :initial_values => Dict(
      :u => U_inlet,
      :j => VectorValue(0.0,0.0,0.0),
      :p => 0.0,
      :φ => 0.0,
    ),
)

#Call the steady state driver

xh,Ω = SteadyState(;
  title = "channel_test",
  path = "./results",
  backend = :mpi,
  np = (2, 2, 1),
  modelGen = Model,
  Ha = Ha,
  N = Ha^2/Re,
  Bfield = B,
  u_inlet = U_inlet,
  source = VectorValue(0.0,0.0,0.0),
  mesh2vtk = false,
  solver = solver,
  convection = :newton,
#  solve = false,
)

GridapMHD.post_process(xh, Ω, B, "./results", "channel_test")
