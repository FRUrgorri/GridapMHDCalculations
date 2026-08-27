using GridapMHDCalculations
using GridapMHDCalculations.models: map_Roberts, channel_model
using Gridap

function Run_test_channel(;
    b::Real = 1.0,         #Channel aspect ratio
    L::Real = 4.0,         #Channel lenght ratio
    Ha::Real = 10,         #Hartmann number
    Re::Real = 1,          #Reynolds number
    nX::Integer = 6,       #Mesh cells X direction
    nY::Integer = 6,       #Mesh cells Y direction
    nZ::Integer = 12,       #Mesh cells Z direction
  )

  #Define the boundary fields
  U_inlet((x,y,z))=VectorValue(0.0,0.0,GridapMHDCalculations.u_parabolic(b)(x,y))
  B((x,y,z))=VectorValue(0.0,1.0,0.0)

  #Define the Gridap model 

  Model = channel_model(
                (nX,nY,nZ);
                b = b,
                L = L,
                mesh_map = map_Roberts(b,Ha)
                )

  #Define solver (direct solver MUMPS in H1Hdiv formulation)
  """
  solver_direct = Dict(
    :solver => :petsc,
    :matrix_type    => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type    => Vector{PetscScalar},
    :petsc_options  => "-snes_monitor -ksp_error_if_not_converged true \\
                        -ksp_converged_reason -ksp_type preonly -pc_type lu \\
                        -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 0 \\
                        -mat_mumps_icntl_28 1 -mat_mumps_icntl_29 2 -mat_mumps_icntl_4 3 \\
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
  """
  #Define the FE spaces parameters
  FE_spaces=Dict(:order_u => 2, :order_j => 2, :fluid_disc => :Qk_dPkm1, :current_disc => :H1)

  #Call the steady state driver

  SteadyState(;
    title = "channel_test",
    path = "./results_test",
  #  backend = :sequential,
  #  np = (2, 2, 1),
    modelGen = Model,
    Ha = Ha,
    N = Ha^2/Re,
    Bfield = B,
    u_inlet = U_inlet,
    source = VectorValue(0.0,0.0,0.0),
    mesh2vtk = false,
    solver = :julia,
    convection = :newton,
    fespaces = FE_spaces,
    post_process = GridapMHDCalculations.post_process_basic,
    order_pp = max(FE_spaces[:order_u],FE_spaces[:order_j]),
  #  solve = false,
  )
end
