function Run_test_channel(np::NTuple{3,Integer};  #Number of processes per multigrid level
  backend::Union{Nothing,Symbol} = :sequential, #Backend for the mpi calculation
  b::Real = 1.0,         #Channel aspect ratio
  L::Real = 4.0,         #Channel lenght ratio
  Ha::Real = 10,         #Hartmann number
  Re::Real = 1,          #Reynolds number
  nX::Integer = 6,       #Mesh cells X direction
  nY::Integer = 6,       #Mesh cells Y direction
  nZ::Integer = 12,      #Mesh cells Z direction
  solve::Bool = true     #Solve the channel
)

 #Define geometry and mesh
 geo = channel_geom(b,L)
 mesh = channel_mesh((nX,nY,nZ),map_Roberts(b,Ha))

 #Define the boundary fields
 U_inlet((x,y,z))=VectorValue(0.0,0.0,u_parabolic(b)(x,y))

 tags = BC_tags(["inlet","walls"],["inlet","outlet","walls"])
 values = BC_values([U_inlet])
 bounds = BC(tags,values)

 #Define the Gridap model

 mounted_insulated_channel = insulated_channel(geo, mesh, bounds)

 #Define the dimensionless numbets
 numbers = Dimensionless_numbers(;Ha=Ha,Re=Re)

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

 #Call the steady state driver

 out = SteadyState(mounted_insulated_channel, numbers;
         title = "channel",
         path = "./data/tests/mpi",
         backend = backend,
         np = np,
         solver = :julia,
         convection = :newton,
         fespaces = FEspaces_options(:Qk_dPkm1,:H1),
         post_process = [writeFields_vtk,gradp_check],
         solve = solve, 
         )

 println("-----------------------------")
 println("Numerial pressure gradient at the outlet:")
 println(out[2])
 println("-----------------------------")

 kp_Shercliff=kp_shercliff_cartesian(b,Ha)

 println("-----------------------------")
 println("Analitical pressure gradient:")
 println(kp_Shercliff)
 println("-----------------------------")

 return out[2]  
end
