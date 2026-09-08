function Run_test_channel_blocks(;
    b::Real = 1.0,         #Channel aspect ratio
    L::Real = 4.0,         #Channel lenght ratio
    Ha::Real = 10,         #Hartmann number
    Re::Real = 1,          #Reynolds number
    nX::Integer = 6,       #Mesh cells X direction
    nY::Integer = 6,       #Mesh cells Y direction
    nZ::Integer = 12,      #Mesh cells Z direction
  )

  #Define the boundary fields
  U_inlet((x,y,z))=VectorValue(0.0,0.0,u_parabolic(b)(x,y))
  B((x,y,z))=VectorValue(0.0,1.0,0.0)

  #Define the Gridap model 

  Model = channel_model(
                (nX,nY,nZ);
                b = b,
                L = L,
                mesh_map = map_Roberts(b,Ha)
                )

  #Define solver (direct solver MUMPS in H1Hdiv formulation)
  
 
 solver_blocks = Dict(
        :solver => :h1h1blocks,
        :niter => 10,        #This I think it is the maximum iteration of the external Kirilov solver (FGMRES)
        :niter_ls => 3,     #This I think it is the maximum iteration in each step of the loop. Over this loop there is a Kirilov solver (FGMRES)
        :matrix_type    => SparseMatrixCSC{Float64,Int},
        :vector_type    => Vector{Float64},
        :block_solvers => [:julia,:julia,:julia],
        :petsc_options  => "-ksp_monitor -ksp_error_if_not_converged true -ksp_converged_reason",
        :initial_values => Dict(
            :u => U_inlet,
            :j => VectorValue(0.0,0.0,0.0),
            :p => 0.0,
            :φ => 0.0,
            ),
        )
  
  #Define the FE spaces parameters
  FE_spaces=Dict(:order_u => 1, :order_j => 2, :fluid_disc => :RT, :current_disc => :H1)

  #Call the steady state driver

  kp = SteadyState(;
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
    solver = solver_blocks,
    convection = :newton,
    fespaces = FE_spaces,
    post_process = pp_gradp_check,
   # solve = false,
    
  )

  println("-----------------------------")
  println("Numerial pressure gradient at the outlet:")
  println(kp)
  println("-----------------------------")

  kp_Shercliff=kp_shercliff_cartesian(b,Ha)

  println("-----------------------------")
  println("Analitical pressure gradient:")
  println(kp_Shercliff)
  println("-----------------------------")

  return kp  
end