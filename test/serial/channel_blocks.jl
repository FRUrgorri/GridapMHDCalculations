function Run_test_channel_blocks(;
    b::Real = 1.0,         #Channel aspect ratio
    L::Real = 4.0,         #Channel lenght ratio
    Ha::Real = 10,         #Hartmann number
    Re::Real = 1,          #Reynolds number
    nX::Integer = 6,       #Mesh cells X direction
    nY::Integer = 6,       #Mesh cells Y direction
    nZ::Integer = 12,      #Mesh cells Z direction
    solve::Bool = true     #Solve the problem
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
  
 
 solver_blocks = Dict(
        :solver => :h1h1blocks,
        :niter => 1,       #This are the maximum iteration of the non-linear Newton-Raphson solver
        :niter_ls => 1,     #This is the maximum iterations of external Kirilov solver loop (FGMRES)
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
  FE_spaces=Dict(:order_u => 1, :order_j => 1, :fluid_disc => :RT, :current_disc => :H1)

  #Call the steady state driver

  out = SteadyState(mounted_insulated_channel, numbers;
          title = "channel",
          path = "./data/tests/serial",
          solver = :julia,
          convection = :newton,
          fespaces = FE_spaces,
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