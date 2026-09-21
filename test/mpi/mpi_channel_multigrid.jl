function Run_test_multigrid(np::NTuple{3,Integer};
    backend::Union{Nothing,Symbol} = :sequential,  #Number of processes per multigrid level
    b::Real = 1.0,            #Channel aspect ratio
    L::Real = 4.0,            #Channel lenght ratio
    Ha::Real = 10,             #Hartmann number
    Re::Real = 1,             #Reynolds number
    nX::Integer = 12,          #Mesh cells X direction
    nY::Integer = 12,          #Mesh cells Y direction
    nZ::Integer = 24,          #Mesh cells Z direction
    map::Function = identity, #Mesh map function
    nrefs::Integer = 2,       #Refinement factor
    levels::Integer = 2,      #Refinement levels
    ζ::Real = 10.0,           #Augmented Lagrangian
    μ_BC::Real = 2.0,         #Penalty parameter for the no_slip BC in the HdivH1 and HdivHdiv formulation 
    solve::Bool = true,       #Solve the test
    )

    #Define geometry and mesh
    geo = channel_geom(b,L)
     
    nX_c, nY_c, nZ_c = round.(Int,(nX,nY,nZ)./(nrefs*(levels-1)))  #Coarser multigrid level
    mesh = channel_mesh((nX_c,nY_c,nZ_c),map,levels,nrefs)

    #Define the boundary fields
    U_inlet((x,y,z))=VectorValue(0.0,0.0,u_parabolic(b)(x,y))

    tags = BC_tags(["inlet","walls"],["inlet","outlet","walls"])
    values = BC_values([U_inlet])
    bounds = BC(tags,values)

    #Define the Gridap model 
    mounted_insulated_channel = insulated_channel(geo, mesh, bounds)
    
    #Define the dimensionless numbets
    numbers = Dimensionless_numbers(;Ha=Ha,Re=Re)

    #Define the FE formulation
    FE_spaces = Dict(:order_u => 1, :fluid_disc => :RT,
                     :order_j => 0, :current_disc => :H1
                    )

    #Define multigrid solver 
    solver_multigrid = Dict(
        :solver => :h1h1blocks,
        :niter => 1,        #This are the maximum iteration of the non-linear Newton-Raphson solver
        :niter_ls => 1,     #This is the maximum iterations of external Kirilov solver loop (FGMRES)  
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

    out = SteadyState(mounted_insulated_channel, numbers;
        title = "channel_multigrid",
        path = "./results/tests/mpi",
        solver = solver_multigrid,
        backend = backend,
        np = np,
        convection = :none,
        fespaces = FE_spaces,
        solve = solve, 
        ζ = ζ,
        μ_BC = μ_BC,
        post_process = [writeFields_vtk,gradp_check,noSlip_check],
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
  
  println("Average the velocity components at the channel wall:")
  println(out[3][1])
  println(out[3][2])
  println(out[3][3])
  println("-----------------------------")

  return out[2], out[3]  
end
