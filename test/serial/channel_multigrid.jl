function Run_test_multigrid(;
    b::Real = 1.0,            #Channel aspect ratio
    L::Real = 4.0,            #Channel lenght ratio
    Ha::Real = 1,             #Hartmann number
    Re::Real = 1,             #Reynolds number
    nX::Integer = 6,          #Mesh cells X direction
    nY::Integer = 6,          #Mesh cells Y direction
    nZ::Integer = 8,         #Mesh cells Z direction
    nrefs::Integer = 2,       #Refinement factor
    levels::Integer = 2,      #Refinement levels
    ζ::Real = 10.0,           #Augmented Lagrangian
    μ_BC::Real = 100.0,         #Penalty parameter for the no_slip BC in the HdivH1 and HdivHdiv formulation 
    map::Function = identity  #Mesh map function
    )

    #Define the boundary fields
    U_inlet((x,y,z))=VectorValue(0.0,0.0,u_parabolic(b)(x,y))
    B((x,y,z))=VectorValue(0.0,1.0,0.0)

    #Define the Gridap model 

    #Coarser multigrid level
    nX_c, nY_c, nZ_c = round.(Int,(nX,nY,nZ)./(nrefs*(levels-1)))  

    Model = channel_model(
                (nX_c,nY_c,nZ_c),       # Number of cells in the coarse level
                levels;                 # Number of multigrid levels
                nrefs = nrefs,          # Refinement factor per level
                b = b,
                L = L,
                mesh_map = map
                )

    #Define the FE formulation
    FE_spaces = Dict(:order_u => 1, :fluid_disc => :RT,
                     :order_j => 1, :current_disc => :H1
                    )

    #Define multigrid solver 
    solver_multigrid = Dict(
        :solver => :h1h1blocks,
        :niter => 1,        #This I think it is the maximum iteration of the gmg (geometric multigrid) internal loop. Over this loop there is a Kirilov solver (FGMRES)
        :niter_ls => 2,     #This I think it is the maximum iterations of the most external Kirilov solver loop (FGMRES) (not counting NR solver if there is convection) 
        :matrix_type    => SparseMatrixCSC{Float64,Int},
        :vector_type    => Vector{Float64},
        :block_solvers  => [:gmg, :petsc_cg_jacobi, :petsc_gmres_amg],
#        :block_solvers => [:julia,:julia,:julia],
        :petsc_options  => "-ksp_monitor -ksp_error_if_not_converged true -ksp_converged_reason",
        :initial_values => Dict(
            :u => U_inlet,
            :j => VectorValue(0.0,0.0,0.0),
            :p => 0.0,
            :φ => 0.0,
            ),
        )



        #Call the steady state driver

    kp, u_wall = SteadyState(;
        title = "channel_multigrid_test",
        path = "./results_test",
    #    backend = :sequential,
    #    np = (2,2,2),
        modelGen = Model,
        Ha = Ha,
        N = Ha^2/Re,
        Bfield = B,
        u_inlet = U_inlet,
        source = VectorValue(0.0,0.0,0.0),
        ζ = ζ,
        μ_BC = μ_BC,
        mesh2vtk = false,
        solver = solver_multigrid,
        convection = :none,
        fespaces = FE_spaces,
        post_process = pp_Noslip_check,
 #       solve = false,
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

  println("Average value of the velocity components at the channel wall:")
  println(u_wall[1])
  println(u_wall[2])
  println(u_wall[3])
  println("-----------------------------")
  
   return kp, u_wall 
end
