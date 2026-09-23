using DrWatson
@quickactivate "GridapMHDCalculations" #This macro activates the project meaning that it is not necessary to do add --project. I guess it increasses reproducibility...

using GridapMHDCalculations
using GridapMHDCalculations: u_parabolic, outlet_U, outlet_J, inlet_p, wall_φ, noSlip_check
using GridapMHDCalculations.models
using GridapMHDCalculations.models: insulated_channel, channel_geom, channel_mesh
using Gridap
using SparseArrays, SparseMatricesCSR

###########Inputs########

Nxy = [8,16,24,32,40]               #Crossectional cells
Nz = [8,16]                        #Axial cells
 
Ha = 10
Re = 1
b = 1
L = 4

ζ = 10                                         #Augmented Lagrangian
μ_BC = [2, 6, 10, 25, 50, 100]                 #Penalty parameter for the no_slip BC in the HdivH1 and HdivHdiv formulation 
map_function = [identity,map_Roberts(b,Ha)]    #Mesh map function
mg_levels = [2,3,4]                            #Multigrid levels
nrefs = 2                                      #Refinement level

#Build the dictionaries
params = @dict Nxy Nz Ha Re b L ζ μ_BC map_function mg_levels nrefs
params_list = dict_list(params)

#Define the functions for running the analysis

function Run_mg_channel(dict::Dict{Symbol,Any},path::String,title::String,np::NTuple{3,Integer})

  #Unpack from the input dict
  @unpack Nxy, Nz, Ha, Re, b, L, ζ, μ_BC, map_function, mg_levels, nrefs = dict

  #Build geometry and mesh
  geo = channel_geom(b,L)
  Nxy_c, Nz_c = round.(Int,(Nxy,Nz)./(nrefs*(mg_levels-1)))  
  mesh = channel_mesh((Nxy_c,Nxy_c,Nz_c),map_function,mg_levels,nrefs)

  #Define the boundary fields
  U_inlet((x,y,z))=VectorValue(0.0,0.0,u_parabolic(b)(x,y))

  tags = BC_tags(["inlet","walls"],["inlet","outlet","walls"])
  values = BC_values([U_inlet])
  bounds = BC(tags,values)

  #Define the mounted Gridap model

  mounted_insulated_channel = insulated_channel(geo, mesh, bounds)

  #Define the dimensionless numbets
  numbers = Dimensionless_numbers(;Ha=Ha,Re=Re)

      #Define the FE formulation
    FE_spaces = Dict(:order_u => 1, :fluid_disc => :RT,
                     :order_j => 0, :current_disc => :H1
                    )

    #Define multigrid solver 
    solver_mg = Dict(
        :solver => :h1h1blocks,
        :niter => 6,        #This are the maximum iteration of the non-linear Newton-Raphson solver
        :niter_ls => 10,    #This is the maximum iterations of external Kirilov solver loop (FGMRES)  
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

  #Make the simulations
  monitors = SteadyState(mounted_insulated_channel, numbers;
          title = title,
          path = path,
          backend = :mpi,
          np = np,
          solver = solver_mg,
          fespaces = FE_spaces,
          post_process = [outlet_U, noSlip_check, outlet_J, inlet_p, wall_φ], 
          )

    #Build the outputs
    out_dict=copy(dict) 
    out_dict[:Uz_out] = monitors[1][3]
    out_dict[:Ux_wall] = monitors[2][1]
    out_dict[:Uy_wall] = monitors[2][2]
    out_dict[:Uz_wall] = monitors[2][3]
    out_dict[:Jx_out] = monitors[3][1]
    out_dict[:Jy_out] = monitors[3][2]
    out_dict[:p_in] = monitors[4]
    out_dict[:φ_wall] = monitors[5] 


   return out_dict
end     

function Run_mg_analysis(list::Vector{Dict{Symbol, Any}},np)
    
    subfolder = savename(list[1];accesses=(:Ha,:Re))
    folder = savename("ins_channel",list[1];accesses=(:L,:b))
    dir = datadir("channel_mg_analysis", joinpath(folder,subfolder))
    done = isdir(dir) ? Set(readdir(dir)) : String[]

    for d in list

        tag = savename(d,"bson";ignores=("Ha","Re","b","L"))
        title = savename(d ;ignores=("Ha","Re","b","L"))
        
        tag_path = joinpath(dir,tag)
        vtk_path = joinpath(dir,"vtk")     

        tag in done && continue     #Skip to the next iteration if the tag has been computed previously

        try 
            out_dict = Run_mg_channel(d,vtk_path,title,np)
            @tagsave(tag_path,out_dict)
        catch e
            @warn "skipped $tag" exception=(e, catch_backtrace())
        end
    end
end

#Run the analysis
#Run_mg_analysis(params_list)   