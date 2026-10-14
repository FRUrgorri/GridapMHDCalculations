using DrWatson
@quickactivate "GridapMHDCalculations" #This macro activates the project meaning that it is not necessary to do add --project. I guess it increasses reproducibility...

using GridapMHDCalculations
using GridapMHDCalculations: writeFields_vtk
using GridapMHDCalculations.models
using GridapMHDCalculations.models: insulated_channel, channel_geom, channel_mesh
using Gridap
using SparseArrays, SparseMatricesCSR

###########Inputs########

#"Fixed" parameters
Ha = 10
Re = 1
b = 1
L = 4
ζ = 10                                         #Augmented Lagrangian

# Distribution of the parallel cases
@assert length(ARGS) == 1 "Please provide the ranks per multigrid level"
_np = Tuple{Int,Int,Int}(parse.(Int, split(ARGS[1], ",")))       #Ranks per multigrid level

# Coarse mesh
Nxy_c = 8       #Mesh in the coarser level (always the same)
Nz_c =  4

#Multigrid parameters
μ_BC = [25,50,75,100]                            #Penalty parameter for the no_slip BC in the HdivH1 and HdivHdiv formulation 
mapFunction = map_Roberts(b,Ha)                #Mesh map function
mg_levels = 4                                  #Multigrid levels
nrefs = 2                                      #Refinement level

#Build the dictionaries
params = @dict Nxy_c Nz_c Ha Re b L ζ μ_BC mapFunction mg_levels nrefs
params_list = dict_list(params)

function Run_mg_channel(dict::Dict{Symbol,Any},path::String,title::String,np::NTuple{3,Integer})

  #Unpack from the input dict
  @unpack Nxy_c, Nz_c, Ha, Re, b, L, ζ, μ_BC, mapFunction, mg_levels, nrefs = dict

  println("Running $title with parameters: 
        Nxy_c=$(Nxy_c) | Nz_c=$(Nz_c) | mg_levels=$(mg_levels) | nrefs=$(nrefs) 
        Ha=$(Ha) | Re=$(Re) | b=$(b) | L=$(L) 
        ζ=$(ζ) | μ_BC=$(μ_BC) | mapFunction=$(mapFunction)"
    )

  #Build geometry and mesh
  geo = channel_geom(b,L)
  mesh = channel_mesh((Nxy_c,Nxy_c,Nz_c),mapFunction,mg_levels,nrefs)

  #Define the boundary fields
  U_inlet((x,y,z))=VectorValue(0.0,0.0,u_parabolic(b)(x,y))

  tags = BC_tags(["inlet","walls"],["inlet","outlet","walls"])
  values = BC_values([U_inlet])
  bounds = BC(tags,values)

  #Define the mounted Gridap model

  mounted_insulated_channel = insulated_channel(geo, mesh, bounds)

  #Define the dimensionless numbets
  numbers = Dimensionless_numbers(;Ha=Ha,Re=Re)

   

  #Define multigrid solver 
  solver_mg = Dict(
        :solver => :h1h1blocks,
        :niter => 6,        #This are the maximum iteration of the non-linear Newton-Raphson solver
        :niter_ls => 10,    #This is the maximum iterations of external Kirilov solver loop (FGMRES)  
        :matrix_type    => SparseMatrixCSC{Float64,Int},
        :vector_type    => Vector{Float64},
        :block_solvers  => [:gmg, :petsc_cg_jacobi, :petsc_gmres_amg],
        :petsc_options  => "-ksp_monitor -ksp_error_if_not_converged false -ksp_converged_reason",
        :initial_values => Dict(
            :u => U_inlet,
            :j => VectorValue(0.0,0.0,0.0),
            :p => 0.0,
            :φ => 0.0,
        ),
    )

  #Make the simulations
   SteadyState(mounted_insulated_channel, numbers;
          title = title,
          path = path,
          backend = :mpi,
          np = np,
          solver = solver_mg,
          fespaces = FEspaces_options(:RT,:H1),
          ζ = ζ,
          μ_BC = μ_BC,
          post_process = [writeFields_vtk], 
    )
  

   return nothing
end     

function Run_mg_analysis(list::Vector{Dict{Symbol, Any}},np)
    
    subfolder = savename(list[1];accesses=(:Ha,:Re))
    folder = savename("ins_channel",list[1];accesses=(:L,:b))

    for d in list
        map_str="mapFunction=Roberts"
            
        title = savename(map_str, d ;ignores=("Ha","Re","b","L"))
        vtk_path = datadir("channel_mg_analysis", joinpath(folder,subfolder,"vtk"))

        Run_mg_channel(d,vtk_path,title,np)        
        
    end
end

#Run the analysis

Run_mg_analysis(params_list,_np)   