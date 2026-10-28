using DrWatson
@quickactivate "GridapMHDCalculations" 

using GridapMHDCalculations
using GridapMHDCalculations: writeFields_vtk, u_parabolic
using GridapMHDCalculations.models
using GridapMHDCalculations.models: insulated_channel, channel_geom, channel_mesh
using Gridap
using GridapPETSc
using SparseArrays, SparseMatricesCSR

###########Inputs########

#"Fixed" parameters
Ha = 10
Re = 1
b = 1
L = 4

# Distribution of the parallel cases
@assert length(ARGS) == 1 "Please provide the ranks per multigrid level"
_np = Tuple{Int,Int,Int}(parse.(Int, split(ARGS[1], ",")))       #Ranks per multigrid level

# Mesh
Nxy = 30       
Nz =  16
mapFunction = map_Roberts(b,Ha)        

#Build the dictionaries
params = @dict Nxy Nz Ha Re b L mapFunction
params_list = dict_list(params)

function Run_mumps_channel(dict::Dict{Symbol,Any},path::String,title::String,np::NTuple{3,Integer})

  #Unpack from the input dict
  @unpack Nxy, Nz, Ha, Re, b, L, mapFunction = dict

  println("Running $title with parameters: 
        Nxy=$(Nxy) | Nz=$(Nz) | Ha=$(Ha) | Re=$(Re) | b=$(b) | L=$(L) | mapFunction=$(mapFunction)"
    )

  #Build geometry and mesh
  geo = channel_geom(b,L)
  mesh = channel_mesh((Nxy,Nxy,Nz),mapFunction)

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

  #Make the simulations
   SteadyState(mounted_insulated_channel, numbers;
          title = title,
          path = path,
          backend = :mpi,
          np = np,
          solver = solver_direct,
          fespaces = FEspaces_options(:Qk_dPkm1,:RT),
          post_process = writeFields_vtk, 
    )

   return nothing
end     

function Run_mumps_analysis(list::Vector{Dict{Symbol, Any}},np)
    
    subfolder = savename(list[1];accesses=(:Ha,:Re))
    folder = savename("ins_channel",list[1];accesses=(:L,:b))

    for d in list
        map_str="mapFunction=Roberts"
            
        title = savename(map_str, d ;ignores=("Ha","Re","b","L"))
        vtk_path = datadir("channel_mumps_analysis", joinpath(folder,subfolder,"vtk"))

        Run_mumps_channel(d,vtk_path,title,np)        
        
    end
end

#Run the analysis

Run_mumps_analysis(params_list,_np)   
