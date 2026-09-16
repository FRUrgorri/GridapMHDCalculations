using GridapMHDCalculations

"""
  SteadyState(; <keyword arguments>)

Driver that solves an MHD inductionless problem in steady state.

# Arguments
- `backend = nothing`: backend for parallelization. Values: nothing, :sequential or :mpi.
- `np = nothing`: array describing mesh partitioning for parallelization.
- `title = "Solid"`: job title used for saved files.
- `path = "."`: path where saved files are stored.
- `modelGen = nothing, function that generates a Gridap model with the shape (parts, ranks) -> model, tags. See Meshers module for examples 
- `normalization = :mhd` normalization of the problem variables (pressure) :mhd or :cfd
- `Ha = 10.0`: Hartmann number.
- `Re = 1.0`: Reynolds number.
- `N = nothing`: interaction number.
- `convection = true`: toggle for the weak form convective term.
- `B = VectorValue(0.0,1.0,0.0)`: external magnetic field B((x,y,z))/B0.
- `U_inlet = VectorValue(0.0,0.0,1.0)`: inlet velocity field U((x,y,z))/U0
- `source = VectorValue(0.0,0.0,0.0)`: momentum source F((x,y,z))/(j0·B0)
- `solve = true`: toggle to run the solver.
- `solver = :julia`: solver to be used and additional solver parameters.
- `verbose = true`: print time statistics.
- `fespaces = Dict{Symbol,Any}(:order_u => 2, :order_j => 2, :fluid_disc => :Qk_dPkm1, :current_disc => :RT): Dictionay with the fe spaces options
- `μ_BC::Real = 2.0`: Penalty factor of the no slip BC in the HdivH1 and HdivHdiv formulations
- `ζ::Real = 0.0`: Augmented Lagrangian factor for the iterative solver
"""
function SteadyState(args...;
  backend::Union{Nothing,Symbol} = nothing,
  np::Union{Nothing,Integer,NTuple{3,Integer}} = nothing,
  title::String = "MHD_SS",
#  nruns = 1,  #There is no interest a priory to repeate the same computation more than once (other than study scalability)
  path::String = ".",
  kwargs...
)

#  for ir in 1:nruns
#    _title = title*"_r$ir"
    if isa(backend,Nothing)
      @assert isa(np,Nothing)
      info, t, pp_out = _SteadyState(args...;title=title,path=path,kwargs...)
    else
      @assert backend ∈ [:sequential,:mpi]
      @assert !isa(np,Nothing)
      if backend === :sequential
        info, t, pp_out = with_debug() do distribute
          _SteadyState(args...;distribute=distribute,rank_partition=np,title=title,path=path,kwargs...)
        end
      else
        info, t, pp_out = with_mpi() do distribute
          _SteadyState(args...;distribute=distribute,rank_partition=np,title=title,path=path,kwargs...)
        end
      end
 #   end

    info[:np] = np
    info[:backend] = backend
    info[:title] = title
    map_main(t.data) do data
      for (k,v) in data
        info[Symbol("time_$k")] = v.max
      end
      save(joinpath(path,"$title.bson"),info)
    end
  end

  return pp_out
end

function _SteadyState(mounted_model::mounted_models, numbers::Dimensionless_numbers;
  title::String = "MHD_SS",
  path::String = ".",
  distribute::Union{Nothing,AbstractVector} = nothing,
  rank_partition::Union{Nothing,Integer,NTuple{3,Integer}} = nothing,           
  normalization::Symbol = :mhd,                 
  convection::Symbol = :newton,   
  solve::Bool = true,
  solver::Union{Dict,Symbol} = :julia,
  verbose::Bool = true,
  mesh2vtk::Bool = false,
#  μ = 0.0,
  μ_BC::Real = 2.0,
  ζ::Real = 0.0,
  fespaces::Dict{Symbol, Any} = Dict{Symbol,Any}(:order_u => 2, :order_j => 2, :fluid_disc => :Qk_dPkm1, :current_disc => :RT),
  post_process::Union{Nothing,Function,Vector{Function}} = nothing,
)

  info = Dict{Symbol,Any}()

  params = Dict{Symbol,Any}(
    :solve=>solve,
#    :res_assemble=>res_assemble,
#    :jac_assemble=>jac_assemble,
  )

  # Communicator
  if isa(distribute,Nothing)
    @assert isa(rank_partition,Nothing)   "Attempt to perform a serial calculations but partition tupple is given"
    rank_partition = Tuple(fill(1,3))     #Always 3D problems (even FD are computationally 3D)
    distribute = DebugArray
  end
  
  # With the mpi backend, distribute(x;kwargs) is an anonymus function definded inside with_mpi(f;kwargs) (inside PartitionedArrays). 
  # It is defined as x->distribute_with_mpi(x;kwargs) which returns an MPIArray. x is the collection that is distributed, in this case a LinearIndices array (AbstracArray) 
  parts = distribute(LinearIndices((prod(rank_partition),)))
  
  # Timer
  t = PTimer(parts,verbose=verbose)
  params[:ptimer] = t
  tic!(t,barrier=true)

  # Solver
  if isa(solver,Symbol)
    solver = default_solver_params(Val(solver))
  end
  params[:solver] = solver
  
  #Fespaces parameters
  params[:fespaces] = fespaces

  #Unpack inputs
  (;Ha,Re,N) = numbers

  # Reduced quantities
  @assert normalization ∈ [:mhd,:cfd]
  if normalization == :mhd
    α = 1.0/N
    β = 1.0/Ha^2
    γ = 1.0
  else
    α = 1.0
    β = 1.0/Re
    γ = N
  end
  
  #Build the model from mounted model structure
  
  # model = rank_partition == (1,1,1) ? mounted_model() : mounted_model(parts,rank_partition)
  # There is something prevetin the serial call of CartesianDiscreteModel inside GridapMHD.main.jl. It is probably a minnor fix but in the meantime mounted_model() is not callable

  model =  mounted_model(parts,rank_partition)  

 # params[:multigrid] = multigrid  #Ignored in single grid case

  params[:model] = model
  Ω = Interior(model)

  if mesh2vtk
    meshpath = joinpath(path, title*"_mesh")
    mkpath(meshpath)
    writevtk(model, meshpath)
  end
 
  #Prepare the input dictionary
  params[:fluid] = Dict{Symbol, Any}(
    :domain=>nothing, #For the moment, only fluid
    :α=>α,
    :β=>β,
    :γ=>γ,
    :f=>mounted_model.source,
    :B=>mounted_model.B,
    :ζᵤ => ζ,
    :convection=>convection,
    :μ => μ_BC,
  )
   
 """ 
  if (tw_s > 0.0) || (tw_Ha > 0.0)
    params[:fluid][:domain] = "fluid"
    σ_Ω = σ_field(model, Ω, cw_Ha, cw_s, tw_Ha, tw_s)
    params[:solid] = Dict(:domain=>"solid", :σ=>σ_Ω)
  end
"""

  #Unpack BC tags and values
  
  (;tags_U, tags_J, tags_φ) = mounted_model.BCs.tags
  (;U, J, φ) = mounted_model.BCs.values

  # Boundary conditions dictionary
    j_BC = Dict(
      :tags => tags_J,
      :values => J
    )

    u_BC = Dict(
      :tags => tags_U,
      :values => U
    )
  
  params[:bcs] = Dict(
    :u => u_BC,
    :j => j_BC,
#    :thin_wall=>thinWall_params, #TBD
  )

  if fespaces[:current_disc] == :H1
    if isempty(tags_φ)
      params[:bcs][:φ] = Dict(:tags=>[])
      params[:fespaces][:φ_constrain] = :zeromean
    else
      params[:bcs][:φ] = Dict(
        :tags => tags_φ,
        :values => φ
        )
    end
  end

"""
TBD: Allow a more general stabilization (at least a bit)
  # Stabilization method
  if μ > 0
    ĥ = b/nl[1]    See how the cell size is computed in GridapTritium
    params[:bcs][:stabilization] = Dict(:μ=>μ*ĥ, :domain=>"fluid")
  end
"""
  toc!(t,"pre_process")

  # Solve it
  if !uses_petsc(params[:solver])
    xh,fullparams,info = main(params;output=info)
  else
    petsc_options = params[:solver][:petsc_options]
    xh,fullparams,info = GridapPETSc.with(args=split(petsc_options)) do
      xh,fullparams,info = main(params;output=info)
      GridapPETSc.gridap_petsc_gc() # Destroy all PETSc objects
      return xh,fullparams,info
    end
  end
  
  

  #Post process functions
  if isnothing(post_process)
    println("No postprocess actions")
  else
    #Construct the output_info and execute the selected postprocess function
    outputs=output_info(xh,Ω,mounted_model.B,path,title)
    pp_out=exec_post_process(post_process)(outputs)
    toc!(t,"post_process")
  end

  t = fullparams[:ptimer]

  if verbose
    display(t)
  end


# Info about the solution
  info[:model] = mounted_model
  info[:order_u] = fespaces[:order_u]
  info[:order_j] = fespaces[:order_j]
  info[:ncells] = num_cells(model)
  info[:fluid_disc] = fespaces[:fluid_disc]
  info[:current_disc] = fespaces[:current_disc]
  info[:ndofs] = length(get_free_dof_values(xh))
  info[:ndofs_u] = length(get_free_dof_values(xh[1]))
  info[:ndofs_p] = length(get_free_dof_values(xh[2]))
  if fespaces[:current_disc] == :RT
    info[:ndofs_j] = length(get_free_dof_values(xh[3]))
    info[:ndofs_φ] = length(get_free_dof_values(xh[4]))
  else
    info[:ndofs_j] = "No current dofs"
    info[:ndofs_φ] = length(get_free_dof_values(xh[3]))
  end
  info[:normalization] = normalization
  info[:Re] = Re
  info[:Ha] = Ha
  info[:N] = N
  info[:convection] = convection
  info[:ζ] = ζ
  info[:μ_BC] = μ_BC
  
#  info[:cw] = cw
  
#  info[:μ] = μ

  return info, t, pp_out 
end

