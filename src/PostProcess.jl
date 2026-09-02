#File for post process functions

"""
exec_post_process(pp_function;kargs...)

Post process function selector (see SteadyState.jl for the execution example)

#Arguments

-`pp_function: post process function to be executed. It has the positional arguments (xh, Ω, B, path, title) and any key word arguments
-`kargs: key word arguments of pp_function 

"""
function exec_post_process(pp_function::Function)
     (output::output_info) -> pp_function(output)

end

"""
post_process_basic(args)

Most basic postprocess function, it generates the cell fields and write them in a vtk file for paraview. It optionally passes an the fields defined in the out_field

#Arguments

-`output`: output information structure
-`path: path for the output vtk file
-`title: title for the output vtk file

#keyword arguments
-`order_pp`: Interpolation order of the vtk file
-`pass_fields`: Tuple with the names (strings) of the fields that want to be passed by the function

"""

function post_process_basic(output::output_info; pass_fields::Union{Nothing,Tuple} = nothing) 
                              
  if length(output.xh) == 4
    _cellfields = pp_4fields(output)
  elseif length(output.xh) == 3
    _cellfields = pp_3fields(output)
  else
    error("post_process expects 3 or 4 fields, got $(length(xh))")
  end
  writevtk(output.Ω, joinpath(output.path, output.title), order=output.order, cellfields=_cellfields)
  
  #Check if an output field is necessary for further pos-process
  if !isnothing(pass_fields)
    field_dict = Dict(_cellfields)
    cellfields_pass=map(x->field_dict[x],pass_fields)
  else
    cellfields_pass = nothing
  end

  return cellfields_pass
end

function pp_4fields(output::output_info)
  uh, ph, jh, φh = output.xh

  div_jh = ∇·jh
  div_uh = ∇·uh
  grad_p = ∇·ph
  grad_phi = ∇·φh

  cellfields=[
    "uh"=>uh,
    "ph"=>ph,
    "jh"=>jh,
    "phi"=>φh,
    "div_uh"=>div_uh,
    "div_jh"=>div_jh,
    "grad_p"=>grad_p,
    "grad_phi" =>grad_phi,
    "B" => CellField(output.B, output.Ω)
  ]

  return cellfields
end

function pp_3fields(output::output_info)
  uh, ph, φh = output.xh
  
  div_uh = ∇·uh
  grad_p = ∇·ph
  grad_phi = ∇·φh

  jh = uh×output.B - grad_phi 
 #  div_jh = ∇·jh
  
  cellfields=[
    "uh"=>uh,
    "ph"=>ph,
    "phi"=>φh,
    "div_uh"=>div_uh,
    "grad_p"=>grad_p,
    "grad_phi"=>grad_phi,
    "B" => CellField(output.B, output.Ω),
    "jh" => jh,
   #    "div_jh"=>div_jh,  # TBS: Unknown error when writting 
  ]

  return cellfields
end

function pp_gradp_check(output::output_info)
  
  grad_p,=post_process_basic(output; pass_fields=("grad_p",))

  #Computing the pressure axial gradient at the outlet to compare it with an analytical formula
  model=get_model(output.Ω)
  Γ_out = Boundary(model;tags="outlet")
  dΓ= Measure(Γ_out,output.order+1)
  
  return sum(∫(-grad_p)*dΓ)[3]/sum(∫(1.0)*dΓ)
end

function pp_Noslip_check(output::output_info)

  uh, grad_p = post_process_basic(output; pass_fields=("uh","grad_p"))
  
  #Computing the pressure axial gradient at the outlet to compare it with an analytical formula
  model=get_model(output.Ω)

  Γ_out = Boundary(model;tags="outlet")
  dΓ_out= Measure(Γ_out,output.order+1)
  kp= sum(∫(-grad_p)*dΓ_out)[3]/sum(∫(1.0)*dΓ_out)

  #Check the no_slip_condition (using average velocity)
  
  #Collect wall tags
  wall_tags = []
  tag_names=get_tag_names(model)
  map(tag_names) do tag
    if tag in ("insulated","conducting","thin_wall")
      push!(wall_tags,tag)
    end
  end

  Γ_wall = Boundary(model;tags="insulated")
  dΓ_wall= Measure(Γ_wall,output.order+1)
  u_wall = sum(∫(uh)*dΓ_wall)/sum(∫(1.0)*dΓ_wall)


  return kp, u_wall

end

#Utilities
get_model(Ω) = Ω.model

get_tag_names(model::DiscreteModel)=get_face_labeling(model).tag_to_name  #Equivalent to get_face_labeling(model) from Gridap
get_tag_names(model::GridapDistributed.DistributedDiscreteModel) = get_face_labeling(model).labels  #Take the names from the first part
