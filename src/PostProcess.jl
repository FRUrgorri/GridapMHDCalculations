"""
Struct with the output information of every computation

#Fields
 -`xh`:Cell fields
 -`Ω`: Model triangulation
 -`B`: External magnetic field, either a function or a VectorValue
 -`path`: path to the writting folder
 -`title`: title of the output files
 -`order`: order of the output vtk file
"""
struct output_info{C,T}
  xh::C
  Ω::T
  B::Union{Function,VectorValue{3,Float64}}
  path::String
  title::String
  order::Int64
end

#Constructors of the output_info type
output_info(xh,Ω,B,path,title) = output_info(xh,Ω,B,path,title,2)

"""
exec_post_process(pp_function;kargs...)

Post process function selector (see SteadyState.jl for the execution example)

#Arguments

-`pp_function: post process function to be executed. A vector of functions can be passed for a sequential execution. The functions must take an output_info as positional argument

"""
#Single execution
exec_post_process(pp_function::Function) = (output::output_info) -> pp_function(output)

#Sequential execution
exec_post_process(pp_functions::Vector{Function}) = (output::output_info) -> map(f->exec_post_process(f)(output),pp_functions)  


"""
writeFields_vtk(output_info)

Most basic postprocess function, it generates the cell fields and write them in a vtk file for paraview. I

#Arguments

-`output`: output information structure
"""

function writeFields_vtk(output::output_info) 
                              
  cellfields = unpack_fields(output)
  writevtk(output.Ω, joinpath(output.path, output.title), order=output.order, cellfields=cellfields)
  
  return nothing
end

"""
gradp_check(output::output_info)
Compute the pressure gradient at the "outlet" tag 
"""

function gradp_check(output::output_info)

  cellfields = unpack_fields(output)
  grad_p = Dict(cellfields)["grad_p"]   #converting into a dictionary to find the value associated to the field is perhaps not efficient?

  #Computing the pressure axial gradient at the outlet to compare it with an analytical formula
  model=get_model(output.Ω)
  Γ_out = Boundary(model;tags="outlet")
  dΓ= Measure(Γ_out,output.order+1)
  
  return sum(∫(-grad_p)*dΓ)[3]/sum(∫(1.0)*dΓ)
end

"""
noSlip_check(output::output_info)
Check the noslip BC by computing the velocity in the tags  ("insulated","conducting","thin_wall","wall","walls")

"""

function noSlip_check(output::output_info)

  cellfields = unpack_fields(output)
  uh = Dict(cellfields)["uh"] 
   
  #Collect wall tags
  model=get_model(output.Ω)
  wall_tags = String[]
  tag_names=get_tag_names(model)
  map(tag_names) do tag
    if tag ∈ ("insulated","conducting","thin_wall","wall","walls")
      push!(wall_tags,tag)
    end
  end

  Γ_wall = Boundary(model;tags=wall_tags)
  dΓ_wall= Measure(Γ_wall,output.order+1)
  
  return sum(∫(uh)*dΓ_wall)/sum(∫(1.0)*dΓ_wall)

end


#################Utilities#########################

get_model(Ω) = Ω.model

get_tag_names(model::DiscreteModel)=get_face_labeling(model).tag_to_name  
get_tag_names(model::GridapDistributed.DistributedDiscreteModel) = get_tag_names(local_views(model).items[1])

function unpack_fields(output::output_info) 
  if length(output.xh) == 4
    cellfields = unpack_4fields(output)
  elseif length(output.xh) == 3
    cellfields = unpack_3fields(output)
  else
    error("post_process expects 3 or 4 fields, got $(length(xh))")
  end

  return cellfields
end

function unpack_4fields(output::output_info)
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

function unpack_3fields(output::output_info)
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