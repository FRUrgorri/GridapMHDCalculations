#Post processing functions

"""
exec_post_process(args)

Post process function selector (see SteadyState.jl for the execution example)

#Arguments

-`pp_function: post process function to be executed. It has the positional arguments (xh, Ω, B, path, title) and any key word arguments
-`kargs: key word arguments of pp_function 

"""

function exec_post_process(pp_function; kargs...)
 (xh, Ω, B, path, title; kargs...)->pp_function(xh, Ω, B, path, title;kargs...)
end


"""
post_process_basic(args)

Most basic postprocess function, it generates the cell fields and write them in a vtk file for paraview

#Arguments

-`xh : Solution CellFields
-`B : External magnetic field imposed
-`Ω: Model interior
-`path: path for the output vtk file
-`title: title for the output vtk file

"""

function post_process_basic(xh, Ω, B, path, title; order_pp=2)
  if length(xh) == 4
    cellfields = _post_process_4fields(xh, Ω, B)
  elseif length(xh) == 3
    cellfields = _post_process_3fields(xh, Ω, B)
  else
    error("post_process expects 3 or 4 fields, got $(length(xh))")
  end
  writevtk(Ω, joinpath(path, title), order=order_pp, cellfields=cellfields)
  nothing
end

function _post_process_4fields(xh, Ω, B)
  uh, ph, jh, φh = xh[1], xh[2], xh[3], xh[4]
 
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
    "B" => CellField(B, Ω)
  ]
end

function _post_process_3fields(xh, Ω, B)
  uh, ph, φh = xh[1], xh[2], xh[3]
  
  div_uh = ∇·uh
  grad_p = ∇·ph
  grad_phi = ∇·φh

  jh = uh×B - grad_phi 
#  div_jh = ∇·jh
  
  cellfields=[
    "uh"=>uh,
    "ph"=>ph,
    "phi"=>φh,
    "div_uh"=>div_uh,
    "grad_p"=>grad_p,
    "grad_phi"=>grad_phi,
    "B" => CellField(B, Ω),
    "jh" => jh,
#    "div_jh"=>div_jh,  # TBS: Unknown error when writting 
  ]
end
