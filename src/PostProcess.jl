#Post processing functions for the SteadyState driver

"""
post_process(args)

Most basic postprocess function, it only generates the cell fields and write them in paraview

#Arguments

-`xh : Solution CellFields
-`B : External magnetic field imposed
-`Ω: Model interior
-`path: path for the output vtk file
-`title: title for the output vtk file

"""

function post_process(xh, Ω, B, path, title; order=2)
  if length(xh) == 4
    cellfields = _post_process_4fields(xh, Ω, B)
  elseif length(xh) == 3
    cellfields = _post_process_3fields(xh, Ω, B)
  else
    error("post_process expects 3 or 4 fields, got $(length(xh))")
  end
  writevtk(Ω, joinpath(path, title), order=order, cellfields=cellfields)
  nothing
end

function _post_process_4fields(xh, Ω, B)
  uh, ph, jh, φh = xh[1], xh[2], xh[3], xh[4]
 
  div_jh = ∇·jh
  div_uh = ∇·uh
  grad_p = ∇·ph
  
  cellfields=[
    "uh"=>uh,
    "ph"=>ph,
    "jh"=>jh,
    "phi"=>φh,
    "div_uh"=>div_uh,
    "div_jh"=>div_jh,
    "grad_p"=>grad_p,
    "B" => CellField(B, Ω)
  ]
end

function _post_process_3fields(xh, Ω, B)
  uh, ph, φh = xh[1], xh[2], xh[3]
  div_uh = ∇·uh
  grad_p = ∇·ph
  
  cellfields=[
    "uh"=>uh,
    "ph"=>ph,
    "phi"=>φh,
    "div_uh"=>div_uh,
    "grad_p"=>grad_p,
    "B" => CellField(B, Ω)
  ]
end
