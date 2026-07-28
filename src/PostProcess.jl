#Post processing functions

"""
post_process(args)

Most basic postprocess function, it generates the cell fields and write them in a vtk file for paraview

#Arguments

-`xh : Solution CellFields
-`B : External magnetic field imposed
-`Ω: Model interior
-`path: path for the output vtk file
-`title: title for the output vtk file

"""

function post_process_basic(xh, Ω, B, path, title; order=2)
  if length(xh) == 4
    uh, ph, jh, φh = xh[1], xh[2], xh[3], xh[4]
  elseif length(xh) == 3
    uh, ph, φh = xh[1], xh[2], xh[3]
    jh =  ∇·φh + uh×B  
  end

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

  writevtk(Ω, joinpath(path, title), order=order, cellfields=cellfields)
  nothing
end

