using DrWatson
@quickactivate "GridapMHDCalculations" #This macro activates the project meaning that it is not necessary to do add --project. I guess it increasses reproducibility...

using GridapMHDCalculations
using GridapMHDCalculations: u_parabolic, outlet_U, outlet_J, inlet_p, wall_φ, writeFields_vtk 
using GridapMHDCalculations.models
using GridapMHDCalculations.models: insulated_channel, channel_geom, channel_mesh
using Gridap

###########Inputs########

Nxy = [4,8,12,16]
Nz = 10
solver = :julia 
Ha = 10
Re = 1
b = 1
L = 4

#Build the dictionaries
params = @dict Nxy Nz Ha Re b L solver
params_list = dict_list(params)

#Define the functions for running the analysis

function Run_channel(dict::Dict{Symbol,Any},path,title)

  #Unpack from the input dict
  @unpack Nxy, Nz, Ha, Re, b, L, solver = dict

  #Build geometry and mesh
  geo = channel_geom(b,L)
  mesh = channel_mesh((Nxy,Nxy,Nz),map_Roberts(b,Ha))

  #Define the boundary fields
  U_inlet((x,y,z))=VectorValue(0.0,0.0,u_parabolic(b)(x,y))

  tags = BC_tags(["inlet","walls"],["inlet","outlet","walls"])
  values = BC_values([U_inlet])
  bounds = BC(tags,values)

  #Define the mounted Gridap model

  mounted_insulated_channel = insulated_channel(geo, mesh, bounds)

  #Define the dimensionless numbets
  numbers = Dimensionless_numbers(;Ha=Ha,Re=Re)

  #Define the FE spaces parameters
  FE_spaces=Dict(:order_u => 2, :order_j => 1, :fluid_disc => :Qk_dPkm1, :current_disc => :H1)

  #Make the simulations
  monitors = SteadyState(mounted_insulated_channel, numbers;
          title = title,
          path = path,
          solver = solver,
          convection = :newton,
          fespaces = FE_spaces,
          post_process = [outlet_U, outlet_J, inlet_p, wall_φ, writeFields_vtk], 
          )

    #Build the outputs
    out_dict=Dict{Symbol,Any}(:Nxy=>dict[:Nxy],:Nz => dict[:Nz])   
    out_dict[:Uz_out] = monitors[1][3]
    out_dict[:Jx_out] = monitors[2][1]
    out_dict[:Jy_out] = monitors[2][2]
    out_dict[:p_in] = monitors[3]
    out_dict[:φ_wall] = monitors[4] 


   return out_dict
end     

function Run_mesh_analysis(list::Vector{Dict{Symbol, Any}})
    for (i,d) in enumerate(list)

        file = savename("Mesh",d,"bson";accesses=(:Nxy,:Nz))
        title = savename("Mesh",d ;accesses=(:Nxy,:Nz))

        subfolder = savename(d;accesses=(:Ha,:Re))
        folder = savename("ins_channel",d;accesses=(:L,:b))

        file_path = datadir("mesh_sensitivity",joinpath(folder,subfolder,file))
        vtk = datadir("mesh_sensitivity", joinpath(folder,subfolder),"vtk")

        out_dict = Run_channel(d,vtk,title)
        @tagsave(file_path,out_dict)
    end
end

#Run the analysis
Run_mesh_analysis(params_list)   #TBD: Make the analysis adaptative so the simulation is run only if necessary
