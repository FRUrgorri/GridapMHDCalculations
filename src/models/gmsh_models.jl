"""
 Struct containing the information of a GridapMHD gmsh model
#Fields
 -`geom_file`: Name of the gmsh geometry file located in the src/models/gmsh_files directory
 -`mesh_file`: Name of the gmsh mesh file located in the src/models/gmsh_files directory

 Provide either the geometry file to generate the mesh or the mesh itself. If both are provided, the provided mesh will be used and the geometry file will be ignored. If only the geometry file is provided, the mesh will be generated from it.

"""

struct gmsh_model 
    geom_file::String                #Name of the gmsh geometry file located in the src/models/gmsh_files directory
    mesh_file::String                #Name of the gmsh mesh file located in the src/models/gmsh_files directory   

    function gmsh_model(geom_file::String, mesh_file::String)  #Geom file is ignored if mesh file is provided
        endswith(mesh_file, ".msh")  || throw(ArgumentError("mesh_file field must end with .msh"))
        new(geom_file, mesh_file)
    end
end


function gmsh_model(file::String)

    if endswith(file, ".geo")  

        geo = joinpath(@__DIR__, "gmsh_files", file)
        msh = joinpath(@__DIR__, "gmsh_files", replace(file, ".geo" => ".msh"))

      gmsh.initialize()
        try
            gmsh.option.setNumber("General.Terminal", 1)
            gmsh.open(geo)                   # Loads and executes the .geo geometry
            gmsh.model.mesh.generate(3)       # Creates the 3D mesh
            gmsh.write(msh)  # Saves the mesh to a .msh file
        finally
            gmsh.finalize()
        end
      gmsh_model(file,replace(file, ".geo" => ".msh"))

    elseif endswith(file, ".msh")       
        gmsh_model("",file)

    else
        throw(ArgumentError("files field must be a path to a .geo or .msh file"))
    end
end


"""
Serial and parallel constructors for gmsh_model. The serial constructor takes the path to the mesh file and the geometry file (optional) and returns a gmsh_model instance. The parallel constructor takes a vector of partitions and the path to the mesh file and the geometry file (optional) and returns a gmsh_model instance.
"""
function (model::gmsh_model)() 

    msh = joinpath(@__DIR__, "gmsh_files", model.mesh_file)
    isfile(msh) || throw(ArgumentError("Mesh file $(msh) does not exist. Please provide a valid mesh file or generate the mesh from the geometry file."))

    return GmshDiscreteModel(msh)
end 

function (model::gmsh_model)(parts::Union{AbstractVector,MPIArray})

    msh = joinpath(@__DIR__, "gmsh_files", model.mesh_file)     
    isfile(msh) || throw(ArgumentError("Mesh file $(msh) does not exist. Please provide a valid mesh file or generate the mesh from the geometry file.")) 

    return GmshDiscreteModel(parts,msh)
end