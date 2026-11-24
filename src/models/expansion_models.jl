"""
Abstract supertype for all the expansion (the expansion along the direction of the B field) models to be run with this package
"""
abstract type mounted_expansion_models <:mounted_models end

struct expansion_gmsh_model <: mounted_expansion_models
   model::gmsh_model
   BCs::BC
   B::Union{VectorValue{3,Float64},Function}      
end

#Default instance for expansion_gmsh_model with B field in the y direction
expansion_gmsh_model(model::gmsh_model,BCs::BC) = expansion_gmsh_model(model,BCs,VectorValue(0.0,1.0,0.0))

#Serial and parallel constructors for expansion_gmsh_model. 

function (expansion_model::expansion_gmsh_model)() 
    
    model = expansion_model.model()
    tags = get_tag_names(model)
    false ∈ map(isequal,tags,expansion_model.BCs.tags.tags_J) && error("The defined BC tags do not match the model tags")
    
    return model
end

function (expansion_model::expansion_gmsh_model)(parts::Union{AbstractVector,MPIArray}) 
    
    model = expansion_model.model(parts)
    tags = get_tag_names(model)
    false ∈ map(isequal,tags,expansion_model.BCs.tags.tags_J) && error("The defined BC tags do not match the model tags")
    
    return model
end