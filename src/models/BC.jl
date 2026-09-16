"""
Struct with the tags of the Dirichlet boundary conditions of a channel model. 

#Fields
 -`U`: Names of the velocity boundaries
 -`J`: Names of the electric current boundaries
 -`φ`: Names of the electric potential boundaries
"""

struct BC_tags
    tags_U::Vector{String}
    tags_J::Vector{String}
    tags_φ::Vector{String}
end

BC_tags(tag_U,tag_J) = BC_tags(tag_U,tag_J,String[])

"""
Struct with the values of the Dirichlet boundary conditions of a channel model

#Fields
 -`U`: Values of of the velocity boundaries
 -`J`: Values of the electric current boundaries
 -`φ`: Values of the electric potential boundaries
"""

struct BC_values
    U::Vector{Union{VectorValue{3,Float64},Function}}
    J::Vector{Union{VectorValue{3,Float64},Function}}
    φ::Vector{Union{Float64,Function}}
end

BC_values(U,J) = BC_values(U,J,Float64[])
BC_values(U) = BC_values(U,VectorValue{3,Float64}[])

"""
Struct containing the BC tags and BC values 
It auto-completes with zeroes the missing values

"""

struct BC
    tags::BC_tags
    values::BC_values
    function BC(tags,values)
        dU = length(tags.tags_U)-length(values.U)
        dJ = length(tags.tags_J)-length(values.J)
        dφ = length(tags.tags_φ)-length(values.φ)
    
        new_U= dU !== 0 ? [values.U...,fill(VectorValue(0.0,0.0,0.0),dU)...] : values.U  
        new_J= dJ !== 0 ? [values.J...,fill(VectorValue(0.0,0.0,0.0),dJ)...] : values.J  
        new_φ = dφ !== 0 ? [values.φ...,fill(0.0,dφ)...] : values.φ  
        
        new_values = BC_values(new_U,new_J,new_φ)

        new(tags,new_values)
    end
end