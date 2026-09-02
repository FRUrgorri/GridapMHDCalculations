#Data structures to work with the project

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
  order::Integer
end

#Constructors of the output_info type

output_info(xh,Ω,B,path,title) = output_info(xh,Ω,B,path,title,2)

"""

"""

"""
struct Dirichlet_tag_names 
    U::Vector{String}
    J::Vector{String}
    φ::Vector{String}
    let allowed = ["insulated", "inlet", "outlet", "conducting", "thin_wall"]
        Dirichlet_tag_names(U,J,φ) = map(U) do tag 
            tag ∈ allowed ? new(U,J,φ) : error("Tag name is not allowed as boundary name:", tag) 
        end
        Dirichlet_tag_names(U,J,φ) = map(J) do tag 
            tag ∈ allowed ? new(U,J,φ) : error("Tag name is not allowed as boundary name:", tag) 
        end
        Dirichlet_tag_names(U,J,φ) = map(φ) do tag 
            tag ∈ allowed ? new(U,J,φ) : error("Tag name is not allowed as boundary name:", tag) 
        end
    end
end
"""
