#Data structures to work with the project

"""
Supertype for all the models to be run with this package
"""
abstract type mounted_models end

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
Struct with the dimensionless numbers. The constructor checks if the numbers are consistent according to the equations: N = Ha²/Re

#Fields
 -`Ha`: Hartmann number
 -`Re`: Reynolds number
 -`N`: Nusselt number
"""

struct Dimensionless_numbers{T<:Real}
    Ha::T
    Re::T
    N::T
    Dimensionless_numbers{T}(Ha,Re,N) where{T<:Real} = N ≈ Ha^2/Re ? new(Ha,Re,N) : ArgumentError("Inconsisten choice of dimensionless numbers: require N = Ha²/Re")
end
Dimensionless_numbers(Ha,Re,N) = Dimensionless_numbers{Float64}(Ha,Re,N)    #Default method for no explicit type definition, it promotes integers to float by default
Dimensionless_numbers() = Dimensionless_numbers(1.0,1.0,1.0)


function Dimensionless_numbers(; Ha=nothing, Re=nothing, N=nothing) 

    n = count(!isnothing, (Ha, Re, N))
    if n == 2
        if isnothing(N)
            N = Ha^2 / Re
        elseif isnothing(Re)
            Re = Ha^2 / N
        else  # Ha missing
            Ha = sqrt(N * Re)
        end
    elseif n == 3
        Ha^2 / Re ≈ N || throw(ArgumentError(
            "Inconsisten choice of dimensionless numbers: require N = Ha²/Re"))
    else
        throw(ArgumentError("provide two of Ha, Re, N (or all three consistently)"))
    end
    return Dimensionless_numbers(Ha, Re, N)
end