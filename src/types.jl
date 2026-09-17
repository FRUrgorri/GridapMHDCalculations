#Data structures to work with the project

"""
Supertype for all the models to be run with this package
"""
abstract type mounted_models end

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


function Dimensionless_numbers(; Ha=nothing, Re=nothing, N=nothing) #Method based on keyword arguments 

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