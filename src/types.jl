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

"""
Struct with the FE spaces options. The constructor checks for only the following discretizations:
    For fluid: :Qk_dPkm1 (u in H1 and p in L2) and :RT (u in Hdiv and p in L2)
    For current: :H1 (φ in H1 with no j as independent variable) and :RT (j in Hdiv and φ in H1) 
 #Note
  GridapMHD accepts more discretizations, they could be included in this struct in future if necessary

 #Fields
 -`fluid_disc`: Symbol for the fluid discretization
 -`current_disc`: Symbol for the current discretization
 -`order_u`: Nusselt number
 -`order_j`: 
"""

struct FEspaces_options
    fluid_disc::Symbol
    current_disc::Symbol
    order_u::Integer
    order_j::Integer
    
    function FEspaces_options(fluid_disc,current_disc,order_u,order_j) 
        fluid_disc ∈ (:Qk_dPkm1, :RT) || throw(ArgumentError("Fluid discretization not implemented, only:Qk_dPkm1 (u in H1 and p in L2) and :RT (u in Hdiv and p in L2)"))
        current_disc ∈ (:H1, :RT) || throw(ArgumentError("Current discretization not implemented, only:H1 (φ in H1 with no j as variable) and :RT (j in Hdiv and φ in H1)"))
        new(fluid_disc,current_disc,order_u,order_j)
    end
end

#Default orders

FEspaces_options() = FEspaces_options(:Qk_dPkm1,:RT,2,1)

function FEspaces_options(fluid_disc,current_disc)
    if fluid_disc == :RT && current_disc == :H1
        order_u = 1
        order_j = 0
    elseif fluid_disc == :Qk_dPkm1 
        order_u = 2
        order_j = 1
    else
        throw(ArgumentError("No default FE spaces orders defined for that combination of discretization. Specify the orders"))
    end
    FEspaces_options(fluid_disc,current_disc,order_u,order_j)
end