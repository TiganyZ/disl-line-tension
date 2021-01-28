module InteractionTypes

using Parameters
using Dierckx

export Solutes, ConcSolutes, C_Lorentzian, H_Lorentzian, lorentzian, dlorentzian, ddlorentzian


@with_kw struct C_Lorentzian{T <: AbstractFloat}
   γ::Float64 = 0.37126788430432184
   r₀::Float64 = 6.514357682207945
end

@with_kw struct H_Lorentzian{T <: AbstractFloat}
   A::Float64 = 390./1000
   # Making the below acceptable to take units of burgers vector
   r₀::Float64 = (√6 * 2.87 / 3) / ( √3 / 2 * 2.87 )
end

struct Solutes{T <: AbstractFloat}
    interact::Bool
    interaction_type::Union{C_Lorentzian{T},H_Lorentzian{T}}
    positions::Array{Float64,2}
end

struct ConcSolutes{T <: AbstractFloat}
    interact::Bool
    interaction_type::Union{C_Lorentzian{T},H_Lorentzian{T}}
    convert_sitelabel::Function
    conc_func::Dierckx.Spline1D
    ref_conc_sum::Float64
    T::Float64
end


# ----- CARBON LORENTZIAN ----- #

function lorentzian(d::C_Lorentzian, x)
    sq = 1. + ( x / d.γ / d.r₀ )^2
    
    l = 1. / ( π * d.γ * sq )
    #    println(" >> Lorentzian: $l")
    return l
end

function dlorentzian(d::C_Lorentzian, x)
    sq = 1. + ( x / d.γ / d.r₀ )^2
    
    return -(2 * x / (d.γ * d.r₀)^2 ) / (π * d.γ) * (1/sq^2)    

end

function ddlorentzian(d::C_Lorentzian, x)
    sq  = 1. + ( x / (d.γ * d.r₀) )^2
    dsq =    2 * x / (d.γ * d.r₀)^2
    ddsq =   2     / (d.γ * d.r₀)^2

    u  = -  dsq / (π * d.γ)
    du = - ddsq / (π * d.γ)

    v  = 1/sq^2
    dv = (-2/sq^3) * dsq
    
    return  u*dv + v*du 

end



# ----- ITAKURA HYDROGEN LORENTZIAN ----- #

function lorentzian(d::H_Lorentzian, x)

    # Factor of 1000 is to keep it in eV 
    l = d.A / ( 1. + 2*( x / d.r₀ )^2 )
    #    println(" >> Lorentzian: $l")
    return l
end

function dlorentzian(d::H_Lorentzian, x)
    

    return -( 4*x / d.r₀^2 ) * d.A / ( 1. + 2*( x / d.r₀ )^2 )^2
end

function ddlorentzian(d::H_Lorentzian, x)
    sq  = ( 1. + 2*( x / d.r₀ )^2 )
    dsq =    4 * x / ( d.r₀)^2
    ddsq =   4     / ( d.r₀)^2

    u  = -  dsq 
    du = - ddsq 

    v  =  d.A/sq^2
    dv = (-2 * d.A/sq^3) * dsq
    
    return  u*dv + v*du 

end
end
