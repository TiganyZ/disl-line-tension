module DislocationTypes

using Parameters
using InteractionTypes
using PeierlsPotential
import TrapSites: ConcSolutes, get_interaction_energy
# using Distributed

using StaticArrays
import ForwardDiff
import Base.LinAlg: gradient



export Disl_line, Dislocation, Analytic_disl_line, energy, gradient, hessian

struct Dislocation{T <: AbstractFloat} 
    K::T
    σ::StaticArrays.SArray{Tuple{3,3},T,2,9}
    b::StaticArrays.SArray{Tuple{3},T,1,3}
end

struct Disl_line
    d::Dislocation
    potential::Potential
    solutes::Union{Solutes,ConcSolutes}
end


function get_position_from_dofs(x)
    xy_pos = reshape( x, ceil(Int64, size(x,1)/2), 2 )'

    N = size(xy_pos,2)

    positions = zeros(eltype(x), 3, N)
    positions[1:2,:] .= xy_pos
    return positions
end


function energy( d::Disl_line, x)
    return energy_function_of_positions( d, x, d.potential.ΔEₚ )
end

function gradient( V::Disl_line, x)
    return ForwardDiff.gradient(y->energy(V,y), x)
end

function hessian( V::Disl_line, x)
    return ForwardDiff.hessian(y->energy(V,y), x)
end


# function gradient(d::Disl_line, x)
#     return vcat(energy_function_of_positions( d, x, d.∂ΔEₚ[1], true, 1 ),
#                 energy_function_of_positions( d, x, d.∂ΔEₚ[2], true, 2 ))
# end

function get_interaction_energy(solutes::Solutes, j, Pⱼ, derivative=false, hessian=false, direction=1, direction2=1)
    E_int = 0.0
    b_mag = 2.87 * √3 / 2 
            
    for idx in 1:size(solutes.positions,2)

        if abs.( j*b_mag .- solutes.positions[3,idx] ) .< 1e-3

            dv = Pⱼ[1:2] .- solutes.positions[1:2,idx]
            dist = norm( dv  )

            if derivative
                dd = get_dd(direction, direction2, false)
                #         meV       eV/b
                E_int += 1000. * dlorentzian(solutes.interaction_type,  dist / b_mag)*dd( dv[1]/b_mag, dv[2]/b_mag )  
            elseif hessian
                dd = get_dd(direction, direction2, true)
                E_int += 1000. *  ddlorentzian(solutes.interaction_type,  dist / b_mag)*dd( dv[1]/b_mag, dv[2]/b_mag )  
            else                           
                E_int += 1000. *   lorentzian(solutes.interaction_type,  dist / b_mag) 
            end

        end
    end
    
    return E_int
end

function get_dd(direction, direction2, hessian=false)

    if hessian
        if direction == 1
            if direction2 == 1
                dd = (x,y) -> 1/√(x^2+y^2) -x^2 / (x^2 + y^2)^(3/2)
            else
                dd = (x,y) -> -y*x / (x^2 + y^2)^(3/2)
            end
        else
            if direction2 == 2
                dd = (x,y) -> 1/√(x^2+y^2) -y^2 / (x^2 + y^2)^(3/2)
            else
                dd = (x,y) -> -y*x / (x^2 + y^2)^(3/2)
            end
        end
    else
        if direction == 1
            dd = (x,y) -> x/√(x^2+y^2)
        else
            dd = (x,y) -> y/√(x^2+y^2)
        end
    end
    return dd
end

function energy_function_of_positions(d::Disl_line, x, f::Function)
    # positions = get_position_from_dofs(x)
    # N = size( positions, 2 ) 

    N = ceil(Int64, size(x,1)/2)

    
    # Conversion to units of Å
    E_line = 0.0 # zeros(eltype(x), N+1 )    

    # E_line = @parallel (+) for j in 1:N
    #     energy_point(d, x, j, N)
    # end
 
    
    for j in 1:N
        E_line += energy_point(d, x, j, N)
    end
    
    return E_line
end



function get_energy_quantities(x, j, N)
    conv = √2 / 3 * 2.87
    # Can access only x at the indices and make static vectors from them
    # > x,y at indices j, j+N
    xi = j
    yi = j+N

    Pⱼ  = SVector( x[xi]*conv, x[yi]*conv, 0.0 )
    if j == N
        Pₖ = SVector( (x[xi] - x[xi-1])*conv + Pⱼ[1], (x[yi] - x[yi-1])*conv + Pⱼ[2] , 0)
    else    
        Pₖ = SVector( x[xi+1]*conv, x[yi+1]*conv, 0 ) #conv .* positions[:,j + 1]
    end
    ΔPⱼₖ = ( Pⱼ - Pₖ )
    b_mag = 2.87 * √3 / 2 

    if j == 1
        mag =  sqrt(  ((x[xi+1] - x[xi])*conv)^2 +  ((x[yi+1] - x[yi])*conv)^2 +  b_mag^2 )
        l = SVector( (x[xi+1] - x[xi])*conv/mag, (x[yi+1] - x[yi])*conv/mag, b_mag/mag )
    else
        mag =  sqrt(  ((x[xi] - x[xi-1])*conv)^2 +  ((x[yi] - x[yi-1])*conv)^2 +  b_mag^2 )
        l = SVector(  (x[xi] - x[xi-1])*conv/mag, (x[yi] - x[yi-1])*conv/mag, b_mag/mag) 
    end
    
    return ΔPⱼₖ, Pⱼ, l
end


# function get_energy_quantities(positions, j, N)
#     conv = √2 / 3 * 2.87
#     if j == N
#         Pⱼ  = @SVector conv .* positions[:,j]
#         Pₖ = @SVector conv .* (positions[:,j] - positions[:,j - 1]) + Pⱼ
#     else    
#         Pⱼ  = @SVector conv .* positions[:,j]
#         Pₖ = @SVector conv .* positions[:,j + 1]
#     end
#     ΔPⱼₖ = @SVector ( Pⱼ - Pₖ )
#     b_mag = 2.87 * √3 / 2 

#     if j == 1
#         l = @SVector conv .* (positions[:,j+1] - positions[:,j]) .+ [ 0, 0., b_mag] ./ norm(conv .* (positions[:,j+1] - positions[:,j]) .+ [ 0, 0., b_mag])
#     else    
#         l = @SVector conv .* (positions[:,j] - positions[:,j-1]) .+ [ 0, 0., b_mag] ./ norm(conv .* (positions[:,j] - positions[:,j-1]) .+ [ 0, 0., b_mag])
#     end
    
#     return ΔPⱼₖ, Pⱼ, l
# end

function energy_point(D::Disl_line, x, j, N, detail=false)
    ΔPⱼₖ, Pⱼ, l = get_energy_quantities(x, j, N)    
    conv = √2 / 3 * 2.87
    if D.solutes.interact
        if isa(D.solutes, ConcSolutes)
            if detail
                # in simulation       ip2_x      = √3
                # In trap site stuff, ip2_trap_x = (2/6. * √2 * 2.87 * √3)
                # Therefore, to go from sim to trap site in angstrom 
                E_int = get_interaction_energy(D.solutes, conv * Pⱼ[1:2], false)
            else
                E_int = get_interaction_energy(D.solutes, Pⱼ[1:2], true)
            end
        else
            E_int = get_interaction_energy(D.solutes, j, Pⱼ, false, false)
        end
    end

    E_line = D.d.K/2. * ( ΔPⱼₖ ⋅ ΔPⱼₖ )  + D.potential.ΔEₚ( (Pⱼ[1:2]./conv)... ) + 1000*cross( D.d.σ * D.d.b, l ) ⋅ Pⱼ - E_int

    if detail
        return [ D.d.K/2. * ( ΔPⱼₖ ⋅ ΔPⱼₖ ),  
                 D.potential.ΔEₚ( (Pⱼ[1:2]./conv)... ), 
                 1000*cross( D.d.σ * D.d.b, l ) ⋅ Pⱼ,
                 -E_int, E_line ]
    else
        return E_line
    end
end


# function gradient_chunk

function gradient_interaction_energy(S::ConcSolutes, core_position, direction, finite=true, s=1e-3, order=:second)
    # ForwardDiff.gradient(y->get_interaction_energy(S,y), core_position)
  
    # > When ForwardDiff gets to the concentration spline, it
    # > necessarily stops, as the Dierckx spline implementation is
    # > based on fortran.
    # >> Solution: Use crude finite differences
    
    if finite
        h = direction == 1 ? [ s, 0 ] : [ 0, s ]
    
        if order == :second
            return (get_interaction_energy(S,core_position + h) - get_interaction_energy(S,core_position - h)) / (2*s)
        else
            return (get_interaction_energy(S,core_position + h) - get_interaction_energy(S,core_position)) / s
        end
    else
        return ForwardDiff.gradient(y->get_interaction_energy(S,y), core_position) 
    end
end


function gradient_point(D::Disl_line, x, j, N, direction)
    ΔPⱼₖ, Pⱼ, l = get_energy_quantities(x, j, N)    

    #  from internal units to Å
    conv = √2 / 3 * 2.87
    b_mag = √3/2 * 2.87
    
    if D.solutes.interact

        if isa(D.solutes, ConcSolutes)
            E_int = gradient_interaction_energy(D.solutes, conv * Pⱼ[1:2], direction)
        else
            E_int = get_interaction_energy(D.solutes, j, Pⱼ, true, false, direction) / (b_mag / conv) #/ b_mag / conv
        end
        # This is in meV/b therefore multiply by 1/b -> 1/(√3/2 * a)
        #        println("Interaction grad = $E_int meV/")
    end
    #              meV / Å
    #   i <--> j <--> k are the points on the line
    if j > 1
        ΔPᵢⱼ, Pᵢ, lᵢ = get_energy_quantities(x, j-1, N)    
        elastic = D.d.K * ( ΔPⱼₖ[direction]*conv - ΔPᵢⱼ[direction]*conv )
    else
        elastic = D.d.K * ( ΔPⱼₖ[direction]*conv  )
    end
    
    E_line =  elastic + D.potential.∂ΔEₚ[direction]( (Pⱼ[1:2]./conv)... ) + 1000*cross( D.d.σ * D.d.b, l )[direction]*conv - E_int

    
    return E_line
end

function hessian_point(D::Disl_line, x, j, k, N, xi, xj)
    ΔPⱼₖ, Pⱼ, l = get_energy_quantities(x, j, N)    
    conv = √2 / 3 * 2.87
    b_mag = √3/2 * 2.87

    E_int = 0.0
    if D.solutes.interact
        if j == k 
            E_int = get_interaction_energy(D.solutes, j, Pⱼ, false, true, xi, xj)  / (b_mag / conv)^2
        end
    end
    conv = √2 / 3 * 2.87

    elastic = 0.0
    
    if (k == j)
        if j == 1 || j == N
            elastic = D.d.K * conv^2
        else
            elastic = 2 * D.d.K * conv^2
        end
    end
    
    if (k == j-1) elastic = -D.d.K * conv^2 end
    if (k == j+1) elastic = -D.d.K * conv^2 end

    peierls = 0.0

    if j == k
        peierls += D.potential.hess[xi,xj]( (Pⱼ[1:2]./conv)... )
    end
    
    E_line = elastic + peierls - E_int
    # println("Hessian j,k = $j, $k  : elastic = $elastic, peierls = $peierls, int = -$E_int")
    return E_line
end



function construct_gradient(d::Disl_line, x, write)
    N = ceil(Int64, size(x,1)/2)

    # E_line = SharedArray{Float64,1}(2N)
    # xi = SharedArray{Float64,1}(x)


    # # Think about whether to do a shared array each of the data points. 
    # # > Data references would have to be thought of 
    # # > sdata() ??
    # # use @parallel, time it??

    # @parallel for j in 1:2N
    #     E_line[j] = gradient_point(d, x, (j-1)%N + 1, N, ceil(Int64, j/N))
    # end

    E_line = zeros(2N)
     #    xi = SharedArray{Float64,1}(x)

    if write 
        for j in 1:2N
            E_line[j] = gradient_point(d, x, (j-1)%N + 1, N, ceil(Int64, j/N))
        end
    end
    #E_line = SVector{2N}( [gradient_point(d, x, (j-1)%N + 1, N, ceil(Int64, j/N)) for j in 1:2N] )
    return E_line
end

# ##-- These functions show some of the behaviour
# function pmap(f, lst)
#     np = nprocs()  # determine the number of processes available
#     n = length(lst)
#     results = Vector{Any}(n)
#     i = 1
#     # function to produce the next work item from the queue.
#     # in this case it's just an index.
#     nextidx() = (idx=i; i+=1; idx)
#     @sync begin
#         for p=1:np
#             if p != myid() || np == 1
#                 @async begin
#                     while true
#                         idx = nextidx()
#                         if idx > n
#                             break
#                         end
#                         results[idx] = remotecall_fetch(f, p, lst[idx])
#                     end
#                 end
#             end
#         end
#     end
#     results
# end

# @everywhere function myrange(q::SharedArray)
#     idx = indexpids(q)
#     if idx == 0 # This worker is not assigned a piece
#         return 1:0, 1:0
#     end
#     nchunks = length(procs(q))
#     splits = [round(Int, s) for s in Linspace(0,size(q,2),nchunks+1)]
#     1:size(q,1), splits[idx]+1:splits[idx+1]
# end

# @everywhere function advection_chunk!(q, u, irange, jrange, trange)
#     @show (irange, jrange, trange)  # display so we can see what's happening
#     for t in trange, j in jrange, i in irange
#         q[i,j,t+1] = q[i,j,t] + u[i,j,t]
#     end
#     q
# end

# @everywhere advection_shared_chunk!(q, u) =
#            advection_chunk!(q, u, myrange(q)..., 1:size(q,3)-1)

function get_xi_xj(k, N, left=true)
    if left
        if k > N
            xi, xj = 2, 1
        else
            xi, xj = 1, 1
        end
    else
        if k > N
            xi, xj = 2, 2
        else
            xi, xj = 1, 2
        end
    end
    return xi, xj
end
function get_hessian_row!(d, k, E_line, N, x)

    left=true
    xi, xj = get_xi_xj(k, N, left)
    for j in 1:N
        E_line[j, k] = hessian_point(d, x, j, k, N, xi, xj)
    end

    left=false
    xi, xj = get_xi_xj(k, N, left)
    for j in 1:N
        E_line[N+j, k] = hessian_point(d, x, j, k, N, xi, xj)
    end

end

function construct_hessian(d::Disl_line, x)
    # positions = get_position_from_dofs(x)
    # N = size( positions, 2 ) 

    N= ceil(Int64, size(x,1)/2)
    # Conversion to units of Å
    E_line = zeros(eltype(x), 2*N, 2*N )    

    conv = √2 / 3 * 2.87
    b_mag = √3/2 * 2.87


    for k in 1:2N
        if k > N xk = 2 else xk = 1 end

        for j in 1:2N
            if j > N xj = 2 else xj = 1 end

            if abs(k - j) <= 1
                E_line[j,k] = hessian_point(d, x, (j-1)%N +1, (k-1)%N +1, N, xj, xk)
            end
        end
    end

    return E_line
end
 


#--------ANALYTIC DISLOCATION LINE---------#

@with_kw type Analytic_disl_line
   h::Float64 = √3
   E₀::Float64 = 0.0
   ΔE::Float64 = 100.
   b::Float64 = 1.
   τₚ::Float64 = π * ΔE / ( h * b )
   pot_type::String = "Sinusoidal"
end


function energy(d::Analytic_disl_line, x , get_arr=false)
    # Using Dorn and Rajnak Sinusoidal
    
    # Uₖ = ( 2^(3/2) / π ) * h * (E₀ * ΔE)^(0.5)
    # > h  = distance between peierls valleys
    # > ΔE = Difference between E₀ and maximal Eₚ.
    # > E₀ = Initial energy in Peierls valley

    xy_pos = reshape( x, ceil(Int64, size(x,1)/2), 2 )

    N = size(xy_pos,1)

    positions = zeros(N,3)
    positions[:,1:2] .= xy_pos

    E_line =  zeros( N )

    
    if d.pot_type == "Sinusoidal"
    
        f = (d, y) -> d.E₀ + (d.ΔE / 2) * ( 1 - cos( 2*π*y / d.h) ) 
        
    elseif d.pot_type == "Eshelby"
    
        f = (d, y) -> d.E₀ + 16 * d.ΔE  * (y/d.h)^2 * ( 1 - y/ d.h )^2 
    
    end

    E_line = populate_energy_array!(d, f, N, E_line, xy_pos[1:N,1])


    println(" \n ---->>>>####    Energy   = ", sum(E_line) ," meV\n")
    if get_arr
        return E_line
    else
        return sum(E_line)
    end
    
end

function gradient( d::Analytic_disl_line, x )
    # Get the line derivative in the x direction
    # > positions: for core position each slab, and it is an array of 1-D arrays
    # > l is vector parallel to dislocation line
    # > ΔEₚ is the peierls energy function
    positions = reshape( x, ceil(Int64, size(x,1)/2), 2 )    

    N = size(positions,1)
    dE = zeros(size(positions,1),2)

    # println("\n--->>>###   Positions   ###<<<--- ", d.pot_type)
    # println(positions, "\n")
    
    
    if d.pot_type == "Sinusoidal"
    
        df = (d, y) ->  (d.ΔE * π / d.h) * ( sin( 2*π*y / d.h) ) 
        
    elseif d.pot_type == "Eshelby"
    
        df = (d, y) ->  16 * d.ΔE  *  2 / d.h * (  (y/d.h) * ( 1 - y/ d.h )^2 +  -(y/d.h)^2 * ( 1 - y/ d.h ) )
    
    end


    direction = 1
    dE[:,direction] = populate_energy_array!(d, df, N, dE[:,direction], positions[1:N,1])
    dE[:,2] .= 0.

    return vcat(dE...)
    
end


function populate_energy_array!(d, f, N, E_line, y)
    # In this case, we are taking y to be as in Caillard
    # > This means x → y, which is the axis along which the peierls potential varies
    
    for j in 1:N
        println("yi = ", y[j])
        E_line[j] = f(d, y[j])
    end
    return E_line
end


end
