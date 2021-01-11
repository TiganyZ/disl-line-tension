module McCleanIsotherm

include("interaction_types.jl")

using Plots
using NLsolve
using Dierckx
using InteractionTypes

export get_concentration_distance_dependence

"""
This solves for the equilibrium carbon/solute concentration within a
thermodynamical mean-field model.

From Ventelon 2015, Dislocation core reconstruction induced by carbon
segregation in bcc iron, 10.1103/PhysRevB.91.220102

Average carbon-dislocation interaction energy given as
> E_int(Cₖ) = E⁰_int + ΔE_easy_hard / Cₖ + CₖV_CC
>> Cₖ = Average carbon concentration on the dislocation core

> Cₖ / ( 1 - Cₖ ) = C_bulk / (1 - C_bulk) * exp( - E_seg( Cₖ / kbT ) )   (2)
> E_seg(Cₖ) = E_int + Cₖ ∂E_int/∂Cₖ
            = E⁰_int + 2CₖV_CC

>> C_bulk = Carbon concentration in the matrix
>> C_nom  = Nominal concentration of carbon atoms per iron atom

Nₖ sites = ρV/b
>> ρ  = dislocation density
>> N₀ = 6V/a³

∴ by matter conservation
>> N₀C_bulk + NₖCₖ = N₀C_nom / 3    (3)

This can be used to solve (2) self-consistently to obtain the carbon
concentration on the dislocation lines as a function of the nominal
carbon concentration and dislocation density

Inputs:
-- C_nom : Nominal carbon concentration
-- ρ     : Dislocation density

Outputs:
-- Cₖ    : Concentration of carbon around dislocation line
-- C_bulk: Carbon concentration around bulk

"""

# export get_concentration_distance_dependence

struct Seg{T <: AbstractFloat}
    E⁰::T
    V_CC::T
    T::T
    kb::T
end

struct DislSystem{T <: AbstractFloat}
    ρ::T
    V::T
    a::T
    b::T
    C_nom::T
end


N₀(D::DislSystem) = 6*D.V / D.a^3
Nₖ(D::DislSystem) = D.ρ * D.V / D.b


function E(p::Seg, Cₖ::Float64)
    global verb
    if verb > 2 println("E Seg: E⁰       = ", p.E⁰) end
    if verb > 2 println("E Seg: 2Cₖ*V_CC = ", 2Cₖ*p.V_CC) end
    if verb > 0 println("E Seg: Total    = ", p.E⁰ + 2Cₖ*p.V_CC) end
    return p.E⁰ + 2Cₖ*p.V_CC
end

function Cₖ_frac(p::Seg, C_bulk, Cₖ)
    global verb
    if verb > 0
        in_exp = - E(p, Cₖ) / ( p.kb * p.T )
        expon  = exp( in_exp )
        coeff  = C_bulk / (1 - C_bulk)
        if verb > 2 println("Cₖ frac exp: kbT    = $(( p.kb * p.T ))  ") end        
        if verb > 2 println("Cₖ frac exp: in_exp = $in_exp  ") end
        if verb > 2 println("Cₖ frac exp:    exp = $expon  ") end
        if verb > 2 println("Cₖ frac exp:  coeff = $coeff  ") end
        if verb > 0 println("Cₖ frac exp:  total = $(coeff*expon)  ") end
    end
    return C_bulk / (1 - C_bulk) * exp( - E(p, Cₖ) / ( p.kb * p.T ) )
end


function get_Cₖ(p::Seg, C_bulk, Cₖ)
    global verb
    if verb > 0
        if verb > 2 println("Cₖ Seg from frac:  Cₖ_frac = $(Cₖ_frac(p, C_bulk, Cₖ))   ") end
        println("Cₖ Seg from frac:  Cₖ      = $(Cₖ_frac(p, C_bulk, Cₖ) / (1. + Cₖ_frac(p, C_bulk, Cₖ)))   ")
    end
    Cₖf = Cₖ_frac(p, C_bulk, Cₖ)
    if isinf( Cₖf )
        return 1.0
    else
        return Cₖf / (1. + Cₖf)
    end
end


function C_bulk_frac(p::Seg, C_bulk, Cₖ)
    global verb
    if verb > 0
        in_exp =  E(p, Cₖ) / ( p.kb * p.T )
        expon  = exp( in_exp )
        coeff  = C_ₖ / (1 - Cₖ)
        if verb > 2 println("Cₖ frac exp: kbT    = $(( p.kb * p.T ))  ") end                
        if verb > 2 println("C_bulk frac exp: in_exp = $in_exp  ") end
        if verb > 2 println("C_bulk frac exp:    exp = $expon  ") end
        if verb > 2 println("C_bulk frac exp:  coeff = $coeff  ") end
        println("C_bulk frac exp:  total = $(coeff*expon)  ") 
    end
    C_bulkf = C_bulk_frac(p, C_bulk, Cₖ)
    if isinf( C_bulkf )
        return 1.0
    else
        return C_bulkf / (1. + C_bulkf)
    end    
end

function get_C_bulk(p::Seg, C_bulk, Cₖ)
    global verb
    if verb > 0
        if verb > 2 println("C_bulk Seg from frac:  C_bulk_frac = $(C_bulk_frac(p, C_bulk, Cₖ))   ") end
        println("C_bulk Seg from frac:  C_bulk      = $(C_bulk_frac(p, C_bulk, Cₖ) / (1. + C_bulk_frac(p, C_bulk, Cₖ)))   ")
    end
    return C_bulk_frac(p, C_bulk, Cₖ) / (1. + C_bulk_frac(p, C_bulk, Cₖ))
end


function get_C_bulk(D::DislSystem, Cₖ)
    global verb
    if verb > 0
        if verb > 2 println("get_C_bulk cons.: D.C_nom / 3        = ", D.C_nom / 3) end
        if verb > 2 println("get_C_bulk cons.: - Nₖ(D)/N₀(D) * Cₖ = ", - Nₖ(D)/N₀(D) * Cₖ) end
        println("get_C_bulk cons.: C_bulk             = ", ( N₀(D)*D.C_nom / 3 - Nₖ(D)*Cₖ     ) / N₀(D))        
    end
    return ( N₀(D)*D.C_nom / 3 - Nₖ(D)*Cₖ     ) / N₀(D)
end

function get_Cₖ(D::DislSystem, C_bulk)
    global verb
    if verb > 0
        if verb > 2 println("get_Cₖ cons.:   N₀(D)/Nₖ(D) * D.C_nom / 3 = ",  N₀(D)*D.C_nom / 3  / Nₖ(D)) end
        if verb > 2 println("get_Cₖ cons.: - N₀(D)/Nₖ(D) * C_bulk    = ", - N₀(D)/Nₖ(D) * C_bulk) end
        println("get_Cₖ cons.:  Cₖ   = ", ( N₀(D)*D.C_nom / 3 - N₀(D)*C_bulk ) / Nₖ(D)) 
    end
    return ( N₀(D)*D.C_nom / 3 - N₀(D)*C_bulk ) / Nₖ(D)
end



####----     Solver for the fixed point solutions    ----####

function sc_iteration!(resvec, C, S, D)
    Cₖ     = get_Cₖ(S, C[2], C[1])
    C_bulk = get_C_bulk(D, Cₖ)
    
    resvec[1] = Cₖ     - C[1]
    resvec[2] = C_bulk - C[2]
end

function solve_nl(D::DislSystem, S::Seg)
    C⁰_bulk = D.C_nom/3/2 
    f_sc! = (resvec, C) -> sc_iteration!(resvec, C, S, D)

    x0 = [ get_Cₖ(D, C⁰_bulk), C⁰_bulk]
    res = nlsolve( f_sc!, x0; method=:anderson, iterations=10_000, m=20, beta=1e-5 )
    return res
end

###---      Obtaining distance dependence of self-consistent concentration      ---###

function get_concentration_data(C_nomi, ρi, energies, T, V_CC_ventelon)
    Å =  1e-10
    a = 2.87Å
    b = √3/2 * a
    kb = 0.000086173324 # eV/K

    all_Cₖ = zeros(Float64, length(T), length(energies))
    D = DislSystem( ρi, 1., a, b, C_nomi )
    
    for (en,E⁰) in enumerate(energies)
        for (tn,Ti) in enumerate(T)
            S = Seg(E⁰, V_CC_ventelon, Ti, kb )
            # C = solve_nl(D, S, D.C_nom/3/2 ).zero[1]
            # if tn > 3:
            #     if abs( C - mean(all_Cₖ[tn-3:tn, en, ρn, cn]) ) 
            all_Cₖ[tn, en] = solve_nl(D, S).zero[1]
        end
    end
    return all_Cₖ
end

function clean_concentration_data(all_Cₖ, energies, T)
    Ckk = copy(all_Cₖ)
    for (en,E⁰) in enumerate(energies)
        for tn in length(T):-1:1 
            if length(T)-tn >= 2
                grad   = ( Ckk[tn+1, en] -  Ckk[tn  , en]  ) / (T[tn+1]-T[tn])
                gradp1 = ( Ckk[tn+2, en] -  Ckk[tn+1, en]  ) / (T[tn+1]-T[tn])                    
                if abs(gradp1 - grad) > 0.001
                    Ckk[tn, en] .= Ckk[tn+1, en] 
                end
            end
        end
    end
    return Ckk
end

####----     Interpolate between the energy points    ----####

function get_concentration_distance_dependence(C_nom, ρ, solute_int)
    # > We have the binding energy of carbon as a function of distance from the dislocation core.
    # > One can sample from this distribution at multiple points and obtain the concentration as a function of distance at a certain temperature.
    # > Operating temperature is 320K

    # We now
    # > Obtain the binding energy at a particular distance and dislocation density
    # > Get data for the whole temperature range
    # > Clean data
    # > Extract data at T = 320°K
    # > Find concentration as a function of distance from the dislocation core (ignoring intersite interactions)
    
    global verb
    const verb = 0
    
    const Å =  1e-10
    const a = 2.87Å
    const b = √3/2 * a
    const kb = 0.000086173324 # eV/K
    
    const maxiter = 1000
    const tol = 1e-3
    
    const E⁰_ventelon = -0.84  # eV
    const V_CC_ventelon = 0.21 # eV
    
    const Tf = 1200
    T = collect(Float64, 200:10:Tf) 
    b_mag = 2.87 * √3 / 2
    # Want range which is up to 2.5 burgers vectors
    distances = collect(0:0.1:2.5) .* b_mag
    energies = [ -lorentzian( solute_int, d/b_mag ) for d in distances ]
    
    all_Cₖ = get_concentration_data(C_nom, ρ, energies, T, V_CC_ventelon)
    Ckk = clean_concentration_data(all_Cₖ, energies, T)

    # Now interpolate the cleaned data at T = 320
    op_temp = 320.
    tidx = find( x -> x==op_temp, T  )

    # We now have the concentrations as a function of distance, as there is a one-to-one correspondence between the energy and the distance
    # > We can interpolate this function as a spline and return the function.
    # > Would make sense to have this in interaction types for clarity and brevity, where we can export the type and define the gradient as such. 

    # This is defined where the distances are in angstrom, so no conversion
    xi = vcat(Ckk[ tidx, : ]...)
    S  = Spline1D(distances, xi, w = ones(length(xi)), k = 3, bc = "error")
    dS = x -> derivative(S, x)



    return S, dS, Ckk, T, distances, energies

end


function plot_conc_vs_temperature(T, Ckk, energies, plot_type)

    if plot_type == "Easy"
        labels = ["E1 ", "E2 ", "E3 ", "E4 ", "E5 ", "E6 ", "E7 ", "E8 ", "E9 ", "E10"]
    elseif plot_type == "Hard"
        labels = [   "H1 ", "H2 ", "H3 ", "H4 ", "H5 ", "H6 ", "H7 "]
    end


    pyplot( xlims = (0, Tf),
            ylims = (0, 1),
            size=(800,600),
            xticks=collect(0:200:Tf),
            yticks=collect(0:0.2:1.2), legend=false)

   fnt = Plots.font( "Helvetica", 30 )
   default(titlefont = fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

   colors = palette(:tab10)
   ls = [:solid, :dash, :dot, :dashdot]

   ###---   For the easy core
    for i in 1:length(energies[1:10])
        if i == j == 1
            pp = plot(T, Ckk[:,i,j,k], label = labels[i], linestyle=ls[j], linewidth=3, show=false, linecolor=colors[i])
            xlabel!("T[K]")
            ylabel!("Cd")
        elseif j == 1
            plot!(T, Ckk[:,i,j,k], label = labels[i], linestyle=ls[j], linewidth=3, show=false, linecolor=colors[i])
            # elseif i == length(energies) && j == length(ρ)
            #     break
        else
            plot!(T, Ckk[:,i,j,k], linestyle=ls[j], linewidth=3, show=false, label="", linecolor=colors[i])
        end
    end
end


function make_dir_if_nonexist(dirname)
    if !isdir(dirname)
        mkdir(dirname)
    end
end

function collect_concentrations(C_nom, ρi, solute_int)
    # Hardcoding the length the concentrations == length distances + 1...
    C = zeros(length(C_nom), 26)

    for (i,C_nomi) in enumerate(C_nom)
        S, dS, Ckk, T, distances, energies = get_concentration_distance_dependence(C_nomi, ρi, solute_int)
        tidx = find( x -> x==op_temp, T  )
        C[i,:] .= vcat(Ckk[tidx,:]...)
        println("ρ = ", ρi)
        println("C_nom = $(floor(Int64,C_nomi)) = ", C_nomi)
        println("C = ", C[i,:])
    end

    return C, distances
end

function plot_concentration_distance_dependence()
    Tf = 1200
    T = collect(Float64, 0:10:Tf)
    ρ = Float64[ 1.e12, 1.e14, 1.e15, 0.5e16 ]
    C_nom = [ 10., 100., 433., 1000.] ./ 1e6

    # Using the lorentzian for the interaction between solute and dislocation
    solute_int = C_Lorentzian{Float64}()

    for ρi in ρ
        C, distances = collect_concentrations(C_nom, ρi, solute_int)

        Å =  1e-10
        a = 2.87Å
        b = √3/2 * a
        b_mag = 2.87 * √3 / 2

        # Plot the concentrations found at operating temperatures for all sites around the core
        # > Plot it against the energy index

        pyplot( xlims = (0, maximum(distances./b_mag)),
                ylims = (0, 1),
                size=(800,600),
                xticks=collect(0:0.5:2.5),
                yticks=collect(0:0.2:1), legend=false)

        fnt = Plots.font( "Helvetica", 30 )
        default(titlefont = fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)


        xlabel!("Distance [b]")
        ylabel!("Cd")

        # Plot of the distance dependence of the in Burger's vectors

        for (i,C_nomi) in enumerate(C_nom)
            if i == 1
                p = plot( distances./b_mag, C[i,:], label="$(floor(Int64, C_nomi*1e6)) appm", lw=3 )
            elseif i == length(C_nom)
                break
            else
                plot!(distances./b_mag, C[i,:], label="$(floor(Int64, C_nomi*1e6)) appm", lw=3 )
            end
        end
        plot!(distances./b_mag, C[end,:], label="$(floor(Int64, C_nomi*1e6)) appm", lw=3 )

        path="figures"
        cd("figures")

        ρ_string = @sprintf("%.2E", ρi)
        output="concentration_vs_solute_distance_$(ρ_string)_appm"
        png(output)
        cd("..")
    end
end



# plot_concentration_distance_dependence()
end
