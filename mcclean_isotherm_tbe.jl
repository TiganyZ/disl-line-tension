using DSP
using Plots
using NLsolve
using Dierckx
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


function get_concentration_data_all(C_nom, ρ, energies, T, V_CC_ventelon)
    all_Cₖ = zeros(Float64, length(T), length(energies), length(ρ), length(C_nom))
    for (cn, C_nomi) in enumerate(C_nom)
        for (ρn, ρi) in enumerate(ρ)
            D = DislSystem( ρi, 1., a, b, C_nomi )
            
            for (en,E⁰) in enumerate(energies)
                for (tn,Ti) in enumerate(T)
                    S = Seg(E⁰, V_CC_ventelon, Ti, kb )
                    # C = solve_nl(D, S, D.C_nom/3/2 ).zero[1]
                    # if tn > 3:
                    #     if abs( C - mean(all_Cₖ[tn-3:tn, en, ρn, cn]) ) 
                    all_Cₖ[tn, en, ρn, cn] = solve_nl(D, S).zero[1]
                end
            end
        end
    end
    return all_Cₖ
end

function clean_concentration_data_all(all_Cₖ, C_nom, ρ, energies, T)
    
    Ckk = copy(all_Cₖ)
    for (cn, C_nomi) in enumerate(C_nom)
        for (ρn, ρi) in enumerate(ρ)
            for (en,E⁰) in enumerate(energies)
                for tn in length(T):-1:1 #(tn,Ti) in enumerate(T)
                    if length(T)-tn >= 2
                        grad   = ( Ckk[tn+1, en, ρn, cn] -  Ckk[tn  , en, ρn, cn]  ) / (T[tn+1]-T[tn])
                        gradp1 = ( Ckk[tn+2, en, ρn, cn] -  Ckk[tn+1, en, ρn, cn]  ) / (T[tn+1]-T[tn])                    
                        if abs(gradp1 - grad) > 0.001
                            Ckk[tn, en, ρn, cn] .= Ckk[tn+1, en, ρn, cn] # ( Ckk[tn+1, en, ρn, cn] + Ckk[tn, en, ρn, cn] ) / 2
                        end
                    end#solve_nl(D, S, D.C_nom/3/2 ).zero[1]
                end
            end
        end
    end
    return Ckk
end


function get_concentration_data(C_nomi, ρi, energies, T, V_CC_ventelon)
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
        for tn in length(T):-1:1 #(tn,Ti) in enumerate(T)
            if length(T)-tn >= 2
                grad   = ( Ckk[tn+1, en] -  Ckk[tn  , en]  ) / (T[tn+1]-T[tn])
                gradp1 = ( Ckk[tn+2, en] -  Ckk[tn+1, en]  ) / (T[tn+1]-T[tn])                    
                if abs(gradp1 - grad) > 0.001
                    Ckk[tn, en] .= Ckk[tn+1, en] # ( Ckk[tn+1, en, ρn, cn] + Ckk[tn, en, ρn, cn] ) / 2
                end
            end#solve_nl(D, S, D.C_nom/3/2 ).zero[1]
        end
    end
    return Ckk
end

####----     Interpolate between the energy points    ----####

function get_concentration_distance_dependence(C_nom, ρ)
    # > We have the binding energy of carbon as a function of distance from the dislocation core.
    # > One can sample from this distribution at multiple points and obtain the concentration as a function of distance at a certain temperature.
    # > Operating temperature is 320K
    include("interaction_types.jl")
    using InteractionTypes

    solute_int = C_Lorentzian{Float64}()
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
    energies = [ lorentzian( solute_int, d/b_mag ) for d in distances ]
    
    all_Cₖ = get_concentration_data(C_nom, ρ, energies, T, V_CC_ventelon)
    Ckk = clean_concentration_data(all_Cₖ, energies, T)

    # Now interpolate the cleaned data at T = 320
    op_temp = 320.
    tidx = find( x -> x==op_temp, T  )

    # We now have the concentrations as a function of distance, as there is a one-to-one correspondence between the energy and the distance
    # > We can interpolate this function as a spline and return the function.
    # > Would make sense to have this in interaction types for clarity and brevity, where we can export the type and define the gradient as such. 

    conc = Ckk[ tidx, : ]
    
    # This is defined where the distances are in angstrom, so no conversion 
    S  = Spline1D(distances, conc, w = ones(length(x)), k = 3, bc = "error")
    S' = x -> derivative(S, x)

    return S, S'
    
    

function get_concentration_distance_dependence()
    # > We have the binding energy of carbon as a function of distance from the dislocation core.
    # > One can sample from this distribution at multiple points and obtain the concentration as a function of distance at a certain temperature.
    # > Operating temperature is 320K
    include("interaction_types.jl")
    using InteractionTypes

    solute_int = C_Lorentzian{Float64}()

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
    ρ = Float64[ 1.e12, 1.e14, 1.e15, 0.5e16 ] 
    C_nom = [ 10., 20, 30, 40, 50, 60, 70, 80, 90, 100., 200, 300, 433., 500, 1000.] ./ 1e6
    
    b_mag = 2.87 * √3 / 2
    # Want range which is up to 2.5 burgers vectors
    distances = collect(0:0.1:2.5) .* b_mag
    energies = [ lorentzian( solute_int, d/b_mag ) for d in distances ]
    
    all_Cₖ = get_concentration_data(C_nom, ρ, energies, T, V_CC_ventelon)
    Ckk = clean_concentration_data(all_Cₖ, C_nom, ρ, energies, T)

    # Now interpolate the cleaned data at T = 320
    op_temp = 320.


    # Once this is done can store array of splines which interpolate between the data and query for particular one in simulation
    
    
    
    
end


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
T = collect(Float64, 0:10:Tf) 
ρ = Float64[ 1.e12, 1.e14, 1.e15, 0.5e16 ] 
C_nom = [ 10., 100., 433., 1000.] ./ 1e6


labels = ["E1 ", "E2 ", "E3 ", "E4 ", "E5 ", "E6 ", "E7 ", "E8 ", "E9 ", "E10",
          "H1 ", "H2 ", "H3 ", "H4 ", "H5 ", "H6 ", "H7 "]

# New values for the interaction energy E⁰, as detailed in the itakura paper.
# > In units of meV 
energies = [
       0.001 *  -775.17     ,
       0.001 *  -792.585    ,
       0.001 *  -139.236    ,
       0.001 *  -234.001    ,
       0.001 *  -791.454    ,
       0.001 *  -603.242    ,
       0.001 *  -388.045    ,
       0.001 *  -177.539    ,
       0.001 *  -683.0      ,                                                                                                           
       0.001 *  -67.0       ,
       0.001 *  -0.881*1000, #1291.3     ,
       0.001 *  -698.331    ,
       0.001 *  -467.286    ,
       0.001 *  -315.768    ,
       0.001 *  -409.328    ,
       0.001 *  114.302     ,
       0.001 *  -819.0     ]



all_Cₖ = zeros(Float64, length(T), length(energies), length(ρ), length(C_nom))


for (cn, C_nomi) in enumerate(C_nom)
    for (ρn, ρi) in enumerate(ρ)
        D = DislSystem( ρi, 1., a, b, C_nomi )

        for (en,E⁰) in enumerate(energies)
            for (tn,Ti) in enumerate(T)
                S = Seg(E⁰, V_CC_ventelon, Ti, kb )
                # C = solve_nl(D, S, D.C_nom/3/2 ).zero[1]
                # if tn > 3:
                #     if abs( C - mean(all_Cₖ[tn-3:tn, en, ρn, cn]) ) 
                all_Cₖ[tn, en, ρn, cn] = solve_nl(D, S).zero[1]
            end
        end
    end
end


###--- Clean the data

function smooth(y, win_len=11, win_method=2)
    # This function requires the DSP package to be installed
    # 1: flat
    # 2: hanning
    # 3: hamming ...
    if win_len%2==0
        win_len+=1 # only use odd numbers
    end
    if win_method == 1
        w=ones(win_len)
    elseif win_method==2
        w=DSP.hanning(win_len)
    elseif win_method==3
        w=DSP.hamming(win_len)
    end
    
    if win_len<3
        return y
    elseif length(y)<win_len
        return y
    else
        y_new = [2*y[1]-flipdim(y[1:win_len],1); y[:]; 2*y[end]-flipdim(y[end-win_len:end],1)]
        y_smooth = conv(y_new, w/sum(w))
        ind = floor(Int, 1.5*win_len)
  return y_smooth[1+ind:end-ind-1]
    end

end # end of function

Ckk = copy(all_Cₖ)
for (cn, C_nomi) in enumerate(C_nom)
    for (ρn, ρi) in enumerate(ρ)
        for (en,E⁰) in enumerate(energies)
            for tn in length(T):-1:1 #(tn,Ti) in enumerate(T)
                if length(T)-tn >= 2
                    grad   = ( Ckk[tn+1, en, ρn, cn] -  Ckk[tn  , en, ρn, cn]  ) / (T[tn+1]-T[tn])
                    gradp1 = ( Ckk[tn+2, en, ρn, cn] -  Ckk[tn+1, en, ρn, cn]  ) / (T[tn+1]-T[tn])                    
                    if abs(gradp1 - grad) > 0.001
                        Ckk[tn, en, ρn, cn] .= Ckk[tn+1, en, ρn, cn] # ( Ckk[tn+1, en, ρn, cn] + Ckk[tn, en, ρn, cn] ) / 2
                    end
                end#solve_nl(D, S, D.C_nom/3/2 ).zero[1]
            end
        end
    end
end




# C_nom = [ 10., 100., 433., 1000.] ./ 1e6
pyplot( xlims = (0, Tf),
        ylims = (0, 1),
        size=(800,600),
        xticks=collect(0:200:Tf),
        yticks=collect(0:0.2:1.2), legend=false)

fnt = Plots.font( "Helvetica", 30 )
default(titlefont = fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

colors = palette(:tab10)
ls = [:solid, :dash, :dot, :dashdot]

easy_plots = []
###---   For the easy core
for k in 1:length(C_nom)
    for j in 1:length(ρ)
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
        if j == length(ρ)
            title!( @sprintf("%4d appm", floor(Int64, C_nom[k]*1e6)) )
            for jj in 1:length(ρ)
                plot!(T, Ckk[:,1,jj,k], linestyle=ls[jj], linewidth=3, show=false, label=@sprintf("ρ = %1.1e", ρ[jj]), linecolor=colors[1])
            end
            break
        end
    end #; plot!(T, Ckk[:,end,end,k], linestyle=ls[end], linewidth=3, show=true, label="", linecolor=colors[10], legend=false);
    output="mcclean_isotherm_rep_easy_$(floor(Int64, C_nom[k]*1e6))_appm_nl_nt"
    png(output)
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
for k in 1:length(C_nom)
    for j in 1:length(ρ)
        for i in 1:length(energies[11:end])
            if i == j == 1
                pp = plot(T, Ckk[:,10+i,j,k], label = labels[10+i], linestyle=ls[j], linewidth=3, show=false, linecolor=colors[i])
            xlabel!("T[K]")
                ylabel!("Cd")
            elseif j == 1
                plot!(T, Ckk[:,10+i,j,k], label = labels[10+i], linestyle=ls[j], linewidth=3, show=false, linecolor=colors[i])
            else
                plot!(T, Ckk[:,10+i,j,k], linestyle=ls[j], linewidth=3, show=false, label="", linecolor=colors[i])
            end
        end
        if j == length(ρ)
            title!( " ")
            for jj in 1:length(ρ)
                plot!(T, Ckk[:,11,jj,k], linestyle=ls[jj], linewidth=3, show=false, label=@sprintf("ρ = %1.1e", ρ[jj]), linecolor=colors[1])
            end
            break
        end
    end
    output="mcclean_isotherm_rep_hard_$(floor(Int64, C_nom[k]*1e6))_appm_nl_nt"
    png(output)
end







    for j in 1:length(ρ)
        for i in 1:length(energies[11:end])
            if i == j == 1
                pp = plot( label = labels[10+i], linestyle=ls[j], linewidth=3, show=false, linecolor=colors[i])
                xlabel!("T[K]")
                ylabel!("Cd")
            elseif j == 1
                plot!(label = labels[10+i], linestyle=ls[j], linewidth=3, show=false, linecolor=colors[i])
            else
                plot!( linestyle=ls[j], linewidth=3, show=false, label="", linecolor=colors[i])
            end
        end
        if j == length(ρ)
            title!( " ")
            for jj in 1:length(ρ)
                plot!(linestyle=ls[jj], linewidth=3, show=false, label=@sprintf("ρ = %1.1e", ρ[jj]), linecolor=colors[1])
            end
            break
        end
    end
# pyplot( xlims = (0, Tf),
#         ylims = (0, 1),
#         size=(800,600),
#         xticks=collect(0:100:Tf),
#         yticks=collect(0:0.2:1))

# fnt = Plots.font( "Helvetica", 16 )                                                                                                 
# default(titlefont = fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)


# #pyplot()
# counter=0
# p = Plots.Plot{Plots.PyPlotBackend}[]
# ls = [:solid, :dash, :dot, :dashdot]

# for k in 1:length(C_nom)
#     for j in 1:length(ρ)
#         for i in 1:length(energies)
#             if i == j == 1
#                 pp = plot(T, all_Cₖ[:,i,j,k], label = labels[i], linestyle=ls[j], linewidth=2, show=false)
#             elseif j == 1
#                 plot!(T, all_Cₖ[:,i,j,k], label = labels[i], linestyle=ls[j], linewidth=2, show=false)
#             elseif i == length(energies) && j == length(ρ)
#                 plot!(T, all_Cₖ[:,i,j,k], linestyle=ls[j], linewidth=2, show=true)
#             else
#                 plot!(T, all_Cₖ[:,i,j,k], linestyle=ls[j], linewidth=2, show=false)
#             end
#         end
#     end
#     title!( @sprintf("%4d appm", floor(Int64, C_nom[k]*1e6)) )
#     append!(p, pp)
# end
# p

# for i in 1:length(energies)
#     if i == j == 1
#         pp = plot(T, all_Cₖ[:,i,j,k], label = labels[i], linestyle=ls[j], linewidth=2, show=false)
#     elseif j == 1
#         plot!(T, all_Cₖ[:,i,j,k], label = labels[i], linestyle=ls[j], linewidth=2, show=false)
#     elseif i == length(energies) && j == length(ρ)
#         plot!(T, all_Cₖ[:,i,j,k], linestyle=ls[j], linewidth=2, show=true)
#     else
#         plot!(T, all_Cₖ[:,i,j,k], linestyle=ls[j], linewidth=2, show=false)
#     end
# end


# open("all_concentration_data_for_sites.dat", "w") do io
#     write(io, "# Site   T[K]    Ck     dislocation_density    C_nominal    ---###   Contains data for all sites" ) 
#     for k in 1:length(C_nom)
#         for j in 1:length(ρ)
#             for i in 1:length(energies)
#                 for h in 1:length(T)
#                     write(io, "$(label[i]) $(T[h]) $(all_Cₖ[i,h,j,k]) $(ρ[j]) $(C_nom[k]) ")
#                 end
#             end
#         end
#     end
# end

# j = 1
# for i in 1:length(energies)
#         if i == j == 1
#             pp = plot(T, Ckk[:,i,j,k], label = labels[i], linestyle=ls[j], linewidth=2, show=false)
#         elseif j == 1
#             plot!(T, Ckk[:,i,j,k], label = labels[i], linestyle=ls[j], linewidth=2, show=false)
#                elseif i == length(energies) && j == length(ρ)
#             plot!(T, Ckk[:,i,j,k], linestyle=ls[j], linewidth=2, show=true)
#         else
#             plot!(T, Ckk[:,i,j,k], linestyle=ls[j], linewidth=2, show=false)
#         end
#     end
# end




