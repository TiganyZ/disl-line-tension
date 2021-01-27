module TrapSites

include(trap_sites.jl)
include(interaction_types.jl)
include(mcclean_isotherm_conc_dist.jl)
using McCleanIsotherm
using InteractionTypes


""" This file combines the concentrations derived from the McClean
Isotherm and uses the trap sites defined such that there is a
consistent definition of total concentration across all dislocation
core positions, which do not change even when there is a different
number of interactions (sites) which interact with the core at any
given time


> Need to constrain total concentration. Can either use the easy core or hard core for reference configuration.

> We run into a problem with a different number of sites around a dislocation,

  || The /effective/ concentration of a particular site /decreases/, if the number of sites around the new core /increases/. ||

> How do we rectify this?

> We also have the problem of a many-to-one mapping of sites in going from the easy to hard core...

> Perhaps for the many-to-one mapping, we see how many are in this mapping and then divide by the number based on the distance, just like we moderate the trap site position itself based on distance to the core.

> Count keys which go to one particular position, then use this to scale the concentration by the number of sites mapped to that position.

> How will the trap sites move as the core moves?
>> Initially, we can change the trap site positions based on where the core is along the line from the easy core to the hard core.
>> Based on the distance between the core position and the Ei and H positions, we can moderate the distance according to

    d_Ei_H = Distance between the cores
    d_H    = Distance from current core position to hard core
    d_Ei   = Distance from current core position to easy core


    Can project each of these lines and distances on the d_Ei_H line
    From this we can get an effective distance between each of the cores


Process:
        > Defined sites at a given core position will have their concentration defined by the distance dependence from
        >>  get_concentration_distance_dependence_splines(C_nom, ρ, solute_int) from McClean Isotherm

        > Sum over all concentrations of sites which surround the core.

        > Use this sum to scale with different sites.
        >>



    > Need to be clear with what's happening here
    > We get the trap sites based on how far they are from the two end points which define a region (an easy and hard core position)
    > We define the backwards paths the same way
    > Each of these define one position "decaying" to another, as such when there are many positions going to one we must scale
    >
 """

# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Defining the types          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////

abstract type AbstractRegion end

struct Ei_H <: AbstractRegion end
struct Ef_H <: AbstractRegion end
struct H_Ei <: AbstractRegion end
struct H_Ef <: AbstractRegion end


struct ConcSolutes{T <: AbstractFloat}
    interact::Bool
    interaction_type::Union{InteractionTypes.C_Lorentzian{T},InteractionTypes.H_Lorentzian{T}}
    convert_sitelabel::Function
end

# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Defining the trap sites          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////

function get_paths(core_position)
    Ei_H_paths, H_Ei_paths, H_Ef_paths, Ef_H_paths = trap_site_paths()
    Ei_H_isolated, H_Ei_isolated, H_Ef_isolated, Ef_H_isolated = isolated_trap_sites()

    midpoint = (1/6. * √2 * 2.87 * √3)
    if core_position < midpoint
        forward  = Ei_H(), Ei_H_paths, Ei_H_isolated
        backward = H_Ei(), H_Ei_paths, H_Ei_isolated
    else
        forward  = H_Ef(), H_Ef_paths, H_Ef_isolated
        backward = Ef_H(), Ef_H_paths, Ef_H_isolated
    end
    return forward, backward
end



# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Defining the trap sites          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////


function get_trap_site_position(p, convert_sitelabel, initial_sitelabel, final_sitelabel)
    initial_position = convert_sitelabel(initial_sitelabel)
    final_position   = convert_sitelabel(final_sitelabel)

    return (initial_position*p + final_position*(1-p))
end

function push_trap_site!(positions, convert_sitelabel, trap_paths, p)
    for (k,v) in trap_paths
        push!(positions, get_trap_site_position(p,convert_sitelabel,k,v))
    end
end

function get_all_trap_positions(forward, backward, core_position, convert_sitelabel)
    region_forward, trap_path_forward, isolated_forward = forward
    region_backward, trap_path_backward, isolated_backward = backward

    positions = []
    push_trap_site!(positions, convert_sitelabel, trap_path_forward, get_proportion(region_forward, core_position))
    push_trap_site!(positions, convert_sitelabel, trap_path_backward, get_proportion(region_backward, core_position))

    for isolated in [isolated_forward, isolated_backward]
        for sitelabel in isolated
            push!(positions, convert_sitelabel(sitelabel))
        end
    end
    references = vcat(keys(trap_path_forward)...,keys(trap_path_backward)...,
                      keys(isolated_forward)..., keys(isolated_backward)...)
    return hcat(positions...)', references
end



# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>         Definining concentrations          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////


function concentrations(trap_positions, trap_references, scaling, core_position, conc_func)
    C_tot = zeros(size(trap_positions,1))
    for i in 1:size(trap_positions,1)
        C_tot[i] += conc_func( norm(core_position - trap_positions[i,:]))
    end
    return C_tot
end

function dinv_conc_sum(trap_positions, core_positions, conc_func, dconc_func)
    # f = 1/g = 1 / ∑ Cᵢ
    # df = - g' / g² = - g' * f²
    g  = sum(concentrations(trap_positions, core_position, conc_func))
    dg = sum(concentrations(trap_positions, core_position, dconc_func))
    return - dg / g^2
end


get_proportion(region::Ei_H, core_position) = (1. - core_position[1] / (1/6. * √2 * 2.87 * √3))
get_proportion(region::H_Ei, core_position) =       core_position[1] / (1/6. * √2 * 2.87 * √3)
get_proportion(region::H_Ef, core_position) = (1. - (core_position[1] - 1/6. * √2 * 2.87 * √3) / (1/6. * √2 * 2.87 * √3))
get_proportion(region::Ef_H, core_position) = (core_position[1] - 1/6. * √2 * 2.87 * √3) / (1/6. * √2 * 2.87 * √3)

get_dproportion(region::Union{Ei_H,H_Ef}, core_position, direction) = direction == 1 ? - 1.0 / (1/6. * √2 * 2.87 * √3) : 0.0
get_dproportion(region::Union{H_Ei,Ef_H}, core_position, direction) = direction == 1 ?   1.0 / (1/6. * √2 * 2.87 * √3) : 0.0


# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Scaling during core motion          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////

function scale_one_to_many_interaction(p, initial_sitelabel, final_sitelabel, core_path_forward, core_path_backward)
    # Find all labels in the interaction

    initial_duplicates =  sum([k == initial_sitelabel for (k,v) in core_path_forward  ])
    initial_duplicates += sum([v == initial_sitelabel for (k,v) in core_path_backward ])

    final_duplicates =  sum([v == final_sitelabel for (k,v) in core_path_forward  ])
    final_duplicates += sum([k == final_sitelabel for (k,v) in core_path_backward ])

    # Want to scale from 1 to 1/n_duplicates
    return (p*1./initial_duplicates + (1. - p)/final_duplicates)

end

function get_scaling_for_all_sites(core_position, forward, backward, references)
    region_forward, trap_path_forward, isolated_forward = forward
    region_backward, trap_path_backward, isolated_backward = backward

    p = get_proportion(region_forward, core_position)

    scaling = Dict{SiteLabel,Float64}()
    for (k,v) in trap_path_forward
        scaling[k] = scale_one_to_many_interaction(p, k, v, trap_path_forward, trap_path_backward)
    end

    for (k,v) in trap_path_backward
        scaling[k] = scale_one_to_many_interaction(p, k, v, trap_path_backward, trap_path_forward)
    end


    for (k,v) in isolated_forward
        scaling[k] = get_proportion(region_forward, core_position)
    end

    for (k,v) in isolated_backward
        scaling[k] = get_proportion(region_backward, core_position)
    end

    scale_array = [ scaling[k] for k in references ]


    return scale_array
end

function convert_conc_to_partial_occupancies!(concs, energies, T, ref_conc_sum)
    # Remember, the degeneracy factor in the in Maxwell-Boltzmann
    # statistics come for sites which will have the same energy but
    # they are distinguishable by other means. Does this apply here?

    # The whole scaling malarkey isn't useful to us at all then: the partial occupancies will sort it out for us and give the correct statistics?

    kb = 0.000086173324 # eV/K
    T = 320.0 # K
    Z = sum( exp(-Ei/(kb*T)) for Ei in energies)
    return [ concs[i] * exp(-Ei/(kb*T)) / Z / ref_conc_sum for Ei in energies]
end


# > Note: Scaling by projection won't work so well. It is better that
# > is just measures the x coordinate as the dislocation traverses
# > from one Peierls valley to the next.

# > Perhaps create Isolated type, for those sites that scale to zero, and Normal, for those that interact normally

function get_single_interaction_energy(solutes::ConcSolutes, Pⱼ, scale, position, derivative=false)
    E_int = 0.0
    b_mag = 2.87 * √3 / 2

    dv = Pⱼ[1:2] - position
    dist = norm( dv )


    if derivative
        dd = get_dd(direction)
                #         meV       eV/b
        E_int += scale * 1000. * dlorentzian(solutes.interaction_type,  dist / b_mag)*dd( dv[1]/b_mag, dv[2]/b_mag )
    else
        E_int += scale * 1000. *  lorentzian(solutes.interaction_type,  dist / b_mag)
    end

    return E_int
end


function get_interaction_energy_array(solutes::ConcSolutes, Pⱼ, scaling, positions, derivative=false)
    return [ get_single_interaction_energy(solutes, Pⱼ, s, p, derivative) for (s,p) in zip(scaling,positions)]
end

function get_total_conc_interaction_energy(solutes::ConcSolutes, core_position, derivative=false)
    E_int = 0.0

    forward, backward  = get_paths(core_position)
    positions, concs = get_position_and_scaled_concentration(core_position, forward, backward, convert_sitelabel, conc_func, ref_conc_sum)

    # What we might actually want to do here, is scale the energy by
    # the scaling factor. This will mean the effective interaction is
    # more reasonable, and it enforces a more realistic occupancy

    for i in 1:size(positions,1)
        E_int += concs[i] * get_single_interaction_energy(solutes, Pⱼ, concs[1,:], positions[i,:], derivative=derivative)
    end
    return E_int
end


function get_position_and_scaled_concentration(core_position, forward, backward, convert_sitelabel, conc_func, ref_conc_sum)

    positions, references = get_all_trap_positions(forward, backward, core_position, convert_sitelabel)
    concs = concentrations(positions, core_position, conc_func)
    #    scaling = get_scaling_for_all_sites(core_position, forward, backward, references)

    energy_array = get_interaction_energy_array(solutes::ConcSolutes, Pⱼ, scaling, positions, derivative=false)
    convert_conc_to_partial_occupancies!(concs, energies, T, ref_conc_sum)
    # Using the easy core as the reference
    # scale_concentrations!(concs, scaling, references, ref_conc_sum)

    # Now we have scaled the concentrations we work out the interaction energy
    # Use the lorentzian
    # As per Ivo's email, this should just be the scaled concentrations multiplied by the interaction energy of that site
    return positions, concs
end

function get_reference_concentration(forward, backward, convert_sitelabel, conc_func)
    core_position = zeros(2)
    positions, references = get_all_trap_positions(forward, backward, core_position, convert_sitelabel)
    concentrations = concentrations(positions, core_position, conc_func)
    return sum(concentrations)
end


function interaction_with_proportions(f, df, conc_f, dconc_f, region, initial_sitelabel, final_sitelabel, convert_sitelabel)
    overall = get_overall_factor(region, core_position)
    p       = get_proportion(region, core_position)
    dp      = get_dproportion(region, core_position)

    trap_position = get_trap_site_position(p, initial_sitelabel, final_sitelabel)

    conc, dconc = get_concentration(trap_position, core_position, conc_f, dconc_f)

    # Function will just be the interaction energy function, as it has been worked out
    factor = scale_one_to_many_interaction(p, final_sitelabel, core_paths)

    # Return the function and derivative
    # Factor is constant so only multiplicative factor which depends on x is the concentration
    # Can just keep it in the energy, as we have already resolved it in the interaction in dislocation types
    # Or actually explicitly pyt it in, need integration with Solute type.
    # h = f*(g*p)
    # h' = f*(g*p)' + f'*(g*p)
    # h' = f*(g'*p + g*p') + f'*(g*p)

    fres  =  f( norm(core_position - site_position) )
    dfres = df( norm(core_position - site_position) )

    ret = fres * (dconc*p + dp*conc) + dfres*(conc*p)
    return fres
end

function get_interaction_energy(solutes::Solutes, Pⱼ, derivative=false)
    E_int = 0.0
    b_mag = 2.87 * √3 / 2

    for idx in 1:size(solutes.positions,2)

        if abs.( j*b_mag .- solutes.positions[3,idx] ) .< 1e-3

            dv = Pⱼ[1:2] .- solutes.positions[1:2,idx]
            dist = norm( dv  )

            if derivative
                dd = get_dd(direction)
                #         meV       eV/b
                E_int += 1000. * dlorentzian(solutes.interaction_type,  dist / b_mag)*dd( dv[1]/b_mag, dv[2]/b_mag )
            else
                E_int += 1000. *   lorentzian(solutes.interaction_type,  dist / b_mag)
            end

        end
    end

    return E_int
end

function get_dd(direction)
    if direction == 1
        dd = (x,y) -> x/√(x^2+y^2)
    else
        dd = (x,y) -> y/√(x^2+y^2)
    end
    return dd
end



function initialise_interactions()
    convert_sitelabel_to_pos = convert_sitelabel_to_pos_function()
end





function proportion_upon_projection_line(core, core_position, side)
    alat = √2 * 2.87
    lengths = [ √(3) * alat, alat]

    if core == "Ei"
        reference_core_position = zeros(2)
        diff = core_position - reference_core_position
    elseif core == "Ef"
        reference_core_position = [ 2/6.*lengths[1], 0.]
        diff = core_position - reference_core_position
    elseif core == "H"
        reference_core_position = [ 1/6.*lengths[1], 1/6.*lengths[2] ]
        diff = reference_core_position - core_position
    end

    if side == "left"
        E_H_diff = [ 1/6.*lengths[1],  1/6.*lengths[2] ]
    else
        E_H_diff = [ 1/6.*lengths[1], -1/6.*lengths[2] ]
    end

    proportion = (diff * E_H_diff') / norm(E_H_diff)

    return proportion
end



end
