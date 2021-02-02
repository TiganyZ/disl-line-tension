module TrapSites

include("trap_sites.jl")
import InteractionTypes: ConcSolutes, C_Lorentzian, H_Lorentzian, lorentzian, dlorentzian

export convert_sitelabel_to_pos_function, get_reference_concentration, get_paths, convert_sitelabel_to_pos_function,  get_interaction_energy

""" This file combines the concentrations derived from the McClean
Isotherm and uses Maxwell-Boltzmann statistics such that one can get the full interaction energy which is modified by the occupancy of the sites. The trap sites defined such that there is a
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



# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Defining the trap sites          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////

function get_paths(core_position)
    Ei_H_paths, H_Ei_paths, H_Ef_paths, Ef_H_paths = trap_site_paths()
    Ei_H_isolated, H_Ei_isolated, H_Ef_isolated, Ef_H_isolated = isolated_trap_sites()

    midpoint = (1/6. * √2 * 2.87 * √3)
    if core_position[1] < midpoint
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
    return hcat(positions...)
end


function get_references(forward, backward)
    return vcat(keys(forward[2])..., keys(backward[2])...,
                forward[3]..., backward[3]...)
end

# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>         Definining concentrations          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////


function concentrations(trap_positions, core_position, conc_func)
    C_tot = zeros(eltype(trap_positions), size(trap_positions,2))
    for i in 1:size(trap_positions,2)
        C_tot[i] += conc_func( norm(core_position - trap_positions[:,i]))
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

    # @show region_forward
    # @show region_backward

    pf = get_proportion(region_forward, core_position)
    pb = get_proportion(region_backward, core_position)

    scaling = Dict{SiteLabel,eltype(core_position)}()
    for (k,v) in trap_path_forward
        scaling[k] = scale_one_to_many_interaction(pf, k, v, trap_path_forward, trap_path_backward)
    end

    for (k,v) in trap_path_backward
        scaling[k] = scale_one_to_many_interaction(pb, k, v, trap_path_backward, trap_path_forward)
    end

    for k in isolated_forward
        scaling[k] = pf
    end

    for k in isolated_backward
        scaling[k] = pb
    end

    # @show references 
    return [scaling[k] for k in references]
end

# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>         Definining references          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////


function get_position_and_scaled_concentration(solutes::ConcSolutes, core_position, scaling, forward, backward)

    positions = get_all_trap_positions(forward, backward, core_position, solutes.convert_sitelabel)
    references = get_references(forward, backward)
    concs =  concentrations(positions, core_position, solutes.conc_func)

    energies = get_interaction_energy_array(solutes, core_position, positions)
    occupancies = partial_occupancies(references, scaling, energies, solutes.T, solutes.ref_conc_sum)

    return positions, occupancies
end

function get_reference_concentration(forward, backward, convert_sitelabel, conc_func)
    core_position = zeros(2)
    positions = get_all_trap_positions(forward, backward, core_position, convert_sitelabel)
    concs = concentrations(positions, core_position, conc_func)
    return sum(concs)
end


# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>         Definining occupancies          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////

function partial_occupancies(references, scaling, energies, T, ref_conc_sum)
    # Remember, the degeneracy factor in the in Maxwell-Boltzmann
    # statistics come for sites which will have the same energy but
    # they are distinguishable by other means. Does this apply here?

    # Perhaps make function which modifies the true concentration, as one can imagine a particle taking a path to the other dislocation smoothly.

    kb = 0.000086173324 # eV/K
    #    e_ref = energies
    Z = sum( exp(-Ei/(kb*T)) for Ei in energies)
    #    C = sum(concs) * 
    #    norm =  ref_conc_sum / C
    return [ ref_conc_sum * scaling[i] * exp(-energies[i]/(kb*T)) / Z  for i in 1:length(energies)]
end

# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>         Interaction energy          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////


function get_single_interaction_energy(solutes::ConcSolutes, core_position, occupancy, position)
    E_int = 0.0
    b_mag = 2.87 * √3 / 2

    dist = norm( core_position - position )

    # Save meV conversion to later 
    #                       eV/b
    E_int += occupancy *  lorentzian(solutes.interaction_type,  dist / b_mag)
    
    return E_int
end


function get_interaction_energy_array(solutes::ConcSolutes, core_position, positions)
    N = size(positions,2)
    energies = zeros(eltype(core_position), N)
    
    for i in 1:N
        p = positions[:,i]
        energies[i] = get_single_interaction_energy(solutes, core_position, 1.0, p)
    end
    return energies
end



function get_interaction_energy(solutes::ConcSolutes, core_position)
    E_int = 0.0

    forward, backward  = get_paths(core_position)
    references = get_references(forward, backward)
    scaling = get_scaling_for_all_sites(core_position, forward, backward, references)
    positions, occupancy = get_position_and_scaled_concentration(solutes, core_position, scaling, forward, backward)
    write_trap_positions_images(positions, occupancy)

    for i in 1:size(positions,2)
        E_int += get_single_interaction_energy(solutes, core_position, occupancy[i], positions[:,i])
    end
    return E_int * 1000 
end


function write_trap_positions_images(positions, occupancy)
    file_ext = "trap_positions_occupancy"
    all_data = vcat(positions,occupancy')'
    
    mode = "a"
    
    open(file_ext, mode) do io 
        println(io,"") 
    end
    
    write_object_to_file(all_data, file_ext, mode)
end


function write_object_to_file(object, filename,  mode)
    open( filename,  mode) do io
        writedlm( io,  object)
    end
end


end
