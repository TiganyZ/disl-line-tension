include(trap_sites.jl)
include(mcclean_isotherm_conc_dist.jl)

using TrapSites
using McCleanIsotherm

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

 """



# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Defining the actual interaction space          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////


function get_active_sites(site_to_site, site_to_positions, core_position)
end


# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Scaling during core motion          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////


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
