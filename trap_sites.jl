module TrapSites

import Base
using Plots

""" This file describes the trap site occupation and the
redistribution of the trap site probabilities

The spline function which describes the concentration dependence comes
from McCleanIsotherm. This function will be used to define the solute
concentration in each dislocation segment along the line.

We need to ascertain all the sites along the line.

> Use McClean Isotherm self-consistent solution such that we can
> find the probability of occupancy of a particular trap site by
  χᵢ = Cₖ, from the self-consistent McClean Isotherm

> The /total/ carbon occupancy is constant, i.e. ∑χᵢ = const.
> This means that carbon will redistribute among the traps within the plane perpendicular to the dislocation.

> For slow glide, and the maintenence of equilibrium, then as
> the dislocation moves between the two Peierls valleys we can
> define an occupation probability
   χᵉᵢ(x) = χₜ exp{ -Eᵢ(x)/kT } / ∑ᵢexp{ -Eᵢ(x)/kT }
> Where Eᵢ(x) is the interaction energy given by a lorentzian function among all the trap sites.

> Can assume that the nominal carbon concentration is at 433 appm, the solubility limit of 0.02wt% ferrite.

> Will have the easy and hard core have different distribution of sites.
> Immediately we run into a problem as around each core there is a different number of sites allowed around each one.

> Criticism: This is not realistic as the sites around the hard core
  would actually be in a different position if the core was allowed to
  relax


    Action of getting sites are rotating by 0, 2π/3, and 4π/3 then a reflection of all of them
    Sites are partitioned into 6 sectors, which are different for easy and hard core
    >> Easy <<
      4  1
    2      5
      6  3

    >> Hard <<
      3  6
    5      2
      1  4

    Tuples control the type of trap site and the sector in which it lies.

    One can get the positions from the dictionaries

 """


struct SiteLabel
    site::String
    section::Int64
end

struct Isolated
end

struct Normal
end

struct Ei_H
    type::Union{Normal,Isolated}
end

struct H_Ei
    type::Union{Normal,Isolated}
end

struct H_Ef
    type::Union{Normal,Isolated}
end

struct Ef_H
    type::Union{Normal,Isolated}
end


# Note: Even though this type is apparently immutable, there is a problem with indexing the dictionaries with the keys as these types.
#     > This is apparently only a problem with mutable structs, yet this isn't...

# >> Workaround by overloading Base.hash() and ==()

import Base.==, Base.convert, Base.hash

==(a::SiteLabel, b::SiteLabel) = a.site == b.site && a.section == b.section
Base.hash(a::SiteLabel, h::UInt) = hash(a.site, hash(a.section, h))


function Base.convert(L::SiteLabel, D::Dict{SiteLabel, Array{Float64,1}})
    # Convert SiteLabel to position
    return D[L]
end



# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Defining site positions for cores          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////

function create_lattice()

    abcc = 2.87 # * 1.88971616463207
    alat = √2 * abcc

    # plat = reshape([17.32050807568877, 0.0, 0.0 ,0.0 ,16.0 ,0.0 ,0.0 ,0.0, 0.6123724356957945] * alat, (3,3) )

    # For radial configuration
    q = √(3./8.) # This should work for bcc...

    lengths = [ √(3)  1.]

    plat = [[√(3)  0];
            [0.    1] ]

    bcc_unit_cell = [[  0.0   0.   ];
                     [  1/3   0.   ];
                     [  2/3.  0.   ];
                     [  0.5   0.5  ];
                     [  5/6.  0.5  ];
                     [  1/6.  0.5  ]] .* lengths


    # Build lattice of the unit cell to plot

    nxyz = (2, 2)


    lattice = zeros( size(bcc_unit_cell)[1]*( (2*nxyz[1]+1)*(2*nxyz[2]+1) ), 2 )

    xlims = (-4*√3, 4*√3)
    ylims = (-5, 5)

    xlims = (-Inf, Inf)
    ylims = (-Inf, Inf)


    luc = size(bcc_unit_cell)[1]
    n = 0
    for j in -nxyz[1]:nxyz[1]
        for k in -nxyz[2]:nxyz[2]
            for i in 1:size(bcc_unit_cell)[1]

                pos = ( bcc_unit_cell[i,:] + plat[1,:]*j + plat[2,:]*k ) * alat

                if xlims[1] < pos[1]
                    if xlims[2] > pos[1]
                        if ylims[1] < pos[2]
                            if ylims[2] > pos[2]
                                n += 1
                                lattice[n,:] .=  pos
                                println("$n:  position = $bcc_unit_cell[i,:] -> ", pos/alat)
                            end
                        end
                    end
                end


            end
        end
    end

    # centre_vector = 0.5*( plat[1,:]*(Int64(round(nxyz[1]/2))-1) + plat[2,:]*(Int64(round(nxyz[2]/2))-1) )
    # Make everything centred on the initial easy core position

    Ei_core_position = [ 1/6.* √(3) * alat, 1/6.*alat ]

    return lattice .- Ei_core_position'
end


function get_octahedral_position_dict(oct_sites, core)
    abcc = 2.87 # * 1.88971616463207
    alat = √2 * abcc

    # plat = reshape([17.32050807568877, 0.0, 0.0 ,0.0 ,16.0 ,0.0 ,0.0 ,0.0, 0.6123724356957945] * alat, (3,3) )

    # For radial configuration
    q = √(3./8.) # This should work for bcc...

    lengths = [ √(3) * alat, alat, q  * alat]

    if core == "Ei"
        core_position = zeros(2)
        labels = ["Ei1","Ei2","Ei3", "Ei4", "Ei5", "Ei6", "Ei7", "Ei8", "Ei9", "Ei10"]
    elseif core == "Ef"
        core_position = [ 2/6.*lengths[1], 0.]
        labels = ["Ef1","Ef2","Ef3", "Ef4", "Ef5", "Ef6", "Ef7", "Ef8", "Ef9", "Ef10"]
    elseif core == "H"
        core_position = [ 1/6.*lengths[1], 1/6.*lengths[2] ]
        labels = ["H1", "H2",     "H3", "H4", "H5", "H6", "H7"]
    end

    labels = [ SiteLabel(l,i) for i in 1:6 for l in labels]

    site_pos_dict = Dict{SiteLabel, Array{Float64,1}}( labels[i] => oct_sites[i,:] + core_position for i in 1:length(labels))
    return site_pos_dict
end

function get_octahedral_positions(core_type)
    abcc = 2.87 # * 1.88971616463207
    alat = √2 * abcc

    # plat = reshape([17.32050807568877, 0.0, 0.0 ,0.0 ,16.0 ,0.0 ,0.0 ,0.0, 0.6123724356957945] * alat, (3,3) )

    # For radial configuration
    q = √(3./8.) # This should work for bcc...

    lengths = [ √(3) * alat, alat, q  * alat]

    if  core_type == "easy"
        core_position = [36.3218504,  33.14680888]
        relative_core_position = [ 1/6.*lengths[1], 1/6.*lengths[2] ]
        octahedral_positions =[[37.54770758  33.85046096];
                                [37.50566447  34.41147207];
                                [36.92444498  35.52976750];
                                [38.12738250  35.54441951];
                                [36.32784633  36.51662541];
                                [38.68239733  36.53429325];
                                [38.09544321  37.50275944];
                                [36.92344511  37.51493620];
                                [39.28127528  35.50886674];
                                [40.51130213  35.56288039]]   .- core_position'

    elseif core_type == "hard"
        core_position = [37.493523, 33.82327437]
        relative_core_position = [ 2/6.*lengths[1], 2/6.*lengths[2] ]
        octahedral_positions = [[  37.49967197  33.82188364];
                                 [  35.73715334  31.44012169];
                                 [  34.57273597  29.41591641];
                                 [  35.75702326  29.40559471];
                                 [  36.89155022  29.38012684];
                                 [  37.49563378  30.34278051];
                                 [  33.40750846  31.46033380]] .- core_position'

    end

    return octahedral_positions, core_position, relative_core_position
end



function get_all_oct_sites(oct_sites, centre)
    N = size(oct_sites)[1]

    angles = [0., 2π/3., 4π/3]

    sites = zeros( length(angles)*N, 2 )

    for a in eachindex(angles)
        rotation = rotation_matrix(angles[a])

        println("\n Rotation matrix = ")
        @show rotation
        for i in 1:N
            index = i + (a-1)*N
            sites[index,:] = rotate_about_centre( rotation, centre, oct_sites[i,:] )
        end
    end

    return sites
end

function get_all_sites_for_core(core_type)
    oct_sites_initial, core_position, relative_core_position = get_octahedral_positions(core_type)

    ###---   Add also a REFLECTION as well as a rotation
    oct_sites  = get_all_oct_sites(oct_sites_initial, [0,0.])

    # reflect this oct,  -0.5*[√(2/3)*2.87 0.0]
    reflected_oct = zeros(size(oct_sites))
    reflected_oct = reflect_data_about_line!(1, oct_sites, reflected_oct, [0., 1.] )

    oct_sites  = vcat( oct_sites, reflected_oct)
    return oct_sites
end


function create_trap_labels_to_positions_dict()
    Ei_dict = get_octahedral_position_dict(get_all_sites_for_core("easy"), "Ei")
    Ef_dict = get_octahedral_position_dict(get_all_sites_for_core("easy"), "Ef")
    H_dict  = get_octahedral_position_dict(get_all_sites_for_core("hard"), "H")

    return merge(Ei_dict, Ef_dict, H_dict)
end


# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Transformations to get site positions          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////


function rotation_matrix(angle)
    rotation = [[cos(angle) -sin(angle)];
                [sin(angle)  cos(angle)] ]
    return rotation
end


function rotate_about_centre(rotation, centre, position)
    println("\n Rotating about $centre, position = $position")
    final_pos = centre + rotation * ( position - centre )
    println("  >                final position = $final_pos")
    return centre + rotation * ( position - centre )
end


function reflect_data_about_line!(index::Int64, data::Array{Float64,2}, all_data::Array{Float64,2}, line::Array{Float64,1}, translation = [0.,0.]; no_idx = [])
    trans = translation
    for i in range( 1, size(data)[1] )
        if i in no_idx
            continue
        else
            V = data[i,:] - trans
            rv = reflect_2D( V, line )

            all_data[index, : ] = rv + translation
            index += 1
        end
    end
    return all_data
end


function reflect_2D( v::Array{Float64, 1}, l::Array{Float64, 1} )
    # v is vector being reflected
    # l is line in reflection plane
    ref_v = 2 * ( v' * l / (l' * l) )*l - v
    return ref_v

end


# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Definition of the trap site paths          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////

function reflect_site_label(sitelabel)
    if sitelabel.section <= 3
        new_sitelabel = SiteLabel(sitelabel.site, sitelabel.section + 3)
    else
        new_sitelabel = SiteLabel(sitelabel.site, sitelabel.section - 3)
    end
    return new_sitelabel
end

function change_to_alternate_core(sitelabel)
    if contains(sitelabel.site, "i")
        new_site = replace(sitelabel.site, "i" => "f")
    elseif contains(sitelabel.site, "f")
        new_site = replace(sitelabel.site, "f" => "i")
    else
        new_site = sitelabel.site
    end
    return SiteLabel(new_site, sitelabel.section)
end


function change_mapping_to_other_core(sitelabel)
    return reflect_site_label( change_to_alternate_core(sitelabel) )
end

function convert_site_site_dict_to_other_core(dict)
    return Dict( change_mapping_to_other_core(k) => change_mapping_to_other_core(v) for (k,v) in dict  )
end

function isolated_trap_sites()
    Ei_H_isolated = [SiteLabel("Ei10", 2), SiteLabel("Ei10", 3),
                   SiteLabel("Ei7" , 2), SiteLabel("Ei8",  2),  SiteLabel("Ei7", 6), SiteLabel("Ei8", 6)]

    H_Ei_isolated  = [SiteLabel("H7", 3), SiteLabel("H3", 3), SiteLabel("H3", 6),
                      SiteLabel("H3", 6), SiteLabel("H4", 6), SiteLabel("H5", 6),
                      SiteLabel("H3", 2), SiteLabel("H4", 2), SiteLabel("H5", 2),
                      SiteLabel("H7", 4), SiteLabel("H3", 4)                      ]


    Ef_H_isolated = map( change_mapping_to_other_core, Ei_H_isolated )

    H_Ef_isolated = map( change_mapping_to_other_core, H_Ei_isolated )

    return  Ei_H_isolated, H_Ei_isolated, H_Ef_isolated, Ef_H_isolated
end

function trap_site_paths()
    """ These are the paths that I have chosen, not Ivo.

    > I think this is a better solution as he only mapped energies,
    but this is actually trap sites. And it seems more realistic as
    sites which stay in the same position due to being "locked" into
    the Cottrell atmosphere will have a consistent description within
    the model, as the binding energy will tend to that in found in the
    Eshelby limit of using the elastic dipole tensor.

    I can recreate his results or do something which is a bit
    more intuitive, like for the reconstruction of the core.

     > Seems like I will need to distinguish between each of the sites:
     >> e.g. Right tri normal/reflected
     >> Split space into 6
     >>

    >> Each of the actual cores, Hard and Easy have unambiguous trap sites.
    >> For each of the sites we map to the sector also

    > There are some self-mapped sites which can be left to decay to zero
    > But maybe we don't do this for consistency...
    > There is a problem with some residual carbon on the other side of the dislocation as it moves but it should only give a minor effect to the enthalpies
    > I think there could be a larger error in ignoring sites.
    > What evidence is there to support my claim compared to Ivo's
    >> Certainly the E6 site should not decay to zero in some circumstances given the large initial binding energy.

    > Can isolate some atoms to genuinely decay to zero as we move dislocation to hard core position or vice versa
"""
    # >> Easy << >> Hard <<
    #   4  1       3  6
    # 2      5   5      2
    #   6  3       1  4

    Ei_H_paths = Dict{SiteLabel, SiteLabel}( SiteLabel("Ei1",  1) => SiteLabel("H1", 6), # Ei -> H sector 1
                                             SiteLabel("Ei2",  1) => SiteLabel("H1", 6),
                                             SiteLabel("Ei3",  1) => SiteLabel("H1", 6),
                                             SiteLabel("Ei4",  1) => SiteLabel("H1", 6),
                                             SiteLabel("Ei5",  1) => SiteLabel("H2", 3),
                                             SiteLabel("Ei6",  1) => SiteLabel("H2", 6),
                                             SiteLabel("Ei7",  1) => SiteLabel("H2", 6),
                                             SiteLabel("Ei8",  1) => SiteLabel("H2", 3),
                                             SiteLabel("Ei9",  1) => SiteLabel("H2", 6),
                                             SiteLabel("Ei10", 1) => SiteLabel("H6", 6),

                                             SiteLabel("Ei1",  2) => SiteLabel("H2", 5), # Ei -> H sector 2
                                             SiteLabel("Ei2",  2) => SiteLabel("H2", 5),
                                             SiteLabel("Ei3",  2) => SiteLabel("H2", 5),
                                             SiteLabel("Ei4",  2) => SiteLabel("H2", 5),
                                             SiteLabel("Ei5",  2) => SiteLabel("H7", 1),
                                             SiteLabel("Ei6",  2) => SiteLabel("H3", 5),
                                             SiteLabel("Ei7",  2) => SiteLabel("Ei7", 2),
                                             SiteLabel("Ei8",  2) => SiteLabel("Ei8", 2),
                                             SiteLabel("Ei9",  2) => SiteLabel("H4", 5),
                                             SiteLabel("Ei10", 2) => SiteLabel("Ei10", 2), #SiteLabel("Ei10", 2) => SiteLabel("H5", 5),

                                             SiteLabel("Ei1",  3) => SiteLabel("H2", 1), # Ei -> H sector 3
                                             SiteLabel("Ei2",  3) => SiteLabel("H2", 1),
                                             SiteLabel("Ei3",  3) => SiteLabel("H2", 4),
                                             SiteLabel("Ei4",  3) => SiteLabel("H6", 1),
                                             SiteLabel("Ei5",  3) => SiteLabel("H2", 4),
                                             SiteLabel("Ei6",  3) => SiteLabel("H5", 4),
                                             SiteLabel("Ei7",  3) => SiteLabel("H4", 4),
                                             SiteLabel("Ei8",  3) => SiteLabel("H2", 4),
                                             SiteLabel("Ei9",  3) => SiteLabel("H5", 1),
                                             SiteLabel("Ei10", 3) => SiteLabel("Ei10", 3), #SiteLabel("Ei10", 3) => SiteLabel("H5", 1),

                                             SiteLabel("Ei2",  4) => SiteLabel("H2", 5),
                                             SiteLabel("Ei3",  4) => SiteLabel("H2", 3),
                                             SiteLabel("Ei4",  4) => SiteLabel("H6", 3),
                                             SiteLabel("Ei6",  4) => SiteLabel("H5", 3),
                                             SiteLabel("Ei7",  4) => SiteLabel("H4", 3),
                                             SiteLabel("Ei8",  4) => SiteLabel("H2", 4),
                                             SiteLabel("Ei9",  4) => SiteLabel("H5", 5),

                                             SiteLabel("Ei2",  5) => SiteLabel("H1", 6),
                                             SiteLabel("Ei3",  5) => SiteLabel("H1", 6),
                                             SiteLabel("Ei4",  5) => SiteLabel("H1", 6),
                                             SiteLabel("Ei6",  5) => SiteLabel("H2", 2),
                                             SiteLabel("Ei7",  5) => SiteLabel("H2", 2),
                                             SiteLabel("Ei8",  5) => SiteLabel("H2", 4),
                                             SiteLabel("Ei9",  5) => SiteLabel("H2", 2),

                                             SiteLabel("Ei2",  6) => SiteLabel("H2", 1),
                                             SiteLabel("Ei3",  6) => SiteLabel("H2", 1),
                                             SiteLabel("Ei4",  6) => SiteLabel("H2", 1),
                                             SiteLabel("Ei6",  6) => SiteLabel("H3", 1),
                                             SiteLabel("Ei7",  6) => SiteLabel("Ei7",  6),
                                             SiteLabel("Ei8",  6) => SiteLabel("Ei8",  6),
                                             SiteLabel("Ei9",  6) => SiteLabel("H4", 1))

    # Symmetric paths
    # > Note: The sectors are reflected when going back this way

    # >> Easy << >> Hard <<
    #   4  1       3  6
    # 2      5   5      2
    #   6  3       1  4

    # >> Might be worth mapping some of the E7/8/9/10 sites (on the left, which are outside the range of the sampled hard core sites) to zero.
    # >> One can argue by not doing this we are simulating how carbon moves *with* the dislocation.

    H_Ef_paths = Dict{SiteLabel,SiteLabel}( SiteLabel("H1", 1) => SiteLabel("Ef1", 4),
                                            SiteLabel("H2", 1) => SiteLabel("Ef5", 2), # Ef1/Ef2 Lower
                                            SiteLabel("H3", 1) => SiteLabel("H3", 1),
                                            SiteLabel("H4", 1) => SiteLabel("Ef7", 6),
                                            SiteLabel("H5", 1) => SiteLabel("Ef6", 6), # Ef1/Ef2 Left
                                            SiteLabel("H6", 1) => SiteLabel("Ef4", 6),
                                            SiteLabel("H7", 1) => SiteLabel("H7", 1), # Can be mapped to same position

                                            SiteLabel("H2", 2) => SiteLabel("Ef2", 5),
                                            SiteLabel("H3", 2) => SiteLabel("Ef6", 5),
                                            SiteLabel("H4", 2) => SiteLabel("Ef9", 5),
                                            SiteLabel("H5", 2) => SiteLabel("Ef9", 1),
                                            SiteLabel("H6", 2) => SiteLabel("Ef4", 1),
                                            SiteLabel("H7", 2) => SiteLabel("Ef5", 3),

                                            SiteLabel("H2", 3) => SiteLabel("Ef6", 4),
                                            SiteLabel("H3", 3) => SiteLabel("H3", 3), # Can really be mapped outside or to the same site and bulk binding
                                            SiteLabel("H4", 3) => SiteLabel("H4", 3),
                                            SiteLabel("H5", 3) => SiteLabel("H5", 3),
                                            SiteLabel("H6", 3) => SiteLabel("Ef10", 4),
                                            SiteLabel("H7", 3) => SiteLabel("H7", 3),

                                            SiteLabel("H2", 4) => SiteLabel("Ef2", 3),
                                            SiteLabel("H3", 4) => SiteLabel("Ef6", 3),
                                            SiteLabel("H4", 4) => SiteLabel("Ef9", 3),
                                            SiteLabel("H5", 4) => SiteLabel("Ef9", 6),


                                            SiteLabel("H2", 5) => SiteLabel("Ef6", 2),
                                            SiteLabel("H3", 5) => SiteLabel("H3", 5),
                                            SiteLabel("H4", 5) => SiteLabel("H4", 5),
                                            SiteLabel("H5", 5) => SiteLabel("H5", 5),


                                            SiteLabel("H2", 6) => SiteLabel("Ef5", 1),
                                            SiteLabel("H3", 6) => SiteLabel("H3", 6),
                                            SiteLabel("H4", 6) => SiteLabel("Ef7", 1),
                                            SiteLabel("H5", 6) => SiteLabel("Ef6", 1))



    H_Ei_paths = convert_site_site_dict_to_other_core(H_Ef_paths)
    Ef_H_paths = convert_site_site_dict_to_other_core(Ei_H_paths)


    return Ei_H_paths, H_Ei_paths, H_Ef_paths, Ef_H_paths
end

# ////////////////////////////////////////////////////////////////////////////////
# >>>>>>>>>>          Utility functions for manipulating sites          <<<<<<<<<<
# ////////////////////////////////////////////////////////////////////////////////

function convert_site_site_dict(paths, label_to_position)
    return hcat( [ vcat(convert(key, label_to_position), convert(value, label_to_position)) for (key,value) in paths]...)
end

function obtain_trap_mappings()
    label_to_position = create_trap_labels_to_positions_dict()
    paths = trap_site_paths()
    return map( x -> convert_site_site_dict(x, label_to_position), paths)
end


function mask_duplicate_rows(array)
    N = size(array, 1)
    mask = ones(Bool, size(array) )
    checked = []

    for i in 1:N
        pos1 = array[i,:]
        for j in 1:N
            if j != i
                if !( [j,i] in checked)
                    pos2 = array[j,:]
                    if norm( pos1 - pos2 ) < 1e-1
                        # Remove this position
                        mask[i,:] .= [ false, false ]
                    end
                    push!( checked, [i,j])
                end
            end
        end
    end
    return mask
end

function remove_duplicate_rows(array)
    mask = mask_duplicate_rows(array)
    N_mask = sum(mask[:,1])
    new_data = zeros(Float64, (Int64(N_mask), 2) )
    new_data[:,:] .= reshape( array[mask], size(new_data) )
    return new_data
end

function mask_duplicate_columns(array)
    N = size(array, 2)
    mask = ones(Bool, size(array) )
    checked = []

    for i in 1:N
        pos1 = array[:,i]
        for j in 1:N
            if j != i
                if !( [j,i] in checked)
                    pos2 = array[:,j]
                    if norm( pos1 - pos2 ) < 1e-1
                        # Remove this position
                        mask[:,i] .= [ false, false, false, false ]
                    end
                    push!( checked, [i,j])
                end
            end
        end
    end
    return mask
end

function remove_duplicate_columns(array)
    mask = mask_duplicate_columns(array)
    N_mask = sum(mask[1,:])
    new_data = zeros(Float64, (Int64(N_mask), 1) )
    new_data[:,:] .= reshape( array[mask], size(new_data) )
    return new_data
end

function convert_sitelabel_to_pos_function()
    label_to_pos = create_trap_labels_to_positions_dict()
    return x -> convert(x::SiteLabel, label_to_pos)
end

function find_self_mapped_sites(trap_site_paths)
    return filter((k,v) -> k==v, trap_site_paths)
end

function combine_mapped_sites(trap_site_paths)
    # Want the Ei state to have the self-mapped Ef states and vice versa
    # Want the H state ot have all self-mapped Ei/Ef sites
    # trap_site_paths is a tuple of the dictionaries defined above
    # > Ei_H_paths, H_Ei_paths, H_Ef_paths, Ef_H_paths

    self_mapped_Ei_H, self_mapped_H_Ei, self_mapped_H_Ef, self_mapped_Ef_H = map(find_self_mapped_sites, trap_site_paths)

    # Want the Ei state to have the self-mapped Ef states
    Ei_H_paths = merge( trap_site_paths[1], self_mapped_Ef_H )
    Ef_H_paths = merge( trap_site_paths[3], self_mapped_Ei_H )

    H_Ei_paths = merge( trap_site_paths[2], self_mapped_Ei_H, self_mapped_Ef_H )
    H_Ef_paths = merge( trap_site_paths[4], self_mapped_Ei_H, self_mapped_Ef_H )



    return Ei_H_paths, H_Ei_paths, H_Ef_paths, Ef_H_paths
end


function plot_trap_mappings()
    Ei_H_positions, H_Ei_positions, H_Ef_positions, Ef_H_positions = obtain_trap_mappings()

    label_to_position = create_trap_labels_to_positions_dict()
    trap_paths = combine_mapped_sites(trap_site_paths())

    Ei_H_positions, H_Ei_positions, H_Ef_positions, Ef_H_positions =
        map( x -> convert_site_site_dict(x, label_to_position),  trap_paths)
    #    ( convert_site_site_dict(t, label_to_position) for t in trap_paths)
    #        map( x -> convert_site_site_dict(x, label_to_position),  trap_paths)

    abcc = 2.87 # * 1.88971616463207
    alat = √2 * abcc
    lengths = [ √(3) * alat, alat]

    Ef_core_position = [ 2/6.*lengths[1], 0.]
    H_core_position = [ 1/6.*lengths[1], 1/6.*lengths[2] ]

    # Ei_positions = remove_duplicate_rows(get_all_sites_for_core("easy"))
    # Ef_positions = remove_duplicate_rows(get_all_sites_for_core("easy")  .+ Ef_core_position')
    # H_positions  = remove_duplicate_rows(get_all_sites_for_core("hard") .+ H_core_position')


    pyplot( xlims = (-4*√3 + 0.5, 4*√3 - 0.5),
            ylims = (-4*√3 + 0.2, 4*√3 - 0.8),
            size=(600,500),
            xticks=nothing,
            yticks=nothing
            )

    fnt = Plots.font( "Helvetica", 30 )
    default(titlefont = fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

    colors = palette(:tab10)
    ls = [:solid, :dash, :dot, :dashdot]

    Ei_core_position = zeros(2)
    H_core_position = [ 1/6.*lengths[1], 1/6.*lengths[2] ]
    Ef_core_position = [ 2/6.*lengths[1], 0.]


    lattice = create_lattice()

    p = scatter(lattice[:,1],lattice[:,2],
            markershape = :circle,
            markersize  = 20,
            markeralpha = 0.5,
            markercolor = :grey,
            # markerstrokewidth = 2,
            # markerstrokecolor = :black,
            markerstrokestyle = :dot,
            xlims = (-4*√3, 5*√3 ),
            ylims = ( -(4 + 1/3)*√3, (4 + 3/3.)*√3),
            ticks=nothing, axis=nothing
            )



    xi = Ei_H_positions[1,:]
    yi = Ei_H_positions[2,:]
    xf = Ei_H_positions[3,:]
    yf = Ei_H_positions[4,:]

    xi = H_Ei_positions[1,:]
    yi = H_Ei_positions[2,:]
    xf = H_Ei_positions[3,:]
    yf = H_Ei_positions[4,:]

    xi = H_Ef_positions[1,:]
    yi = H_Ef_positions[2,:]
    xf = H_Ef_positions[3,:]
    yf = H_Ef_positions[4,:]

    xi = Ef_H_positions[1,:]
    yi = Ef_H_positions[2,:]
    xf = Ef_H_positions[3,:]
    yf = Ef_H_positions[4,:]

    # scatter!(xi, yi,
    #              markershape = :circle,
    #              markersize  = 10,
    #              markeralpha = 0.5,
    #          markercolor = :blue,
    #          markerstrokewidth = 2,leg=false
    #          # markerstrokecolor = :black
    #              )

    # scatter!(xf, yf,
    #              markershape = :circle,
    #              markersize  = 10,
    #              markeralpha = 0.5,
    #          markercolor = :red,
    #          markerstrokewidth = 2
    #              )

    scatter!(Ei_positions[:,1], Ei_positions[:,2],
                 markershape = :circle,
                 markersize  = 10,
                 markeralpha = 0.5,
             markercolor = :blue,
             markerstrokewidth = 2
                 )


    scatter!(H_positions[:,1], H_positions[:,2],
                 markershape = :circle,
                 markersize  = 10,
                 markeralpha = 0.5,
             markercolor = :green,
             markerstrokewidth = 2
             )

    scatter!(Ef_positions[:,1], Ef_positions[:,2],
                 markershape = :circle,
                 markersize  = 10,
                 markeralpha = 0.5,
             markercolor = :red,
             markerstrokewidth = 2, leg=false
                 )


    dx = xf-xi
    dy = yf-yi

    arrowsize=0.7
    arrowhead=0.5
    colour = repeat(fill(1.,length(xi)), inner=4)
    quiver!(xi,yi, quiver=(dx,dy), leg=false, arrow=(arrowsize,arrowhead), color=:black)#, line_z=colour, c=:viridis)

    xi = H_Ef_positions[1,:]
    yi = H_Ef_positions[2,:]
    xf = H_Ef_positions[3,:]
    yf = H_Ef_positions[4,:]

    dx = xf-xi
    dy = yf-yi
    quiver!(xi,yi, quiver=(dx,dy), leg=false, arrow=(arrowsize,arrowhead), color=:black)#, line_z=colour, c=:viridis)


    scatter!([Ei_core_position[1], Ef_core_position[1]], [Ei_core_position[2], Ef_core_position[2]],
             markershape = :utriangle,
             markercolor = :yellow,
             markersize  = 10,
             markeralpha = 0.9,
             markerstrokewidth = 2,
             markerstrokecolor = :black)

    scatter!([H_core_position[1]], [H_core_position[2]],
             markershape = :dtriangle,
             markercolor = :yellow,
             markersize  = 10,
             markeralpha = 0.9,
             markerstrokewidth = 2,
             markerstrokecolor = :black)



    # >>>>>     H TO EF     <<<<<



    p = scatter(lattice[:,1],lattice[:,2],
            markershape = :circle,
            markersize  = 20,
            markeralpha = 0.5,
            markercolor = :grey,
            # markerstrokewidth = 2,
            # markerstrokecolor = :black,
            markerstrokestyle = :dot,
            xlims = (-4*√3, 4*√3),
            ylims = ( -(4 + 1/3)*√3, (4 + 3/3.)*√3),
            ticks=nothing, axis=nothing
            )


    xi = H_Ef_positions[1,:]
    yi = H_Ef_positions[2,:]
    xf = H_Ef_positions[3,:]
    yf = H_Ef_positions[4,:]



    # scatter!(xi, yi,
    #              markershape = :circle,
    #              markersize  = 10,
    #              markeralpha = 0.5,
    #              markercolor = :purple
    #              )

    # scatter!(xf, yf,
    #              markershape = :circle,
    #              markersize  = 10,
    #              markeralpha = 0.5,
    #              markercolor = :yellow
    #              )


    scatter!(H_positions[:,1], H_positions[:,2],
             markershape = :circle,
             markersize  = 10,
             markeralpha = 0.5,
             markercolor = :blue,
             markerstrokewidth = 2
                 )


    scatter!(Ef_positions[:,1], Ef_positions[:,2],
                 markershape = :circle,
                 markersize  = 10,
                 markeralpha = 0.5,
             markercolor = :red,
             markerstrokewidth = 2
                 )

    dx = xf-xi
    dy = yf-yi

    colour = fill(2.,length(xi))
    quiver!(xi,yi, quiver=(dx,dy), leg=false, arrow=(arrowsize,arrowhead), color=:black)#, line_z=colour, c=:viridis)

    scatter!([Ei_core_position[1], Ef_core_position[1]], [Ei_core_position[2], Ef_core_position[2]],
             markershape = :utriangle,
             markercolor = :deeppink4,
             markersize  = 10,
             markeralpha = 0.9,
             markerstrokewidth = 2,
             markerstrokecolor = :black)

    scatter!([H_core_position[1]], [H_core_position[2]],
             markershape = :dtriangle,
             markercolor = :deeppink4,
             markersize  = 10,
             markeralpha = 0.9,
             markerstrokewidth = 2,
             markerstrokecolor = :black)


    for i in 1:length(labels)
        annotate!(x[i]+0.1, y[i], Plots.text(labels[i], :black, :left, 14) )
    end


end
end
