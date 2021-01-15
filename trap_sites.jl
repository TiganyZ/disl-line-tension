
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

 """

function get_octahedral_position_dict(oct_sites, core)
    abcc = 2.87 # * 1.88971616463207
    alat = √2 * abcc

    # plat = reshape([17.32050807568877, 0.0, 0.0 ,0.0 ,16.0 ,0.0 ,0.0 ,0.0, 0.6123724356957945] * alat, (3,3) )

    # For radial configuration
    q = √(3./8.) # This should work for bcc...

    lengths = [ √(3) * alat, alat, q  * alat]

    if core == "Ei"
        core_position = zeros(3)
        labels = ["Ei1","Ei2","Ei3", "Ei4", "Ei5", "Ei6", "Ei7", "Ei8", "Ei9", "Ei10"]
    elseif core == "Ef"
        core_position = [ 2/6.*lengths[1], 0.]
        labels = ["Ef1","Ef2","Ef3", "Ef4", "Ef5", "Ef6", "Ef7", "Ef8", "Ef9", "Ef10"]
    elseif core == "H"
        core_position = [ 1/6.*lengths[1], 1/6.*lengths[2] ]
        labels = ["H1", "H2",     "H3", "H4", "H5", "H6", "H7"]
    end

    labels = [ (l,i) for i in 1:6 for l in labels]

    site_pos_dict = Dict{String, Array{Float64,1}}( labels[i] => oct_sites[i,:] + core_position for i in 1:length(labels))
    return site_pos_dict
end

function create_trap_labels_to_positions_dict()
    Ei_dict = get_octahedral_position_dict(get_octahedral_positions("easy"), "Ei")
    Ef_dict = get_octahedral_position_dict(get_octahedral_positions("easy"), "Ef")
    H_dict  = get_octahedral_position_dict(get_octahedral_positions("hard"), "H")

    return merge(Ei_dict, Ef_dict, H_dict)
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

function get_all_sites_for_core(core_type)
    oct_sites_initial, core_position, relative_core_position = get_octahedral_positions(core_type)

    ###---   Add also a REFLECTION as well as a rotation
    oct_sites  = get_all_oct_sites(oct_sites_initial, [0,0.])

    # reflect this oct,  -0.5*[√(2/3)*2.87 0.0]
    reflected_oct = zeros(size(oct_sites))
    reflected_oct = reflect_data_about_line!(1, oct_sites, reflected_oct, [0., 1.] )

    oct_sites  = vcat( oct_sites, get_all_oct_sites(reflected_oct, [0,0.]) )
    return oct_sites
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

    Quadrants are the number next to the
"""
    # >> Easy << >> Hard <<
    #   4  1       3  6
    # 2      5   5      2
    #   6  3       1  4

    Ei-H_paths = Dict( ("Ei1",  1) => (  "H1", 6), # Ei -> H sector 1
                      ("Ei2",  1) => (  "H1", 6),
                      ("Ei3",  1) => (  "H1", 6),
                      ("Ei4",  1) => (  "H1", 6),
                      ("Ei5",  1) => (  "H2", 3),
                      ("Ei6",  1) => (  "H2", 6),
                      ("Ei7",  1) => (  "H2", 6),
                      ("Ei8",  1) => (  "H2", 3),
                      ("Ei9",  1) => (  "H2", 6),
                      ("Ei10", 1) => (  "H6", 6),

                      ("Ei1",  2) => (  "H2", 5), # Ei -> H sector 2
                      ("Ei2",  2) => (  "H2", 5),
                      ("Ei3",  2) => (  "H2", 5),
                      ("Ei4",  2) => (  "H2", 5),
                      ("Ei5",  2) => (  "H7", 1),
                      ("Ei6",  2) => (  "H3", 5),
                      ("Ei7",  2) => (  "Ei7", 2),
                      ("Ei8",  2) => (  "Ei8", 2),
                      ("Ei9",  2) => (  "H4", 5),
                      ("Ei10", 2) => (  "H5", 5),

                      ("Ei1",  3) => (  "H2", 1), # Ei -> H sector 3
                      ("Ei2",  3) => (  "H2", 1),
                      ("Ei3",  3) => (  "H2", 4),
                      ("Ei4",  3) => (  "H6", 1),
                      ("Ei5",  3) => (  "H2", 4),
                      ("Ei6",  3) => (  "H5", 4),
                      ("Ei7",  3) => (  "H4", 4),
                      ("Ei8",  3) => (  "H2", 4),
                      ("Ei9",  3) => (  "H5", 1),
                      ("Ei10", 3) => (  "H5", 1),

                      ("Ei1",  4) => (  "H2", 5), # Ei -> H sector 4
                      ("Ei2",  4) => (  "H2", 5),
                      ("Ei3",  4) => (  "H2", 3),
                      ("Ei4",  4) => (  "H6", 3),
                      ("Ei5",  4) => (  "H2", 3),
                      ("Ei6",  4) => (  "H5", 3),
                      ("Ei7",  4) => (  "H4", 3),
                      ("Ei8",  4) => (  "H2", 4),
                      ("Ei9",  4) => (  "H5", 5),
                      ("Ei10", 4) => (  "Ei10", 4),

                      ("Ei1",  5) => (  "H1", 6), # Ei -> H sector 5
                      ("Ei2",  5) => (  "H1", 6),
                      ("Ei3",  5) => (  "H1", 6),
                      ("Ei4",  5) => (  "H1", 6),
                      ("Ei5",  5) => (  "H2", 4),
                      ("Ei6",  5) => (  "H2", 2),
                      ("Ei7",  5) => (  "H2", 2),
                      ("Ei8",  5) => (  "H2", 4),
                      ("Ei9",  5) => (  "H2", 2),
                      ("Ei10", 5) => (  "H6", 2),

                      ("Ei1",  6) => (  "H2", 1), # Ei -> H sector 5
                      ("Ei2",  6) => (  "H2", 1),
                      ("Ei3",  6) => (  "H2", 1),
                      ("Ei4",  6) => (  "H2", 1),
                      ("Ei5",  6) => (  "H7", 1),
                      ("Ei6",  6) => (  "H3", 1),
                      ("Ei7",  6) => (  "Ei7",  6),
                      ("Ei8",  6) => (  "Ei8",  6),
                      ("Ei9",  6) => (  "H4", 1),
                      ("Ei10", 6) => (  "H4", 2) )

    # >> Might be worth mapping some of the E7/8/9/10 sites (on the left, which are outside the range of the sampled hard core sites) to zero.
    # >> One can argue by not doing this we are simulating how carbon moves *with* the dislocation.


    # >> Easy << >> Hard <<
    #   4  1       3  6
    # 2      5   5      2
    #   6  3       1  4

    #  REMEMBER, these are translated sites.

    # If a site is mappee to the same position we can set the criteria

    H-Ef_paths = Dict( ("H1", 1) => ("Ef1", 1),
                       ("H2", 1) => ("Ef5", 2), # Ef1/Ef2 Lower
                       ("H3", 1) => ("H3", 1),
                       ("H4", 1) => ("Ef7", 6),
                       ("H5", 1) => ("Ef6", 6), # Ef1/Ef2 Left
                       ("H6", 1) => ("Ef4", 6),
                       ("H7", 1) => ("Ef8", 2), # Can be mapped to same position

                       ("H1", 2) => ("Ef1", 1),
                       ("H2", 2) => ("Ef2", 5),
                       ("H3", 2) => ("Ef6", 5),
                       ("H4", 2) => ("Ef9", 5),
                       ("H5", 2) => ("Ef9", 1),
                       ("H6", 2) => ("Ef4", 1),
                       ("H7", 2) => ("Ef6", 3),

                       ("H1", 3) => ("Ef1", 1),
                       ("H2", 3) => ("Ef6", 4),
                       ("H3", 3) => ("H3", 3), # Can really be mapped outside or to the same site and bulk binding
                       ("H4", 3) => ("H4", 3),
                       ("H5", 3) => ("H5", 3),
                       ("H6", 3) => ("Ef10", 4),
                       ("H7", 3) => ("H7", 3),

                       ("H1", 4) => ("Ef1", 1),
                       ("H2", 4) => ("Ef2", 3),
                       ("H3", 4) => ("Ef6", 3),
                       ("H4", 4) => ("Ef9", 3),
                       ("H5", 4) => ("Ef9", 6),
                       ("H6", 4) => ("Ef4", 6),
                       ("H7", 4) => ("Ef5", 3),

                       ("H1", 5) => ("Ef1", 1),
                       ("H2", 5) => ("Ef6", 2),
                       ("H3", 5) => ("H3", 5),
                       ("H4", 5) => ("H4", 5),
                       ("H5", 5) => ("H5", 5),
                       ("H6", 5) => ("Ef10", 2),
                       ("H7", 5) => ("H7", 5),

                       ("H1", 6) => ("Ef1", 1), # Sector 6
                       ("H2", 6) => ("Ef5", 1),
                       ("H3", 6) => ("H3", 6),
                       ("H4", 6) => ("Ef7", 1),
                       ("H5", 6) => ("Ef6", 1),
                       ("H6", 6) => ("Ef4", 1),
                       ("H7", 6) => ("H7", 6) )


    return merge(Ei-H_paths, H-Ef_paths)
end


function obtain_trap_mappings()
    trap_label_to_position = create_trap_labels_to_positions_dict()
    trap_site_paths = trap_site_paths()


end
