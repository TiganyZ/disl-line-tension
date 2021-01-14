
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

 """

function get_octahedral_position_dict(core_type)
    abcc = 2.87 # * 1.88971616463207
    alat = √2 * abcc

    # plat = reshape([17.32050807568877, 0.0, 0.0 ,0.0 ,16.0 ,0.0 ,0.0 ,0.0, 0.6123724356957945] * alat, (3,3) )

    # For radial configuration
    q = √(3./8.) # This should work for bcc...

    lengths = [ √(3) * alat, alat, q  * alat]


    easy_core_position = [36.3218504,  33.14680888]
    easy_relative_core_position = [ 1/6.*lengths[1], 1/6.*lengths[2] ]

    hard_core_position = [37.493523, 33.82327437]
    hard_relative_core_position = [ 2/6.*lengths[1], 2/6.*lengths[2] ]


    easy_octahedral_positions =[[37.54770758  33.85046096];
                                [37.50566447  34.41147207];
                                [36.92444498  35.52976750];
                                [38.12738250  35.54441951];
                                [36.32784633  36.51662541];
                                [38.68239733  36.53429325];
                                [38.09544321  37.50275944];
                                [36.92344511  37.51493620];
                                [39.28127528  35.50886674];
                                [40.51130213  35.56288039]]   .- easy_core_position'


    # Assuming here that the hard core is offset from the easy core.
    hard_octahedral_positions = [[  37.49967197  33.82188364];
                                 [  35.73715334  31.44012169];
                                 [  34.57273597  29.41591641];
                                 [  35.75702326  29.40559471];
                                 [  36.89155022  29.38012684];
                                 [  37.49563378  30.34278051];
                                 [  33.40750846  31.46033380]] .- easy_core_position'

    hard_labels = ["H1", "H2",     "H3", "H4", "H5", "H6", "H7"]
    easy_labels = ["E1","E2","E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10"]
    labels = vcat( easy_labels, hard_labels )

    positions = vcat(easy_octahedral_positions, hard_octahedral_positions)
    site_pos_dict = Dict{String, Array{Float64,1}}( labels[i] => positions[i,:] for i in 1:length(labels))
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
    elseif core_type == "hard"
        core_position = [37.493523, 33.82327437]
        relative_core_position = [ 2/6.*lengths[1], 2/6.*lengths[2] ]
    end

    easy_octahedral_positions =[[37.54770758  33.85046096];
                                [37.50566447  34.41147207];
                                [36.92444498  35.52976750];
                                [38.12738250  35.54441951];
                                [36.32784633  36.51662541];
                                [38.68239733  36.53429325];
                                [38.09544321  37.50275944];
                                [36.92344511  37.51493620];
                                [39.28127528  35.50886674];
                                [40.51130213  35.56288039]]   .- core_position'


    hard_octahedral_positions = [[  37.49967197  33.82188364];
                                 [  35.73715334  31.44012169];
                                 [  34.57273597  29.41591641];
                                 [  35.75702326  29.40559471];
                                 [  36.89155022  29.38012684];
                                 [  37.49563378  30.34278051];
                                 [  33.40750846  31.46033380]] .- core_position'



    return easy_octahedral_positions, hard_octahedral_positions, core_position, relative_core_position
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
    oct_sites_initial, hard_oct, easy_core_position, relative_easy_core_position = get_octahedral_positions(core_type)

    ###---   Add also a REFLECTION as well as a rotation
    oct_sites  = get_all_oct_sites(oct_sites_initial, [0,0.])

    # reflect this oct,  -0.5*[√(2/3)*2.87 0.0]
    reflected_oct = zeros(size(oct_sites))
    reflected_oct = reflect_data_about_line!(1, oct_sites, reflected_oct, [0., 1.] )

    oct_sites  = vcat( oct_sites, get_all_oct_sites(reflected_oct, [0,0.]) )
end

function get_labels_for_core(core_type)
    if core_type == "easy"
        return ["E1","E2","E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10"]
    elseif core_type == "hard"
        return ["H1", "H2",     "H3", "H4", "H5", "H6", "H7"]
    end
end

function get_all_labels_for_core(core_type)
    # Action of getting sites are rotating by 0, 2π/3, and 4π/3 then a reflection of all of them
    # Sites are partitioned into 6 sectors, which are labelled from 1-6 going anticlockwise
    # >> Counter <<
    #   2  1
    # 3      6
    #   4  5

    # >> Easy <<
    #   4  1
    # 2      5
    #   6  3

    # >> Hard <<
    #   3  6
    # 5      2
    #   1  4

    # So order of labels is 1, 3, 5, 2, 6, 4
    easy_sector_order = [1, 3, 5, 2, 6, 4]
    hard_sector_order = [4, 6, 2, 5, 3, 1]




end

function trap_site_paths()
    """ These are the paths that Ivo has
    chosen. I can recreate his results or do something which is a bit
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

    paths = [["E1",  1,  "H1", 6], # E -> H sector 1
             ["E2",  1,  "H1", 6],
             ["E3",  1,  "H1", 6],
             ["E4",  1,  "H1", 6],
             ["E5",  1,  "H2", 3],
             ["E6",  1,  "H2", 6],
             ["E7",  1,  "H2", 6],
             ["E8",  1,  "H2", 3],
             ["E9",  1,  "H2", 6],
             ["E10", 1,  "H6", 6],


             ["E1",  2,  "H2", 5], # E -> H sector 2
             ["E2",  2,  "H2", 5],
             ["E3",  2,  "H2", 5],
             ["E4",  2,  "H2", 5],
             ["E5",  2,  "H7", 1],
             ["E6",  2,  "H3", 5],
             ["E7",  2,  "H3", 5],
             ["E8",  2,  "H7", 1],
             ["E9",  2,  "H4", 5],
             ["E10", 2,  "H5", 5],


             ["E1",  3,  "H2", 1], # E -> H sector 3
             ["E2",  3,  "H2", 1],
             ["E3",  3,  "H2", 4],
             ["E4",  3,  "H6", 1],
             ["E5",  3,  "H2", 4],
             ["E6",  3,  "H5", 4],
             ["E7",  3,  "H4", 4],
             ["E8",  3,  "H2", 4],
             ["E9",  3,  "H5", 1],
             ["E10", 3,  "H5", 1],



             ["E1",  4,  "H2", 5], # E -> H sector 4
             ["E2",  4,  "H2", 5],
             ["E3",  4,  "H2", 3],
             ["E4",  4,  "H6", 3],
             ["E5",  4,  "H2", 3],
             ["E6",  4,  "H5", 3],
             ["E7",  4,  "H4", 3],
             ["E8",  4,  "H2", 4],
             ["E9",  4,  "H5", 5],
             ["E10", 4,  "H5", 5],

    # >> Easy << >> Hard <<
    #   4  1       3  6
    # 2      5   5      2
    #   6  3       1  4

             
             ["E1",  5,  "H1", 6], # E -> H sector 5
             ["E2",  5,  "H1", 6],
             ["E3",  5,  "H1", 6],
             ["E4",  5,  "H1", 6],
             ["E5",  5,  "H2", 4],
             ["E6",  5,  "H2", 2],
             ["E7",  5,  "H2", 2],
             ["E8",  5,  "H2", 4],
             ["E9",  5,  "H2", 2],
             ["E10", 5,  "H6", 2],



             ["E1",  6,  "H2", 1], # E -> H sector 5
             ["E2",  6,  "H2", 1],
             ["E3",  6,  "H2", 1],
             ["E4",  6,  "H2", 1],
             ["E5",  6,  "H7", 1],
             ["E6",  6,  "H3", 1],
             ["E7",  6,  "H2", 1],
             ["E8",  6,  "H7", 1],
             ["E9",  6,  "H4", 1],
             ["E10", 6,  "H4", 2],

             # >> Might be worth mapping some of the E7/8/9/10 sites to zero which are outside the range of the sampled hard core sites.
             # >> One can argue by not doing this we are simulating how carbon moves _with_ the dislocation.

             ["E2", 5,  "H1", 6],
             ["E1", 1,  "H2", 1], # E1/E2 Lower
             ["E2", 1,  "H2", 1],
             ["E2", 1,  "H2", 1],
             ["E1", 1,  "H2", 1], # E1/E2 Left
             ["E2", 1,  "H2", 1],
             ["E2", 1,  "H2", 1],
             ["E3", 1,  "H2", 1], # E3/E4 Right
             ["E4", 1,  "H1", 1],



             ["E1", "H1", "E1"], # >> Upper right quadrant with core going to right
             ["E2", "H1", "E2"],
             ["E3", "H1", "E1"], # Because of the reconstruction of the hard core
             ["E4", "H1", "E1"], # > Likewise
             ["E5", "H2", "E6"],
             ["E6", "H2", "E5"],
             ["E7", "H2", "E6"],
             ["E8", "H2", "E6"],
             ["E9", "H2", "E4"],

             ["E1", "H1", "E1"],
             ["E1", "H1", "E1"],
             ["E1", "H1", "E1"],
             ["E1", "H1", "E1"],
             ["E1", "H1", "E1"],
             ["E1", "H1", "E1"],
             ["E1", "H1", "E1"],
             ["E1", "H1", "E1"],
             ]

end
