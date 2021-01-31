using Plots

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


# First read the solute positions and the occupancies from the files
function read_image_solutes_occupancies(file, n, n_points)
    open(file) do f
        for i in eachline(f)
            if ismatch("Image $n", i)
                # Read the solute positons where each point is demarcated by space
                readline()
                # Now read dlm
                pos_occ = []
                for j in 1:n_points
                    push!(pos_occ, read_solute_positions_for_point(f))
                    if j == 1
                        n_solutes = ceil(Int64,length(pos_occ[1])/3)
                    end
                end
            end
        end
    end
    # Array is 3xN, where each row corresponds to x, y, and occ.
    return reshape(hcat(pos_occ...), reshape(3,n_solutes,n_points))
end

function read_solute_positions_for_point(f)
    lines = []
    line = readline(f)
    while line != ""
        push!(lines, parse.(Float64, split(line) ))
        line = readline(f)
    end
    return hcat(lines)
end

function plot_solute_positions(lattice, x, y, occ)
    abcc = 2.87 # * 1.88971616463207
    alat = √2 * abcc
    lengths = [ √(3) * alat, alat]

    Ei_core_position = zeros(2)
    H_core_position = [ 1/6.*lengths[1], 1/6.*lengths[2] ]
    Ef_core_position = [ 2/6.*lengths[1], 0.]

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

    scatter!(x, y,
             markershape = :circle,
             markersize  = 10,
             markeralpha = 0.5,
             # markerstrokewidth = 2,leg=false,
             m = ColorGradient(:grayC),
             zcolor=occ
             # markerstrokecolor = :black
             )

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
    return p
end



file     = ARGS[1]
n        = 15 # ARGS[2]
n_points = 91 # ARGS[3]

pos_occ = read_image_solutes_occupancies(file, n, n_points)

pyplot( xlims = (-4*√3 + 0.5, 4*√3 - 0.5),
        ylims = (-4*√3 + 0.2, 4*√3 - 0.8),
        size=(600,500),
        xticks=nothing,
        yticks=nothing,
        clims =( 0, 1. )
        )

fnt = Plots.font( "Helvetica", 30 )
default(titlefont = fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

colors = palette(:tab10)
ls = [:solid, :dash, :dot, :dashdot]



lattice = create_lattice()

# Filter the points that I actually want to plot
points = [1, 46, 91]

pl = []
l = @layout [ a b c ]

for point in point
    x   = pos_occ[point,:,1]
    y   = pos_occ[point,:,2]
    occ = pos_occ[point,:,3]

    push!(pl,plot_solute_positions(lattice, x, y, occ))
end

plot(pl..., layout=l)
