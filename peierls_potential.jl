module PeierlsPotential

using PyPlot
using DSP
using PyCall
@pyimport scipy.interpolate as si


export get_2D_peierls_profile, Potential

struct Potential{Function} 
    ΔEₚ::Function
    ∂ΔEₚ::Array{Function}
    hess::Array{Function,2}
end


function get_2D_peierls_profile(;pp_type="sdTB")
    if pp_type == "itakura"
        f, dx, dd = itakura_data()
    else
        f, dx, dd = pp_data(pp_type)
    end
    
    return f, dx, dd
end



function get_xy_arrays!(x, y, idx, v_range)
    for i in 1:length(idx)
        x[idx[i]] = v_range[i][1]
        y[idx[i]] = v_range[i][2]    
    end
    return x, y
end


function reflect_data_about_line!(index::Int64, data::Array{Float64,2}, all_data::Array{Float64,2}, line::Array{Float64,1}, translation = [0.,0.]; no_idx = [])
    trans = zeros(3)
    trans[1:2] = translation
    for i in range( 1, size(data)[1] )
        if i in no_idx
            continue
        else
            V = data[i,1:end] - trans
            rv = reflect_2D( V[1:end-1], line )

            all_data[index,1:end-1] = rv + translation
            all_data[index,3] = V[3]
            index += 1
        end
    end
    return all_data
end



function get_points_in_range(data, xlim, ylim)

    xl, xu = xlim
    yl, yu = ylim
    
    xmask = (data[:,1] .>= xl) .& (data[:,1] .<= xu)
    ymask = (data[:,2] .>= yl) .& (data[:,2] .<= yu)
    rows = xmask .& ymask
    
    return reshape( data[ repeat( rows, outer=(1,3) )  ], (sum(rows),3) )

end

function findmat(f, A::AbstractMatrix)
    m,n = size(A)
    out = []
    for i in 1:m, j in 1:n
        f( A[i,j] ) && push!( out, (i,j) )
    end
    out
end



function fft_expansion(ftd, xr, yr, nx, ny)
    # ftd is the data to fourier transform
    # xr and yr are the data ranges

    xl, xu = xr
    yl, yu = yr

    T_x = xu - xl
    T_y = yu - yl

    # ADD shift to x and y when using ft function
    shift =  [ -xl , -yl ]
    
    nfreqs = 3
    
    F = fft( ftd ) # |> fftshift

    # Low pass filter
    sz = size( F )
    center = convert( Int64,  floor( nx / 2 ) + 1 )
    half_width = convert( Int64, nfreqs )

    filter =  zeros( sz )
    filter[ center - half_width : center + half_width, center - half_width : center + half_width  ] .= 1.0
    filter = ifftshift( filter )

    Fn = F .* filter 

    # > Can plot the new potential
    # heatmap( real( new_potential )  )
    # show()

    # Non-zero indices are

    indices = findmat( x -> x==1, filter )

    ffreq_x =  fftfreq(nx) * nx / T_x 
    ffreq_y =  fftfreq(ny) * ny / T_y 

    #// Get the complex exponentials for the FT
    exponentials    = [ (x,y) ->                           exp(im * 2π * ( (x-xl) * ffreq_x[u[1]] + (y-yl) * ffreq_y[u[2]] ) ) for u in indices  ]
    dx_exponentials = [ (x,y) -> im * 2π * ffreq_x[u[1]] * exp(im * 2π * ( (x-xl) * ffreq_x[u[1]] + (y-yl) * ffreq_y[u[2]] ) ) for u in indices  ]
    dy_exponentials = [ (x,y) -> im * 2π * ffreq_y[u[2]] * exp(im * 2π * ( (x-xl) * ffreq_x[u[1]] + (y-yl) * ffreq_y[u[2]] ) ) for u in indices  ]    

    dx_dx_exponentials = [ (x,y) -> ( -4π^2 * ffreq_x[u[1]] * ffreq_x[u[1]] * exp(im * 2π * ( (x-xl) * ffreq_x[u[1]] + (y-yl) * ffreq_y[u[2]] )   )) for u in indices  ]
    dy_dx_exponentials = [ (x,y) -> ( -4π^2 * ffreq_y[u[2]] * ffreq_x[u[1]] * exp(im * 2π * ( (x-xl) * ffreq_x[u[1]] + (y-yl) * ffreq_y[u[2]] )   )) for u in indices  ]    
    dx_dy_exponentials = [ (x,y) -> ( -4π^2 * ffreq_x[u[1]] * ffreq_y[u[2]] * exp(im * 2π * ( (x-xl) * ffreq_x[u[1]] + (y-yl) * ffreq_y[u[2]] )   )) for u in indices  ]
    dy_dy_exponentials = [ (x,y) -> ( -4π^2 * ffreq_y[u[2]] * ffreq_y[u[2]] * exp(im * 2π * ( (x-xl) * ffreq_x[u[1]] + (y-yl) * ffreq_y[u[2]] )   )) for u in indices  ]    

    
    coefficients = [ Fn[ u[1], u[2] ] / (nx * ny) for u in indices ]


    fexp(x,y)    = sum(  real( coefficients[i] *    exponentials[i]( x,y ) ) for i in 1:length(coefficients)  )
    dx_fexp(x,y) = sum(  real( coefficients[i] * dx_exponentials[i]( x,y ) ) for i in 1:length(coefficients)  )
    dy_fexp(x,y) = sum(  real( coefficients[i] * dy_exponentials[i]( x,y ) ) for i in 1:length(coefficients)  )


    dx_dx_fexp(x,y) = sum(  real( coefficients[i] * dx_dx_exponentials[i]( x,y ) ) for i in 1:length(coefficients)  )
    dy_dx_fexp(x,y) = sum(  real( coefficients[i] * dy_dx_exponentials[i]( x,y ) ) for i in 1:length(coefficients)  )
    dx_dy_fexp(x,y) = sum(  real( coefficients[i] * dx_dy_exponentials[i]( x,y ) ) for i in 1:length(coefficients)  )
    dy_dy_fexp(x,y) = sum(  real( coefficients[i] * dy_dy_exponentials[i]( x,y ) ) for i in 1:length(coefficients)  )
    
    plotfft=false
    if plotfft
        it = 0
        for xi in LinSpace( xl, xu , 30)
            for yi in LinSpace( yl, yu, 30 )
                it +=  1
                x[it] = xi 
                y[it] = yi 
                f[it] = fexp( xi , yi  ) 
                # println("it = $(it): x = $(x[it]), y = $(y[it]), ft = $(f[it])" )
            end
        end

        # scatter3D( new_data[:,1], new_data[:,2], new_data[:,3])
        surf( x, y, f, alpha=0.9, cmap = ColorMap("magma") )
        PyPlot.show()
    end

    grad = [ dx_fexp, dy_fexp ]
    
    hessian = [dx_dx_fexp dx_dy_fexp;
               dy_dx_fexp dy_dy_fexp ]
    
    return fexp, grad, hessian
end


function append_if_close!(X,d,i; tol=1e-4)
    if length(X) !== 0 
        if !any( X .- d[i] .< tol)
            append!(X, d[i])
        end
    else
        append!(X, d[i])
    end
    return X  
end

function mask_duplicates(all_interp)
    
    N = size(all_interp, 2)

    mask = ones(Bool, size(all_interp) )
    checked = []
    for i in 1:N
        idx1 = i
        pos1 = all_interp[ idx1, 1:2]
        sorted_idx = nearest_neighbours_to_pos( pos1', all_interp[:,1:2] )

        for j in 2:10
            idx2 = Int64(sorted_idx[j])

            # if idx2 < idx1
            #     continue
            # end

            if !( [idx2, idx1] in checked)
                pos2 = all_interp[ idx2, 1:2 ]
                if norm( pos1 - pos2 ) < 1e-3
                    # Remove this position
                    #                    println("Removing $i, $j")
                    mask[ idx2, :] .= [ false, false, false ]
                    push!( checked, [idx1, idx2])
                end
                
            end
        end
    end

    #    println("\n  Removed $(length(checked)) indices")
    return mask
end



function nearest_neighbours_to_pos( pos, lattice )
    # Broadcast to get subtraction of all points from all others
    N, M = size(lattice)
    return sortperm( reshape( sum( ( lattice .- pos ).^2 , 2 ), (N,) ) )
end



function nearest_neighbours( lattice )
    # Broadcast to get subtraction of all points from all others
    N, M = size(lattice)
    dist = zeros(N,N)
    dist .= reshape( sum( (      reshape(lattice, ( N, 1, M ) )
                              .- reshape(lattice, ( 1, N, M ) ) ).^2, 3 ) , (N,N) )
    # Make self neighbours a large distance away
    dist[ [ CartesianIndex(i,i) for i in 1:N ] ] .= 9e10

    # Sort perm each slice to get all the sorted neighbours
    sorted_idx = zeros(N,N)

    for i in 1:N
        sorted_idx[i,:] .= sortperm( dist[i,:] )
    end
    return sorted_idx
end



function reflect_2D( v::Array{Float64, 1}, l::Array{Float64, 1} )
    # v is vector being reflected
    # l is line in reflection plane
    ref_v = 2 * ( v' * l / (l' * l) )*l - v
    return ref_v
    
end


function pp_data(pp_type)
    # These are uncorrected energies. 
    
     # ΔEₚ = [0                       ,
     #        3.65082532905568403354  ,
     #        9.29815470646174823531  ,
     #        12.92467059693784315456 ,
     #        13.45124433806511119728 ,
     #        25.06252842878757068658 ,
     #        50.68441172749887316207 ,
     #        77.53293787911755489605 ,
     #        90.40503002259714072790 ,
     #        34.89064341817444790128 ,
     #        55.09639757819516543199 ,
     #        75.00662762782350888909 ,
     #        29.91981038926772858294 ,
     #        55.05561588962403666901 ,
     #        33.89145161059813450940  ] 


    #  Pure, unadulterated tbe energies
    ΔEₚ =[0,
          7.28292608920528150000,
          16.0263558372945985000,
          22.2295337279900560000,
          24.8292424859318985000,
          3.01461452046146300000,
          7.13618863556455000000,
          13.0119454023041240000,
          5.41738079889940300000,
          22.1115042972683735000,
          18.1917707129729130000,
          14.0259100506653215000,
          11.4987876919050555000,
          15.1033452805039425000,
          18.5540224239775005000 ]


    if pp_type == "sdTB"
        ΔEₚ =[0,
              7.98660445472177013708  ,
              17.38583648684110470790,
              24.19724779326713327748,
              27.35839739416652913876,
              6.31932233790614794103  ,
              13.62982749656252740248 ,
              22.57877423329679557779,
              17.94145187561050074170 ,
              26.88027044255292307409,
              25.37217688298149369611,
              23.79111856299613622333,
              15.34109098202042633111,
              21.96630168047209769439,
              22.88416331017000554862 ]
    elseif pp_type == "dTB"
        ΔEₚ =[0,
              6.3,
              15.1,
              20.4,
              22.6,
              4.6,
              12.7,
              22.7,
              26.8,
              23.0,
              23.5,
              24.4,
              13.2,
              20.3,
              20.0 ]
    end

    # ΔEₚ = [0    ,
    #         36.65694624150005743962  ,
    #         74.78764899514385721565  ,
    #         110.39417477085900585558 ,
    #         142.41613263213650907903 ,
    #         64.68508378701671282553  ,
    #         130.63061291074666481024 ,
    #         198.48105894590046131036 ,
    #          253.00923890856024992561,
    #          172.42798345683760898504,
    #          201.08376470679133173949,
    #          229.34212824224212898088,
    #          102.70622250950392397572,
    #          168.29243304089817397884,
    #          139.30920069214012218680 ] 

    n_points = length( ΔEₚ )

    x = zeros(n_points)
    y = zeros(n_points)

    #  Easy core is at       0.0,  0.0
    #  Hard core is at sqrt(3)/2,  0.5
    # Split core is at sqrt(3)/2, -0.5

    easy  = [0.0, 0.0]
    hard  = [√(3)/2,   0.5]
    split = [√(3)/2,  -0.5]

    # E -> H: 
    idx_eh =  collect(1:5)
    eh = LinSpace(easy, hard, 5 )
    x, y = get_xy_arrays!(x, y, idx_eh, eh)
    
    # E -> S:  0, 2, 5,  9, 18
    idx_es =  [1, 6, 7, 8, 9]
    es = LinSpace(easy, split, 5 )
    x, y = get_xy_arrays!(x, y, idx_es, es)

    # H -> S:  11--17
    idx_hs =  [5, 10, 11, 12, 9 ]
    hs = LinSpace(hard, split, 5 )
    x, y = get_xy_arrays!(x, y, idx_hs, hs)    

    # Rest
    # P4 == 0.5( P3 + P5 )
    x[13] = 0.5 * ( x[3] + x[7])
    y[13] = 0.5 * ( y[3] + y[7])     

    x69 = LinSpace(x[4], x[8], 4)
    y69 = LinSpace(y[4], y[8], 4)
    
    x[14], y[14]  = x69[2], y69[2]
    x[15], y[15]  = x69[3], y69[3]        
    

    z = ΔEₚ
    
    original_data = zeros(n_points,3)

    original_data[:,1] .= x 
    original_data[:,2] .= y
    original_data[:,3] .= z
    

    all_data= zeros( n_points * 6 * 4 * 2, 3 )
    irr_data = zeros( 110, 3 )
    
    all_data[1:n_points, 1:3] =  original_data
    irr_data[1:n_points, 1:3] =  original_data    
    
    triangle1 = zeros(n_points,3)
    triangle2 = zeros(n_points,3)
    triangle3 = zeros(n_points,3)
    triangle4 = zeros(n_points,3)
    triangle5 = zeros(n_points,3)
    triangle6 = zeros(n_points,3)

    
    triangle1  = original_data
    index = n_points + 1
    
    all_data     = reflect_data_about_line!(index, original_data, all_data, hard)
    triangle2    = all_data[ index: index + n_points - 1, : ]
    index += n_points

    all_data     = reflect_data_about_line!(index, triangle1, all_data, split)
    triangle4    = all_data[index : index + n_points - 1, : ]
    index += n_points

    
    # all_data     = reflect_data_about_line!(index, triangle2, all_data, split, [0.,1])    
    # triangle3    = all_data[index : index + n_points - 1, : ]
    # index += n_points


    all_data = reflect_data_about_line!(index, all_data[1:index-1,1:end], all_data, [0., 1.] )
    index += 3 * n_points

    all_data = reflect_data_about_line!(index, all_data[1:index-1,1:end], all_data, [0., 1.], [ √3 / 2, 0. ]  )
    all_data[index + 6*n_points   : index + 6*n_points*2 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [  √3 / 2,  1.5, 0.0], (1,3) )
    all_data[index + 6*n_points*2 : index + 6*n_points*3 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [  √3 / 2, -1.5, 0.0], (1,3) )

    all_data[index + 6*n_points*3 : index + 6*n_points*4 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [ 3*√3 / 2,  1.5, 0.0], (1,3) )
    all_data[index + 6*n_points*4 : index + 6*n_points*5 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [ 3*√3 / 2, -1.5, 0.0], (1,3) )

    all_data[index + 6*n_points*5 : index + 6*n_points*6 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [ -√3 / 2,  1.5, 0.0], (1,3) )
    all_data[index + 6*n_points*6 : index + 6*n_points*7 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [ -√3 / 2, -1.5, 0.0], (1,3) )
    
    # index += 6 * n_points
    # all_data = reflect_data_about_line!(index, all_data[1:index-1,1:end], all_data, [1., 0.], [ 0, 0.5 ]  )

    # Now mask to remove duplicates    
    mask = mask_duplicates(all_data)
    N_mask = sum(mask[:,1])
    new_data = zeros(Float64, (Int64(N_mask), 3) )
    new_data[:,:] .= reshape( all_data[mask], size(new_data) )


    # new_data[:,2] .+= 2
    # new_data[:,2] .= mod.(new_data[:,2], 3.0)

    
    # >>> Interpolating data by grid interpolation using Scipy Interpolate and griddata

    nx=150
    ny=150
    all_interp =  zeros( nx*ny, 3 )
    ndat = get_points_in_range(new_data, (0, √3), (-2, 1))
    x_grid = LinSpace(minimum(ndat[:,1]), maximum(ndat[:,1]),  nx)
    y_grid = LinSpace(minimum(ndat[:,2]), maximum(ndat[:,2]),  ny)
    X = kron(ones(nx), x_grid)
    Y = kron(y_grid, ones(ny))
    all_interp[:,1] .= X
    all_interp[:,2] .= Y    
    all_interp[:,3] .= si.griddata((ndat[:,1], ndat[:,2]), ndat[:,3], (X, Y) , method="cubic")

    
    # Using the interpolated function as input for Plane wave expansion
    # >> Plane wave expansion is not a good representation of the energy surface, therefore can just use the griddata
    # >> Get a function for the spline interpolation and then pass into the NEB / String Method

    # Create periodic area within some bounds
    # > Array with all the original points in is new_data

    ftd = reshape( all_interp[:,3], (nx,ny) )

    xr = (0, √3)
    yr = (-2, 1)    

    fexp, grad, hessian  = fft_expansion(ftd, xr, yr, nx, ny)
        
    return fexp, grad, hessian
end


function itakura_data()
    points = 0:1:18

    #  Easy core is at       0.0,  0.0
    #  Hard core is at sqrt(3)/2,  0.5
    # Split core is at sqrt(3)/2, -0.5

    easy  = [0.0, 0.0]
    hard  = [sqrt(3)/2,   0.5]
    split = [sqrt(3)/2,  -0.5]

    n_points = 19
    x = zeros(19)
    y = zeros(19)
    
    # E -> H:  0, 1, 3, 6, 10
    idx_eh =  [1, 2, 4, 7, 11]
    eh = LinSpace(easy, hard, 5 )
    x, y = get_xy_arrays!(x, y, idx_eh, eh)
    
    # E -> S:  0, 2, 5,  9, 18
    idx_es =  [1, 3, 6, 10, 19]
    es = LinSpace(easy, split, 5 )
    x, y = get_xy_arrays!(x, y, idx_es, es)

    # H -> S:  11--17
    idx_hs =  collect(11:19) 
    hs = LinSpace(hard, split, 9 )
    x, y = get_xy_arrays!(x, y, idx_hs, hs)    

    # Rest
    # P4 == 0.5( P3 + P5 )
    x[5] = 0.5 * ( x[4] + x[6])
    y[5] = 0.5 * ( y[4] + y[6])     

    x69 = LinSpace(x[7], x[10], 4)
    y69 = LinSpace(y[7], y[10], 4)
    
    x[8], y[8]  = x69[2], y69[2]
    x[9], y[9]  = x69[3], y69[3]        
    
    ΔEₚ = [ 0.0,
            3.2,
           11.5,
           19.2,
           17.6,
           39.9,
           31.1,
           29.9,
           38.7,
           75.2,
           39.3,
           37.0,
           34.8,
           35.3,
           37.9,
           48.4,
           60.7,
           76.9,
            108.9 ]


    # Maybe make it periodic, in x-y plane
    # > This involves reflection about the E-H and E-S lines and then reflection of the whole thing along 

    z = ΔEₚ

    triangle1 = zeros(n_points,3)
    triangle2 = zeros(n_points,3)
    triangle3 = zeros(n_points,3)
    triangle4 = zeros(n_points,3)
    triangle5 = zeros(n_points,3)
    triangle6 = zeros(n_points,3)

    
    original_data = zeros(19,3)
    original_data[:,1] = x
    original_data[:,2] = y
    original_data[:,3] = z    

    all_data= zeros( 19 * 6 * 4 * 2, 3 )
    all_data[1:19, 1:3] =  original_data

    all_data[1:n_points, 1:3] =  original_data
    triangle1  = original_data
    index = n_points + 1

    all_data     = reflect_data_about_line!(index, triangle1, all_data, hard)
    triangle2    = all_data[ index: index + n_points - 1, : ]
    index += n_points

    # all_data     = reflect_data_about_line!(index, triangle2, all_data, split, [0.,1])    
    # triangle3    = all_data[index : index + n_points - 1, : ]
    # index += n_points



    #all_data = reflect_data_about_line!(index, all_data[1:index-1,1:end], all_data, [0., 1.], [ √3 / 2, 0.] )

    all_data     = reflect_data_about_line!(index, triangle1, all_data, split)
    triangle4    = all_data[index : index + n_points - 1, : ]
    index += n_points

    all_data = reflect_data_about_line!(index, all_data[1:index-1,1:end], all_data, [0., 1.] )
    index += 3 * n_points

    # all_data     = reflect_data_about_line!(index, triangle4, all_data, hard, [0., -1])
    # triangle5    = all_data[index : index + n_points - 1, : ]
    # index += n_points

    # all_data     = reflect_data_about_line!(index, triangle5, all_data, split, [0., -1.])
    # triangle6    = all_data[index : index + n_points - 1, : ]
    # index += n_points
    
    # Now reflect all data across y axis
    all_data = reflect_data_about_line!(index, all_data[1:index-1,1:end], all_data, [0., 1.], [ √3 / 2, 0. ]  )

    # # tesselate the first hexagon
    # translated_hexagon = zeros(6*n_points, 3)
    # translated_hexagon[1:index-1, 1:end] = 
    # translated_hexagon[1:index-1, end] = all_data[1:index-1,end] 
    
    all_data[index + 6*n_points   : index + 6*n_points*2 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [  √3 / 2,  1.5, 0.0], (1,3) )
    all_data[index + 6*n_points*2 : index + 6*n_points*3 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [  √3 / 2, -1.5, 0.0], (1,3) )

    all_data[index + 6*n_points*3 : index + 6*n_points*4 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [ 3*√3 / 2,  1.5, 0.0], (1,3) )
    all_data[index + 6*n_points*4 : index + 6*n_points*5 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [ 3*√3 / 2, -1.5, 0.0], (1,3) )

    all_data[index + 6*n_points*5 : index + 6*n_points*6 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [ -√3 / 2,  1.5, 0.0], (1,3) )
    all_data[index + 6*n_points*6 : index + 6*n_points*7 - 1, 1:end] .= all_data[1:index-1,1:end] .+ reshape( [ -√3 / 2, -1.5, 0.0], (1,3) )
    
    # index += 6 * n_points
    # all_data = reflect_data_about_line!(index, all_data[1:index-1,1:end], all_data, [1., 0.], [ 0, 0.5 ]  )

    # Now mask to remove duplicates    
    mask = mask_duplicates(all_data)
    N_mask = sum(mask[:,1])
    new_data = zeros(Float64, (Int64(N_mask), 3) )
    new_data[:,:] .= reshape( all_data[mask], size(new_data) )


    # new_data[:,2] .+= 2
    # new_data[:,2] .= mod.(new_data[:,2], 3.0)

    
    # >>> Interpolating data by grid interpolation using Scipy Interpolate and griddata
    nx=150
    ny=150
    all_interp =  zeros( nx*ny, 3 )
    ndat = get_points_in_range(new_data, (0, √3), (-2, 1))
    x_grid = LinSpace(minimum(ndat[:,1]), maximum(ndat[:,1]),  nx)
    y_grid = LinSpace(minimum(ndat[:,2]), maximum(ndat[:,2]),  ny)
    X = kron(ones(nx), x_grid)
    Y = kron(y_grid, ones(ny))
    all_interp[:,1] .= X
    all_interp[:,2] .= Y    
    all_interp[:,3] .= si.griddata((ndat[:,1], ndat[:,2]), ndat[:,3], (X, Y) , method="cubic")

    
    # Using the interpolated function as input for Plane wave expansion
    # >> Plane wave expansion is not a good representation of the energy surface, therefore can just use the griddata
    # >> Get a function for the spline interpolation and then pass into the NEB / String Method

    # Create periodic area within some bounds
    # > Array with all the original points in is new_data

    ftd = reshape( all_interp[:,3], (nx,ny) )

    xr = (0, √3)
    yr = (-2, 1)    

    fexp, grad, hessian  = fft_expansion(ftd, xr, yr, nx, ny)

    # f = zeros(60*60)
    # dfx = zeros(60*60)
    # dfy = zeros(60*60)
    # x = zeros(60*60)    
    # y = zeros(60*60)
    
    # it = 0
    # for xi in LinSpace(-√3, 3*√3/2 , 60)
    #     for yi in LinSpace( -3, 3, 60 )
    #         it +=  1
    #         x[it] = xi 
    #         y[it] = yi 
    #         f[it] = fexp( xi  , yi )
    #         dfx[it] = dx_fexp( xi  , yi )
    #         dfy[it] = dy_fexp( xi  , yi )             
    #         #            println("it = $(it): x = $(x[it]), y = $(y[it]), ft = $(f[it])" )
    #     end
    # end

    # scatter3D( new_data[:,1], new_data[:,2], new_data[:,3])
    # surf( x, y, dfx, alpha=0.9, cmap = ColorMap("magma") )
    # PyPlot.show()

    return fexp, grad, hessian
    
end


end

