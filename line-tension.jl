# This is the run code for the dislocation line tension model. 
# One can use the Itakura or the tight binding peierls potentials


# addprocs()
# thisDir = dirname(@__FILE__())
# any(path -> path == thisDir, LOAD_PATH) || push!( LOAD_PATH, thisDir )

include("interaction_types.jl")
include("peierls_potential.jl")
include("trap_site_concentrations.jl")
include("dislocation_types.jl")
include("mcclean_isotherm_conc_dist.jl")


using PeierlsPotential
using DislocationTypes
using InteractionTypes
using McCleanIsotherm
using TrapSites
using StaticArrays

using Plots
using Dierckx

using Parameters
using SaddleSearch



function line_tension_model( N, potential="Normal", stress=zeros(3,3), interaction_type="C",
                             C_nom=433./1e6, ρ=1.e15, equilibrium=true, mcclean=true )
    # Define dislocation
    # > N = Number of images for the String method
    # > C_nom is the nominal carbon concentration
    # > ρ is the dislocation density
    # > redistribute switches the redistribution of the sites on or off
    # > mcclean switches the partial occupancies on, otherwise the position of the solute atoms remain fixed.

    global N_iter_total, N_images, name, scale

    file_ext = "$(name)_$(potential)_$(scale)"
    
    
    N_iter_total = 0
    N_images = N
    #    eV → meV  |  1/Å² → 1/b²
    conv = 1000 # * 0.75 * 2.87^2

    b_mag = 2.87 * √3 / 2 

    # fe, dfe, ddfe = get_2D_peierls_profile(pp_type=potential)

    #    interaction_type="C"
    if interaction_type == "C"
        interaction = C_Lorentzian{Float64}()
    elseif interaction_type == "H" 
        interaction = H_Lorentzian{Float64}()
    end
    

    if potential == "itakura"
        K = 0.866 * conv
    elseif potential == "sdTB"
        K = 0.734166 * conv
    elseif potential == "dTB"
        K = 0.734166 * conv
    end

    if mcclean
        # Get the splines (and it's derivative) which determine the concentration of a solute with distance
        C, dC = get_concentration_distance_dependence_splines(C_nom, ρ, interaction)
    else
        C = dC = x -> 0
    end

    interact = true

    if equilibrium
        conv_sitelabel = convert_sitelabel_to_pos_function()
        ref_conc_sum = get_reference_concentration(get_paths(zeros(2))..., conv_sitelabel, C)
        temp = 320.0 # K
        # Test
        # p_solute = ( √2 / 3 * 2.87) .* [√3/2 0.5 (46 * 3/2 * √(3/2))]' 
        # solutes_normal = solutes = Solutes{Float64}( interact, interaction, p_solute )
        # @show solutes_normal
        
        solutes = ConcSolutes{Float64}(interact, interaction, conv_sitelabel, C, ref_conc_sum, temp)
        @show solutes
    else
        p_solute = ( √2 / 3 * 2.87) .* [√3/2 0.5 (46 * 3/2 * √(3/2))]' 
        solutes = Solutes{Float64}( interact, interaction, p_solute )
    end

    
    b = SVector(0., 0., b_mag )

    d = Disl_line( Dislocation{Float64}(K, stress, b),
                   Potential{Function}( get_2D_peierls_profile(;pp_type=potential)... ),
                   solutes )


    eVÅ_to_GPa = 160.21766208
    println( "Resolved shear stress yz (MPa) yz: d.σ[2,3]  ", 1000 * eVÅ_to_GPa * d.d.σ[2,3])

    # Define initial and final positions
    ip1 = [ 0.0, 0.0 ]
    ip2 = [ √3 , 0.0 ]
    l1, l2, l3 = 30, 31, 30
    offset = 0
    

    # -->> Create Images <<--
    # > Create path of dislocations from one peierls valley to the next
    # > These are straight dislocations where the images are straight dislocations going over barrier.
    disl_string = get_kink_images(l1, l2, l3, offset, ip1, ip2, N)

    Ed  = x -> energy(d,x)
    # dEd = x -> ForwardDiff.gradient( y -> energy( d, y ) , x)
    
    dEd = x -> dE(d,x)

    tol = 0.1 
    tolP = 1e-3
    maxnit = 500_000
    verbose = 3


    # Find energy scale, μ, for the preconditioner
    # precon = x->[copy(precond(d, xn)) for xn in x]
    # x_test = disl_string[ceil(Int64,N/2)]
    # μ = get_mu_precond(x_test, dEd, precond(d, x_test))
    # precon = x->[copy(precond(d, xn; μ=μ)) for xn in x]


    # precon = x->[copy(hessprecond(d, xn; stab=0.1)) for xn in x]

    # precon = x -> [ copy( hessprecond(d, xn; stab=0.001) ) for xn in x ]
    #     precon = x->[copy(hessprecond(d,xn)) for xn in x]

    # preconP = SaddleSearch.localPrecon(precon = precon(disl_string),
    #                                    precon_prep! = (P, x) -> precon(disl_string))

    preconI = SaddleSearch.localPrecon(precon = [I], precon_prep! = (P, x) -> P)
    

    for img in 1:N
        mode="a"
        x = disl_string[ img ]
        write_object_to_file(reshape( x,  ceil(Int64, size(x,1)/2), 2 ),
                         "initial_image_positions_final_$(file_ext)_$img", "w")
        write_object_to_file(reshape( x,  ceil(Int64, size(x,1)/2), 2 ),
                         "initial_image_positions_final_$(file_ext)_ALL", mode)
        write_object_to_file(sum(energy(d,x)),
                         "initial_etot_final_$file_ext", mode)
    end


    tol = 1e-3
    path = ODEString(reltol=0.1, tol = tol, maxnit = maxnit,
                     precon_scheme = preconI,
                     path_traverse = serial(),
                     verbose = verbose,
                     fixed_ends = false)

    # Fixed ends is false here, but constraining the gradient to be zero if at the ends, as in SaddleSearch/src/string.jl. 
    #       dE0_temp = [[zeros(x[1])]; [dE(x[i]) for i in direction[2:end-1]];
    #                  [zeros(x[1])]]
    #    cost = length(param) - 2
    # NOT reducing the cost parameter, will see if it affects anything


    X =  Path(disl_string)
    write_line_energies_of_images(d, disl_string, "initial")
    Xout, PATHlog, _ = SaddleSearch.run!(path, Ed, dEd, X )

    write_images_only(d,Xout)    
    write_line_energies_of_images(d,Xout)

    
end

function dE(V,x)
    global N_iter_total, N_images

    img = N_iter_total % ( N_images ) + 1
    write_grad = (img == 1 || img == N_images) ? false : true

    write_images(V,x)
    
    return DislocationTypes.construct_gradient( V, x, write_grad )# DislocationTypes.gradient( V, x )
end


function get_mu_precond(x, g, P, rnn=1.0)
    M = 1e-2 * rnn * UniformScaling(1.0)
    Lx, Ly = √3, 1.
    X = reshape( x,  ceil(Int64, size(x,1)/2),2 )'
    ν = M * vcat( sin.( X[1,:]/Lx ), sin.( X[2,:]/Ly ) )
    println("get_mu_precond: g(x+ν) = ", g(x+ν))
    println("get_mu_precond: g(x)   = ", g(x))
    println("get_mu_precond: ν      = ", ν)
    println("get_mu_precond: ν'Pν   = ", ν'*P*ν)
    μ = ν' * ( g(x+ν) - g(x) ) / ( ν'*P*ν )

    println("get_mu_precond: μ = ", μ)
    return μ
end


function exp_precond(d::Disl_line, X::Matrix; μ=1.0, rcut = 3.5, α=3.0)
    nX = size(X, 2)
    P = zeros(2*nX, 2*nX)
    I = zeros(Int, 2, nX)
    I[:] = 1:2*nX
    
    # Looking at the paper by Packwood (2016), they use A = 3.0, μ=1
    A = 3
    #    μ = 1.0 #/1000
    rnn = 1

    conv = √2 / 3 * 2.87
    b = √3 / 2 * 2.87
    Rij = zeros(3)
    H(Xi)  = ForwardDiff.hessian(  y -> d.ΔEₚ(y...), Xi  )
    #    Hl(Xi) = ForwardDiff.hessian(  y -> DislocationTypes.get_interaction_energy(d, j, Pⱼ, false), Xi  )
    for i = 1:nX, j = i+1:nX

        Ii, Ij = I[:,i], I[:,j]
        # Rij[1:2] .= ( X[:,i] - X[:,j] ) .* conv ./ b
        # Rij[3] = ( j - i )

        rij = abs( j - i )
        
        
        R = norm(Rij)
        if 0 < rij < rcut

            #            a = pphi(rij) * eye(2)
            a = μ*exp( -A*( rij/rnn - 1. ) ) * eye(2)
            # Can change the preconditioning thing to have relevance to peierls potential
            # DislocationTypes.approx_hessian(d, x) 
            #       a = μ * ( H(X[:,i]) .* eye(2) +  d.K * eye(2) + Hl( norm( X[:,i] - d. ) .* conv  )[1]*eye(2)  ) .* exp( -A*( R/rnn - 1. )) # / 1000. # .* exp( -A*( rij/rnn - 1. ) ) #  μ * exp( -A*( rij/rnn - 1. ) ) * eye(2)
            # a = μ * exp( -A*( rij/rnn - 1. ) ) * eye(2)
            P[Ii, Ij] -= a
            P[Ij, Ii] -= a
            P[Ii, Ii] += a
            P[Ij, Ij] += a

        end
    end
    return sparse(P) + μ * 0.1 * speye(2*nX)
end


precond(d, x; μ=1.0) = exp_precond(d, reshape( x,  ceil(Int64, size(x,1)/2), 2 )'; μ = 1.0 )


function hessprecond(d, x; stab=0.0)
   H = Symmetric( DislocationTypes.construct_hessian(d, x))
   D, V = eig(H)
   D = abs.(D) .+ stab
   return V * diagm(D) * V'
end


function write_line_energies_of_images(d,X, prefix="")
    for (i,x) in enumerate(X)
        write_line_energy_of_image(d, x, i, prefix)
    end
end

function write_line_energy_of_image(d, xi, image, prefix)
    global name, potential, scale
    file_ext = prefix*"$(name)_$(potential)_$(scale)"

    N = ceil(Int64, size(xi,1)/2)
    storage = zeros(eltype(xi), N, 3)
    storage[:,1:2] .= reshape( xi, N, 2 )
    storage[:,3  ] .= [ DislocationTypes.energy_point(d, xi, j, N) for j in 1:N ]
    
    write_object_to_file( storage, "line_energies_image_$(image)_$(file_ext)", "w")
end



function write_images_only(V,X)
    global  name, potential, scale
    file_ext = "$(name)_$(potential)_$(scale)"
    #   println(X)
    #     println("Xsize = ", size(X))
    img_count = 0
    for i in 1:length(X)
        x = X[i]
        #        println("xsize = ", size(x))
        #        println(x)

        (img_count == 0? mode = "w" : mode = "a")
        write_object_to_file(reshape( x,  ceil(Int64, size(x,1)/2), 2 ),
                             "image_positions_final_$(file_ext)_$img_count", "w")
        write_object_to_file(reshape( x,  ceil(Int64, size(x,1)/2), 2 ),
                             "image_positions_final_$(file_ext)_ALL", mode)
        
        write_object_to_file(DislocationTypes.energy(V,x),
                             "etot_final_$file_ext", mode)
        img_count += 1
    end

end


function write_images(V,x)
    global N_iter_total, N_images, name, potential, scale
    file_ext = "$(name)_$(potential)_$(scale)"

    img_count = N_iter_total % ( N_images  ) 

    (img_count == 0? mode = "w" : mode = "a")
   
    # Trap positions occupancy
    if N_iter_total > 0 && img_count == 0
        if N_iter_total == N_images
            cp("trap_positions_occupancy", "trap_positions_occupancy_initial"; remove_destination=true)
        else
            cp("trap_positions_occupancy", "trap_positions_occupancy_final"; remove_destination=true)
        end
    end

    open( "trap_positions_occupancy",  mode) do io 
        println(io, "Image $(img_count + 1)")
    end
  

    write_object_to_file(reshape( x,  ceil(Int64, size(x,1)/2), 2 ),
                         "image_positions_final_$(file_ext)_$img_count", "w")
    write_object_to_file(reshape( x,  ceil(Int64, size(x,1)/2), 2 ),
                         "image_positions_final_$(file_ext)_ALL", mode)
    
    write_object_to_file(energy(V,x),
                         "etot_final_$file_ext", mode)

    N_iter_total += 1

end


function write_object_to_file(object, filename,  mode)
    open( filename,  mode) do io
        writedlm( io,  object)
    end
end




function get_kink_images(l1, l2, l3, offset, ip1, ip2, N; simple=false)


    disl0=zeros(l1 + l2 + l3,2)
    disl0[:,1] .= ip1[1]
    disl0[:,2] .= ip1[2]

    disl1=zeros(l1+l2+l3,2)
    disl1[:,1] .= ip2[1]
    disl1[:,2] .= ip2[2]


    dislm =zeros(l1+l2+l3,2)

    offset = 0

    l1s = ceil(Int64, l1/2)
    lspace  = linspace( 0, 1, l1s)
    nlspace = linspace( 1, 0, l1s)
    
    for i in 1:size(dislm,1)
        if i <= l1s
            dislm[i,1] .= ip1[1];  dislm[i,2] .= ip1[2]
        end

        if (i > l1s)
            if i <= l1
                inti = (i - l1s)
                np = nlspace[inti] .* ip1 + lspace[inti] .* ip2
                dislm[i,1] .= np[1];  dislm[i,2] .= np[2]
            else
                if i <= l1 + l2
                    dislm[i,1] .= ip2[1];  dislm[i,2] .= ip2[2]
                end
            end
        end
            
        
        if (i > l1+l2 )
            if i <= l1+l2+l1s
                inti = (i - (l1+l2) )
                np = nlspace[inti] .* ip2 + lspace[inti] .* ip1
                dislm[i,1] .= np[1];  dislm[i,2] .= np[2]
            else 
                dislm[i,1] .= ip1[1];  dislm[i,2] .= ip1[2]
            end
                
        end
    end


    # 120 length vectors
    d0 = vcat(disl0...) # SVector{2*(l1+l2+l3)}( vcat(disl0...) )
    d1 = vcat(disl1...) # SVector{2*(l1+l2+l3)}( vcat(disl1...) )
    dm = vcat(dislm...) # SVector{2*(l1+l2+l3)}( vcat(dislm...) )
    
    disl_string_sm =  [ (1-s) * d0  +  (s * dm)  for s in linspace(0, 1, ceil(Int64,N/2)) ]
    disl_string_mf =  [ (1-s) * dm  +  (s * d1)  for s in linspace(0, 1, ceil(Int64,N/2)) ]
    
    disl_string = [ i < ceil(Int64,N/2)? disl_string_sm[i] :  disl_string_mf[i-floor(Int64,N/2)] for i in 1:N]

    
    #        disl_string = [ (1-s) * d0  +  (s * d1)  for s in linspace(0, 1, N) ]
    
    return disl_string
end


function read_lt_points(file)
    points = readdlm( file, '\t', Float64, '\n' )

    println("\n From file $file \n", points)
    return points
end

function get_points_energy( d, file )
    points = read_lt_points(file)
    E_line = energy_function_of_positions( d, vcat(points...), d.ΔEₚ )

    E_points = zeros( length(E_line), 3 )

    E_points[:,1:2] .= points
    E_points[:,3] .= E_line
    
    write_object_to_file( E_points, "line_energies_$file", "w")
end

function get_line_energies(file, potential)
    conv = 1000 # * 0.75 * 2.87^2

    d = Disl_line( ΔEₚ   = get_2D_peierls_profile(pp_type = potential)[1],
                   ∂ΔEₚ  = get_2D_peierls_profile(pp_type = potential)[2],
                   hess  = get_2D_peierls_profile(pp_type = potential)[3],
                       K = 0.866 * conv,
                       σ = zeros(3,3),
                       b = [0., 0., 1.] )
    
   get_points_energy( d, file ) 
    
end

function get_all_line_energies()
    
    files = ["path_for_results_itakura_90pts_ODEstring_Nimg_15_Econv",
             "path_for_results_sdTB_90pts_ODEstring_Nimg_15_Econv",
             "path_for_results_dTB_90pts_ODEstring_Nimg_15_Econv" ]
    potentials = [ "itakura", "sdTB", "dTB"]

    for (file, potential) in zip(files, potentials)
        get_line_energies(file, potential)
    end

end


function rotate_stress_tensor( σ )
    # Rotate the stress tensor by the simple rotation rules
    
    R =  transpose( [-2.       0.     √2  ;        # {After transpose} maps  < -2  1  1 > -> x
                     1.     -√3      √2  ;         # {After transpose} maps  <  0 -1  1 > -> y
                     1.      √3      √2  ] ./ √6 ) # {After transpose} maps  <  1  1  1 > -> z

    # This is the transformation which changes 
    # σ' = R⋅σ⋅Rᵗ

    σ_new = R * σ * R'
    println("Stress tensor\n", σ, "\nRotated Stress tensor\n",σ_new)
    return σ_new
    
end



# get_all_line_energies()

# The stress tensor should be scaled by the shear modulus, which
# in this model is μ = 50 GPa

global name, potential, scale
# Convert the stress into other units
# 1 eV/Å³ = 160.21766208 GPa
GPa_to_eVÅ³ = 1/160.21766208

μ = 1. * GPa_to_eVÅ³  # GPa

# println(ARGS)

if length(ARGS) != 0
    N = parse(Int, ARGS[1])
    potential = ARGS[2]
    scale = parse(Float64,ARGS[3])
    name = ARGS[4]
    C_nom = parse(Float64,ARGS[6])/1e6
else
    N = 35
    potential = "dTB"
    scale = 1.0
    name = "test_hessian"
    eq = true
end
# Stress in yz is what we want for this coordinate system
stress = (- scale * μ ) .* [ 0  0  0 ;
                             0  0  1 ;
                             0  1  0. ]

interaction_type="C"
# C_nom=433./1e6
ρ=1.e15
equilibrium=true
mcclean=true

if potential == "sdTB"
    a = 126.62
    b = -51.854
elseif potential == "dTB"
    a = 144.666
    b = -63.168
elseif potential == "itakura"
    a = 166.0
    b = -46.67
else
    println("\n WARNING: Coefficients for stress matrix aren't specified.\n")
end

a =  1.
b = -1.
# Stress in yz is what we want for this coordinate system
stress = SMatrix{3,3}( (- scale * μ ) .* [ 0  b  0 ;
                             b  0  a ;
                             0  a  0. ] )

# stress = (- scale * μ ) .* [ 0  0  0 ;
#                              0  0  1 ;
#                              0  1  0. ]


# Itakura does 0 → 1 GPa for the shear stress. 
# so go in increments of 0.01 for 100 points

#  rotate_stress_tensor( stress )

line_tension_model(N, potential, stress, interaction_type, C_nom, ρ, equilibrium, mcclean )
