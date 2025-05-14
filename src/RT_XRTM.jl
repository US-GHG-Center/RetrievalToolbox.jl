"""
$(TYPEDSIGNATURES)

Determines which weighting functions we must ask XRTM to calculate for us, based on the
state vector. Results are stored in a dictionary which contains `String` keys that map to
`Vector{Int}`. The vectors of integers will then correspond to the derivative indices used
by XRTM. The dictionary used here is the `wfunctions_map` field of the
`MonochromaticRTMethod` structure, and is emptied at the start of this function.

# Details

This function interrogates the state vector `sv`, and, depending on the state vector
elements found, will create respective entries in a Dict. Generally, we ask XRTM for these
basic types of weighting functions, which then are later transformed into Jacobians.

1) ∂I/∂τ for each layer
2) ∂I/∂ω for each layer
3) ∂I/∂s (s being the surface kernel amplitude) for each kernel

The dictionary to be returned contains the vectors of integers which determine the
position of the appropriate XRTM weighting functions for the various derivatives. For
example, `d[dI_dTau] = [6,7,8,9]` would mean that the per-layer ∂I/∂τ (of a 4-layer
atmosphere) are stored in the XRTM weighting function indices 6 through 9.

Notable exceptions are aerosol-related Jacobians. Computing per-layer Jacobians for every
aerosol mixture is computationally inefficient. Instead, we compute dedicated weighting
functions for each and every aerosol-related state vector element.

"""
function _create_weighting_function_dictonary!(
    d::Dict{Any, Vector{Int}},
    rt::MonochromaticRTMethod
    )

    # Get handy short-cuts
    scene = rt.scene
    swin = rt.optical_properties.spectral_window
    sv = rt.state_vector

    # Empty out the dictionary!
    empty!(d)

    # Start counting from 1
    current_idx = 1

    #=
        Do we need per-layer ∂I/∂τ, and ∂I/∂ω?
    =#
    need_dI_dTau = false
    need_dI_dOmega = false
    need_dI_dBeta = false

    if any_SVE_is_type(sv, GasLevelScalingFactorSVE)
        @debug "Found GasLevelScalingFactorSVE! Need ∂I/∂τ"
        need_dI_dTau = true
    elseif any_SVE_is_type(sv, GasVMRProfileSVE)
        @debug "Found GasVMRProfileSVE! Need ∂I/∂τ"
        need_dI_dTau = true
    elseif any_SVE_is_type(sv, TemperatureOffsetSVE)
        @debug "Found TemperatureOffsetSVE! Need ∂I/∂τ"
        need_dI_dTau = true
    elseif any_SVE_is_type(sv, SurfacePressureSVE)
        @debug "Found SurfacePressureSVE! Need ∂I/∂τ"
        need_dI_dTau = true
    end

    #=
        ∂I/∂ω is required when there is either Rayleigh scattering or
        aerosols in the atmosphere
    =#
    if findanytype(scene.atmosphere.atm_elements, AbstractRayleighScattering)
        @debug "Found RayleighScattering! Need ∂I/∂ω"
        need_dI_dOmega = true
    end
    if findanytype(scene.atmosphere.atm_elements, AbstractAerosolType)
        @debug "Found some aerosol! Need ∂I/∂ω"
        need_dI_dOmega = true
        need_dI_dBeta = true
    end

    if need_dI_dTau
        # Need one weighting function per atmospheric layer
        d["dI_dTau"] = collect(current_idx:current_idx+scene.atmosphere.N_layer-1)
        current_idx += scene.atmosphere.N_layer
    end

    if need_dI_dOmega
        # Need one weighting function per atmospheric layer
        d["dI_dOmega"] = collect(current_idx:current_idx+scene.atmosphere.N_layer-1)
        current_idx += scene.atmosphere.N_layer
    end

    #=
        Do we need surfaces?
    =#

    if any_SVE_is_type(sv, BRDFPolynomialSVE)
        # For surface polynomials, we only require ∂I/∂s, and can
        # construct the polynomial coefficient derivatives via chain rule
        # ∂I/∂c_i = ∂I/∂s * ∂s/∂c_i
        @debug "Found surface retrieval! Need ∂I/∂s"
        d["dI_dSurface"] = Int[]

        # We need one index per BRDF kernel..
        for kernel in get_surface(rt.scene, swin).kernels
            push!(d["dI_dSurface"], current_idx)
            current_idx += 1
        end
    end


    #=
        Aerosols (only retrieved ones)
    =#

    for sve in sv.state_vector_elements
        if is_aerosol_SVE(sve)
            if !haskey(d, sve.aerosol)

                d[sve.aerosol] = collect(current_idx:current_idx+scene.atmosphere.N_layer-1)
                current_idx += scene.atmosphere.N_layer
            end
        end
    end

    @debug "Require a total of $(current_idx - 1) XRTM weighting functions"

end

function _calculate_radiances_and_jacobians_XRTM!(
    rt::MonochromaticRTMethod,
    observer::Union{SatelliteObserver, UplookingGroundObserver};
    xrtm_in::Union{Nothing, Ptr{Nothing}}=nothing
    )

    #=
        Two-step process for XRTM. First, we calculate radiances and the weighting
        functions (partial derivatives w.r.t. some elemental linearized inputs). Then we
        call the `calculate_rt_jacobian!` functions which turn those so-called weighting
        functions into proper Jacobians to represent partial derivatives w.r.t. state
        vector elements.
    =#


    # If the user supplies one dictionary, we only need to run XRTM with a single
    # configuration.
    if rt.model_options isa AbstractDict
        _calculate_radiances_and_wfs_XRTM!(
            rt, observer, rt.model_options, xrtm_in=xrtm_in)
    end

    # If the user supplies a Vector{} of Dictionaries, simply loop over all
    # options and call the routine for every one of them.
    if rt.model_options isa Vector{<:AbstractDict}

        for model_option in rt.model_options
            _calculate_radiances_and_wfs_XRTM!(
                rt, observer, model_option, xrtm_in=xrtm_in)
        end

    end

    #=
        Calculate Jacobians once the radiances and weighting functions are done.
        (no XRTM needed here, all results should be accessible from within `rt`)

        Reminder: functions `calculate_rt_jacobian!` for some SVE usually turn
        weighting functions into an actual Jacobian, but (usually) need the full
        radiance field calculations to be done. It is hence wise to execute this function
        after all the RT business is finished.

        Note that when using RT acceleration methods such as LSI, these functions should
        be called again! Otherwise the Jacobians will still be those from the initial RT
        calculations and not contain any corrections applied to them!
    =#

    for (i, sve) in enumerate(rt.state_vector.state_vector_elements)
        if calculate_jacobian_before_isrf(sve)
            calculate_rt_jacobian!(rt.hires_jacobians[sve], rt, sve)
        end
    end

end


function create_XRTM(
    rt::MonochromaticRTMethod,
    observer::Union{SatelliteObserver, UplookingGroundObserver},
    options_dict::AbstractDict,
)

    #=
        The XRTM handler has to be created with some fixed sizes and parameters,
        which we derive here from the RT object
    =#


    options = options_dict["options"]
    solvers = options_dict["solvers"]


    #=
        If derivatives are needed/requested, we must keep track of which
        derivatives are calculated at which location.
    =#

    _create_weighting_function_dictonary!(rt.wfunctions_map, rt)
    if "calc_derivs" in options
        n_derivs = map(length, values(rt.wfunctions_map)) |> sum
    else
        n_derivs = 0
    end
    @debug "Requiring $(n_derivs) derivatives from XRTM"



    # Set max coefs to either 1 (no aerosols in scene), or whatever
    # the optical properties tell us to use
    if isnothing(rt.optical_properties.total_coef)
        max_coef = 1
    else
        max_coef = size(rt.optical_properties.total_coef, 1)
    end

    # XRTM wants half-space streams, but we are more used to
    # full-space notation..
    if !(options_dict["streams"] isa Integer)
        @error "XRTM streams needs to be an integer!"
    end
    if (options_dict["streams"] < 2)
        @error "XRTM streams must be > 1!"
    end

    n_quad = options_dict["streams"] // 2

    #=
        We peform vector calculations ONLY if the XRTM option "vector" is passed,
        regardless whether the underlying radiance container is a ScalarRadiance or a
        VectorRadiance.
        The reasoning is as follows. For some applications, we want to calculate the
        scalar multiple scattering radiance and add the contributions to the total vector
        radiance.
    =#
    if !("vector" in options)
        n_stokes = 1
    elseif ("vector" in options) & (rt.hires_radiance isa VectorRadiance)
        n_stokes = 3
    else
        @error "Incompatible XRTM option: if `vector` is set, you must also supply a
                VectorRadiance type."
        return false
    end

    n_layers = rt.scene.atmosphere.N_layer
    n_theta_0s = 1 # Incoming solar zenith angles
    n_kernel_quad = 16

    #=
        NOTE that XRTM requires a surface, regardless if you want to
        do an uplooking scenario or not
    =#

    kernels = String[]

    # Find out, which XRTM surface kernels we need

    swin = rt.optical_properties.spectral_window
    if swin isa BinnedSpectralWindow
        swin = swin.original_window
    end

    if rt.scene.surfaces[swin] isa BRDFSurface
        # A BRDF surface is a vector of BRDF combinations. Note that we keep the ordering
        # as per the order in the vector itself. So the XRTM BRDF kernel `k-1` will be the
        # k'th entry here.
        for kernel in get_surface(rt.scene, swin).kernels
            if kernel isa LambertianPolynomialKernel
                push!(kernels, "lambertian")
            elseif kernel isa RPVPolynomialKernel
                push!(kernels, "rahman")
            end
        end
    else
        @error "XRTM is currently only set up to use `BRDFSurface`-type surfaces!"
        return false
    end



    n_out_levels  = 1 # Number of output levels
    n_out_thetas  = 1 # Number of output viewing zenith angles
    #n_out_phis    = 1 # Number of output viewing azimuth angles


    #=
        Create the XRTM handler object with options and parameters as decided above.
    =#

        xrtm = XRTM.create(
            options,
            solvers,
            max_coef,
            n_quad,
            n_stokes,
            n_derivs,
            n_layers,
            n_theta_0s,
            n_kernel_quad,
            kernels,
            n_out_levels,
            n_out_thetas
            )

    return xrtm

end

"""
$(TYPEDSIGNATURES)

# Details

## Notes on threading

The function was written to support multi-threading, so running julia via e.g.
`julia -t N` to spawn `N` threads can be used to speed up the overall calculation by
spreading the RT calculations over those threads. Rather than creating one XRTM object,
we create `N` objects, and each thread only uses the one XRTM instance that it is assigned
to via its `Threads.threadid()`.

"""
function _calculate_radiances_and_wfs_XRTM!(
    rt::MonochromaticRTMethod,
    observer::Union{SatelliteObserver, UplookingGroundObserver},
    options_dict::AbstractDict;
    xrtm_in::Union{Nothing, Ptr{Nothing}, Vector{Ptr{Nothing}}}=nothing
    )

    atm = rt.scene.atmosphere
    solvers = options_dict["solvers"]
    options = options_dict["options"]

    if isnothing(xrtm_in)
        # Create a new XRTM instance, one per thread
        xrtm_l = [create_XRTM(rt, observer, options_dict)
            for i in 1:Threads.nthreads()]

    else
        # Use the user-supplied XRTM instance instead
        if xrtm_in isa Vector
            xrtm_l = xrtm_in
        else
            xrtm_l = [xrtm_in]
        end
    end

    #=
        Set necessary parameters
    =#

    if "sos" in solvers
        # Successive orders of scattering:
        # orders, maximum tau, radiance minimum tolerance
        for xrtm in xrtm_l
            XRTM.set_sos_params(xrtm, 2, 10.0, 1e10)
        end
    end

    if "pade_add" in solvers
        if haskey(options_dict, "pade_add")

            pade_s = options_dict["pade_add"][1]
            pade_r = options_dict["pade_add"][2]
            @debug "Setting Padé parameters: s = $(pade_s), r = $(pade_r)"
            for xrtm in xrtm_l
                XRTM.set_pade_params(xrtm, pade_s, pade_r)
            end
        else
            @debug "Padé parameters not supplied via `pade_add`. " *
            "XRTM will use its own lookup table."
        end
    end

    # Set the output level(s)

    if observer isa SatelliteObserver
        for xrtm in xrtm_l
            XRTM.set_out_levels(xrtm, Int32[0])
        end
    elseif observer isa UplookingGroundObserver
        # Output at the last level = BOA level = bottom of BOA layer
        # (level indices for XRTM are 0 .. n_layer)
        for xrtm in xrtm_l
            XRTM.set_out_levels(xrtm, Int32[n_layers])
        end
    end

    # Set the output zenith angles
    out_thetas = Float64[]
    push!(out_thetas, rt.scene.observer.viewing_zenith)
    for xrtm in xrtm_l

        XRTM.set_fourier_tol(xrtm, 1e-4)

        XRTM.set_out_thetas(xrtm, out_thetas)

        # Set the isotropic TOA and BOA contributions
        # (F_iso_bot is where e.g., SIF would go)
        XRTM.set_F_iso_top(xrtm, 0.)
        XRTM.set_F_iso_bot(xrtm, 0.)

        # Set solar zenith
        XRTM.set_theta_0(xrtm, rt.scene.solar_zenith)

        # Set solar azimuth
        XRTM.set_phi_0(xrtm, 0.0) #rt.scene.solar_azimuth)
    end

    # Pseudo-spherical approximation can be used on request
    if "psa" in options
        # Set plantary radius and altitude levels
        p_radius = ustrip(rt.scene.atmosphere.altitude_unit, 6371.0u"km")

        for xrtm in xrtm_l
            XRTM.set_planet_r(xrtm, p_radius)
        end

        #=
            Note that the atmosphere object contains the altitude/gravity for the
            meteorological grid only, but not the retrieval/RT grid.  We assume the
            MET grid values are smooth enough such that we can just pick linearly
            interpolated values.
        =#

        altitude_int = linear_interpolation(
            ustrip.(Ref(atm.pressure_unit), atm.met_pressure_levels * atm.met_pressure_unit),
            atm.altitude_levels,
            extrapolation_bc = Line()
            )


        altitude_levels = altitude_int.(atm.pressure_levels)
        for xrtm in xrtm_l
            XRTM.set_levels_z(xrtm, altitude_levels)
        end
    end


    # Set derivatives, if requested. These are simply set once since they stay the same
    # for all spectral points.
    if "calc_derivs" in options
        if haskey(rt.wfunctions_map, "dI_dTau")
            # Set one layer derivative per derivative index
            # NOTE -- this code here assumes that "dI_dTau" really means
            # one derivative per layer, and that they are ordered
            for (lay, idx) in enumerate(rt.wfunctions_map["dI_dTau"])
                # XRTM counts 0-based, e.g. TOA layer has index 0!
                for xrtm in xrtm_l
                    XRTM.set_ltau_l_11(xrtm, lay-1, idx-1, 1.0)
                end
            end
        end

        if haskey(rt.wfunctions_map, "dI_dOmega")
            # Set one layer derivative per derivative index
            # NOTE -- this code here assumes that "dI_dOmega" really means
            # one derivative per layer, and that they are ordered
            for (lay, idx) in enumerate(rt.wfunctions_map["dI_dOmega"])
                # XRTM counts 0-based, e.g. TOA layer has index 0!
                for xrtm in xrtm_l
                    XRTM.set_omega_l_11(xrtm, lay-1, idx-1, 1.0)
                end
            end
        end

        if haskey(rt.wfunctions_map, "dI_dSurface")
            for (k, idx) in enumerate(rt.wfunctions_map["dI_dSurface"])
                for xrtm in xrtm_l
                    XRTM.set_kernel_ampfac_l_1(xrtm, k-1, idx-1, 1.0)
                end
            end
        end
    end

    # Other derivative calculations are given by dedicated functions inside the
    # spectral loop (see _run_XRTM!)

    #=
        Perform XRTM calculations
    =#

    _run_XRTM!(
        xrtm_l,
        rt,
        options_dict
    )

    #=
        Destroy the XRTM handler object and free memory!

        But only if there is no user-supplied one. If the user supplies their own XRTM
        to this function, they are responsible for destroying it!
    =#

    if isnothing(xrtm_in)
        for xrtm in xrtm_l
            XRTM.destroy(xrtm)
        end
    end

end


"""
$(TYPEDSIGNATURES)


"""
function _run_XRTM!(
    xrtm_l::Vector{Ptr{Nothing}},
    rt::AbstractRTMethod,
    options_dict::AbstractDict
)

    N_layer = rt.scene.atmosphere.N_layer
    N_hires = rt.optical_properties.spectral_window.N_hires

    # Grab the number of derivatives calculated by XRTM
    n_derivs = XRTM.get_n_derivs(xrtm_l[1])
    # Grab the number of Stokes elements
    n_stokes = XRTM.get_n_stokes(xrtm_l[1])

    if n_stokes == 1
        n_elem = 1
    elseif n_stokes == 3
        n_elem = 6
    end

    # Due to an incosistency in the coordinate conventions for single scattering and
    # 2OS solvers, the sign of the Stokes U components must be flipped in cases where
    # the azimuth is adjusted to stay inside of (0, 360).
    flip_U = false

    # Set the total phase function expansion coefficients for all layers
    # NOTE!
    # XRTM wants to be supplied with the *full* coef matrix, which needs to have
    # max_coef coefficients; even if you plan to only use a smaller number..
    if isnothing(rt.optical_properties.total_coef)

        n_coef = Int32(1)
        n_coef_arr = [n_coef for i in 1:N_layer]

        # XRTM requires an array, even if we do not have any aerosols or Rayleigh.
        # So we have to set a [1,1,N_layer] array, for example, filled with 1.0 for
        # fully isotropic scattering, while at the same time (hopefully) setting the
        # single-scatter albedo to 0, so the scattering contribution will be zero.
        coef_arr = ones(1, n_elem, N_layer)
        for xrtm in xrtm_l
            XRTM.set_coef_n(xrtm, n_coef_arr, coef_arr)
        end

    else

        #=
            Determine how many phase function coefficients should be used in each layer.
            This has a significant impact on performance.

            General consensus is that for single-scattering, one should use all available
            coefficients, and then use 2*n_quad - 1 coefficients for multiple scattering,
            and apply the Delta-M correction (and use TMS correction for single sc.).

            Further, it is generally not necessary to use all coefficients in every layer.
            For a given layer `l`, most of the scattering might be due to one aerosol
            mixture, which is characterized by a low number of coefficients. In that case
            it would be a convenient approximation to reduce the number of considered
            coefficients, rather than calculate the contributions of all coefficients from
            all aerosols which might have several thousand coefficients.

        =#
        n_coef_arr = _calculate_needed_n_coef(rt)

        # Users might want to override the number of coefs used
        if haskey(options_dict, "n_coefs")
            for i in 1:length(n_coef_arr)
                n_coef_arr[i] = min(options_dict["n_coefs"], n_coef_arr[i])
                #@info "Setting ncoef for $(i): $(n_coef_arr[i])"
            end
        end

        for xrtm in xrtm_l
            XRTM.set_coef_n(
                xrtm,
                n_coef_arr,
                rt.optical_properties.total_coef[:,1:n_elem,:]
                )
        end

        # Coef derivatives have to be set AFTER the coefficients are ingested
        # (for the first time at least)

        if "calc_derivs" in options_dict["options"]

            # The number of coefs might be different for every layer, so we must
            # create a new derivative coef matrix for each layer, sized correctly.
            l_coef_layer = Matrix{Float64}[]
            for l in 1:N_layer
                n_coef = XRTM.get_n_coef(xrtm_l[1], l-1)
                l_coef = zeros(n_coef, n_elem)
                push!(l_coef_layer, l_coef)
            end

            # Loop through aerosol-related SVEs
            for sve in filter(is_aerosol_SVE, rt.state_vector.state_vector_elements)

                # Grab the aerosol from the state vector element
                aer = sve.aerosol

                if !haskey(rt.wfunctions_map, aer)
                    @debug "No weighting function assigment for aerosol $(aer)."
                    continue
                end

                for (l, wf_idx) in enumerate(rt.wfunctions_map[aer])

                    create_aerosol_coef_deriv_inputs!(
                        l_coef_layer[l],
                        rt,
                        aer,
                        l
                    )

                    l_coef_layer[l][1,1] = 0 # This must be zero always

                    for xrtm in xrtm_l
                        # Ingest into each XRTM instance.

                        XRTM.set_coef_l_11(
                            xrtm,
                            l-1,
                            wf_idx-1,
                            l_coef_layer[l]
                        )

                    end # End XRTM list loop
                end # Layer loop

            end # Aerosol loop
        end # End "do we need Jacobians?"


    end

    # Output azimuth angles: [N_output_azimuths, N_output_zeniths]
    out_phis = Matrix{Float64}(undef,1,1)
    if rt.scene.observer isa SatelliteObserver

        vaa = rt.scene.observer.viewing_azimuth + 180
        saa = rt.scene.solar_azimuth

        # Set ϕ to be the relative azimuth
        out_phis[1,1] = vaa - saa

        # Adjust ϕ such that ϕ - ϕ0 is > 0 but < 360
        # (phi_0 is set to 0)
        if (vaa - saa) > 360
            @debug "Subtracting 360 from azimuth!!"
            out_phis[1,1] -= 360.0
            flip_U = true
        end
        if (rt.scene.observer.viewing_azimuth - rt.scene.solar_azimuth) < 0
            @debug "Adding 360 to viewing azimuth!!"
            out_phis[1,1] += 360
            flip_U = true
        end

    elseif rt.scene.observer isa UplookingGroundObserver
        out_phis[1,1] = 0.0
    end

    # Zero out all radiance containers, unless user declares otherwise
    if haskey(options_dict, "add")
        if options_dict["add"] == true
            @debug "Model option -add- found and set to -false-: " *
                "we are *NOT* zero-ing out radiances and " *
                "derivatives, but adding to previous results!"
        else
            rt.hires_radiance[:] .= 0
            for i in 1:n_derivs
                rt.hires_wfunctions[i][:] .= 0
            end
        end
    else
        rt.hires_radiance[:] .= 0
        for i in 1:n_derivs
            rt.hires_wfunctions[i][:] .= 0
        end
    end

    # Take the used spectral window ...
    swin = rt.optical_properties.spectral_window
    # Get the surface attached to it
    surface = get_surface(rt.scene, swin)
    #=
        BRDF kernel parameters can be set here, outside the spectral loop. Only the
        amplitudes are spectrally varying (for now), so we can set the spectrally constant
        kernel parameters in advance, without having to repeatedly call the relevant
        setter function.
    =#

    for (k, kernel) in enumerate(surface.kernels)

        if kernel isa RPVPolynomialKernel
            #=
                Rahman (or RPV) kernel parameters are:
                    0: hotspot parameter
                    1: asymmetry
                    2: anisotropy
            =#

            for xrtm in xrtm_l
                XRTM.set_kernel_params_1(xrtm, k-1, 0, kernel.hotspot)
                XRTM.set_kernel_params_1(xrtm, k-1, 1, kernel.asymmetry)
                XRTM.set_kernel_params_1(xrtm, k-1, 2, kernel.anisotropy)
            end
        end

    end

    #=
        Establish which spectrally dependent XRTM weighting functions we must use!

        Some weighting functions, notably aerosol-related for example, require their
        linearized inputs to be calculated for every spectral point. In order to keep that
        portion of the spectral loop fast, however, we figure out here which ones we need.
        Doing this inside the spectral loop leads to significant performance loss!

    =#
    spec_dep_wfuncs = [sve for sve in keys(rt.wfunctions_map)
        if sve isa AbstractStateVectorElement]


    if length(spec_dep_wfuncs) > 0
        #=
            If we have to calculate aerosol-related Jacobians, we also need the linearized
            inputs for XRTM for the phase function expansion coefficients. To minimize
            allocations, we create empty arrays here and re-use them in the `set_XRTM_wf`
            functions inside the hires spectral loop.

            This will be a list of list of matrices.
            l_coef[t][l] is a 2D matrix for thread `t` for layer `l`.

            Remember, we have to have differently-sized l_coef for every layer due to the
            fact that each layer might have a different number of effective phase function
            coefficients.

        =#
        l_coef_threads = Vector{Matrix{Float64}}[]

        for t in 1:Threads.nthreads()
            l_coef_layers = Matrix{Float64}[]

            for l in 1:rt.scene.atmosphere.N_layer

                n_coef = XRTM.get_n_coef(xrtm_l[1], l-1)
                l_coef = zeros(n_coef, n_elem)

                push!(l_coef_layers, l_coef)

            end

            push!(l_coef_threads, l_coef_layers)
        end

    end


    #=
        We can direct the spectral loop to only evaluate certain points, which is very
        useful if we want to speed up the process using some interpolation method, such as
        the ACOS-type "non-uniform sampling" strategy (NUS).

        Note that the spectral points (for radiances and weighting functions or Jacobians)
        must be recovered elsewhere, as it is out-of-scope for this function.

        Users must present a a boolean array (or integer array) in which the points to be
        evaluated are set to be `true` (or 1), and skipped points are `false` (or 0).
        The array must be stored in the options dictionary under the key `sampling`.

    =#

    if haskey(options_dict, "sampling") & !(swin isa BinnedSpectralWindow)
        # Make sure we don't skip points for a BinnedSpectralWindow type.
        spec_iterator = findall(options_dict["sampling"])
    else
        spec_iterator = 1:N_hires
    end

    # Is this the first XRTM call?
    # Keep in mind this must be a vector since every XRTM instance must do this once in a
    # threaded environment.
    # (some XRTM calls or computations only need to be done once!)
    first_XRTM_call = ones(Bool, Threads.nthreads())

    # Do we calculate weighting functions?
    have_jacobians = "calc_derivs" in options_dict["options"]

    # We need some temporary arrays to do various calculations, unfortunately
    tmp_vec_list = [zeros(N_layer) for x in 1:Threads.nthreads()]

    desc_str = "(Nthread=$(Threads.nthreads())) XRTM loop $(swin) for solver(s): " *
        "$(join(options_dict["solvers"], ", "))"

    # Decide whether we want progress bars for XRTM calculations
    XRTM_PROGRESS = false
    if haskey(ENV, "XRTM_PROGRESS") && (ENV["XRTM_PROGRESS"] == "1")
        XRTM_PROGRESS = true
    end

    prog = Progress(length(spec_iterator);
        dt=0.25, barlen=10, desc=desc_str,
        enabled=XRTM_PROGRESS
    )

    thread_error_flags = zeros(Bool, Threads.nthreads())
    looplock = Threads.SpinLock()

    # Hires spectral loop
    Threads.@threads for i_spectral in spec_iterator

        #=
            Notes on threading

            Inside this loop, there are only a few writing operations. We first start
            with a list of XRTM instances - each to be used by one thread only. The same
            holds for "tmp_vec", which we generate ahead of this loop and then simply grab
            the one we want (corresponding to the thread ID).

            A threaded loop via Threads.@threads cannot be stopped or broken out of. So
            if one thread encounters some error due to bad inputs into XRTM (e.g., the
            surface kernel amplitude is negative), we cannot simply `break` out of it.
            In that case, we flag this thread as having 'errored', and let it continue
            without further evaluating the loop body. Then, after the loop is processed,
            we check for the flags for all threads and proceed as appropriate.

        =#


        if thread_error_flags[Threads.threadid()]
            # An error has been noticed for this thread -> skip evaluating the loop body.
            continue
        end


        # Pick the XRTM instance for this thread!
        if Threads.nthreads() == length(xrtm_l)
            xrtm = xrtm_l[Threads.threadid()]
        else
            xrtm = xrtm_l[1]
        end

        # Pick the tmp vector for this thread!
        tmp_vec = tmp_vec_list[Threads.threadid()]

        # Sun-normalized
        XRTM.set_F_0(xrtm, 1.0)

        #=
            Move total optical depths into XRTM
            (copy contents into tmp_vec and then call function)
        =#
        @views tmp_vec[:] = rt.optical_properties.total_tau[i_spectral,:]
        # Force τ to be [0, ∞)
        @turbo for l in eachindex(tmp_vec)
            tmp_vec[l] = max(1e-10, tmp_vec[l])
        end

        XRTM.set_ltau_n(xrtm, tmp_vec)

        #=
            Move total single-scatter albedo into XRTM
            (copy contents into tmp_vec and then call function)
        =#

        @views tmp_vec[:] = rt.optical_properties.total_omega[i_spectral,:]
        # Force ω to be (0,1]
        @turbo for l in eachindex(tmp_vec)
            tmp_vec[l] = max(0.0, min(0.999999, tmp_vec[l]))
        end

        XRTM.set_omega_n(xrtm, tmp_vec)

        #=
            Move surface kernel factor into XRTM. Here we have to evaluate the
            value for the given spectral index.
        =#

        # Loops through all kernels and evalutes them at their respective spectral
        # point `i_spectral`. Note that the order of BRDF kernels is maintained
        # such that the k'th BRDF kernel belongs to the XRTM kernel `k-1`
        # (due to XRTM's 0-based indexing.)
        for (k, kernel) in enumerate(surface.kernels)

            ampfac = evaluate_surface_at_idx(kernel, i_spectral)

            if (ampfac <= 0) | (ampfac >= 1.0)
                @error "[XRTM] Surface kernel amplitude factor not in (0,1)."
                thread_error_flags[Threads.threadid()] = true
                break
            end

            XRTM.set_kernel_ampfac(xrtm, k-1, ampfac)
        end

        if thread_error_flags[Threads.threadid()]
            continue
        end

        # Let XRTM know which layers involve weighting functions
        if first_XRTM_call[Threads.threadid()] & have_jacobians
            # This function needs to be called only once (as per XRTM manual)
            XRTM.update_varied_layers(xrtm)
        end

        #=
            Calculate radiances!
            Upwelling radiance (I_up),
            downwelling radiance (I_dn),
            upwelling weighting functions (K_up),
            downwelling weighting functions (K_dn)
        =#

        for solver in options_dict["solvers"]
            # Execute XRTM to get radiances and (optionally) weighting functions!

            # Results are:
            # upwelling intensity, downwelling intensity,
            # upwelling weighting functions, downwelling weighting functions

            I_up, I_dn, K_up, K_dn = XRTM.radiance(
                xrtm,
                solver,
                size(out_phis, 1),
                out_phis
            )


            # DEBUG
            #=
            I_up, I_dn, K_up, K_dn = [
                rand(Float64, (Int(n_stokes),1,1,1)),
                rand(Float64, (Int(n_stokes),1,1,1)),
                rand(Float64, (Int(n_stokes),1,1,Int(n_derivs),1)),
                rand(Float64, (Int(n_stokes),1,1,Int(n_derivs),1)),
            ]
            =#

            # Flip the sign of Stokes-U due to convention difference

            if (n_stokes >= 3) && flip_U
                for s in [3]
                    I_up[s,:,:,:] .*= -1
                    I_dn[s,:,:,:] .*= -1
                    K_up[s,:,:,:,:] .*= -1
                    K_dn[s,:,:,:,:] .*= -1
                end
            end

            #=
                Store radiances depending on observation type, add them to the total
                contributions.
                Note: this is a function call that seems unneccesarily expressive w.r.t.
                the number of arguments. Experiments have shown, however, that this is a
                very performance critical part of the monochromatic loop, and Julia is
                seemingly better able to optimize this function when written like it is.

                Note that we have to lock this, so that in case of running with threads,
                we don't trip on the toes of another thread..
            =#

            lock(looplock)
            _add_xrtm_results!(
                rt.hires_radiance, # Radiance object
                rt.hires_wfunctions, # Vector of Radiance objects
                have_jacobians, # Bool (if we calculate weighting functions .. true)
                rt.scene.observer, # The observer type so we know which results to use
                i_spectral, # spectral index
                n_stokes,
                n_derivs,
                I_up, I_dn, K_up, K_dn # The XRTM results
            )
            unlock(looplock)

        end # Solver loop end

        # Establish that we have succesfully called XRTM once!
        first_XRTM_call[Threads.threadid()] = false

        # Update progress bar
        next!(prog)

    end # Spectral loop end

    # For now, set all radiances and jacobians to NaNs if errors were encountered
    if any(thread_error_flags)
        @warn "Errors were encountered in this XRTM run. Setting all to NaN!"
        rt.hires_radiance[:] .= NaN
        for jac in rt.hires_wfunctions
            jac[:] .= NaN
        end
    end


end

"""
$(TYPEDSIGNATURES)

Specialized function to copy XRTM results into respective containers for a
`SatelliteObserver` observer.
"""
function _add_xrtm_results!(
    radiance::Radiance,
    wfunctions,
    have_jacobians::Bool,
    observer::SatelliteObserver,
    i_spectral::Integer,
    n_stokes::Integer,
    n_derivs::Integer,
    I_up::Array,
    I_dn::Array,
    K_up::Array,
    K_dn::Array
)

    # Radiances are always copied
    for s in 1:n_stokes # Loop over Stokes elements
        radiance.S[i_spectral,s] += I_up[s,1,1,1]
    end

    !have_jacobians && return

    # Copy weighting function results
    for i in 1:n_derivs
        for s in 1:n_stokes
            wfunctions[i].S[i_spectral,s] += K_up[s,1,1,i,1]
        end
    end

end

"""
$(TYPEDSIGNATURES)

Specialized function to copy XRTM results into respective containers for a
`UplookingGroundObserver` observer.
"""
function _add_xrtm_results!(
    radiance::Radiance,
    wfunctions,
    have_jacobians::Bool,
    observer::UplookingGroundObserver,
    i_spectral::Integer,
    n_stokes::Integer,
    n_derivs::Integer,
    I_up,
    I_dn,
    K_up,
    K_dn
)

    # Radiances are always copied
    for s in 1:n_stokes # Loop over Stokes elements
        radiance.S[i_spectral,s] += I_dn[s,1,1,1]
    end

    !have_jacobians && return

    # Copy weighting function results
    for i in 1:n_derivs
        for s in 1:n_stokes
            wfunctions[i].S[i_spectral,s] += K_dn[s,1,1,i,1]
        end
    end


end


"""
$(TYPEDSIGNATURES)

Calculates the numnber of phasefunction expansion moments needed for all layers in
this `MonochromaticRTMethod` object `rt`. The threshold specifies the fraction of
optical depth below which the contributions of that scatterer will be ignored.

# Details

This function simply loops through all layers, and calculates the total scattering optical
depth (all aerosols plus Rayleigh). Then, for each scatterer (all aerosols, Rayleigh),
it is determined what fraction of the total scattering optical depth it accounts for. If
the fraction is larger than some threshold (given by the keyword parameter `threshold`),
then the number of needed coefficients for this layer is updated to include those of that
scatterer.
"""
function _calculate_needed_n_coef(rt::MonochromaticRTMethod; threshold=0.05)

    # XRTM wants Int32 arrays for n_coef
    n_coef_array = ones(Int32, rt.scene.atmosphere.N_layer)

    # Short-cuts
    swin = rt.optical_properties.spectral_window
    opt = rt.optical_properties

    # We pick the appropriate spectral index of the window to calculate scattering
    # optical depth.
    i_spectral = get_scattering_index(swin)

    for l in 1:rt.scene.atmosphere.N_layer

        # Calculate the total optical depth due to scattering
        # (for this layer `l`)
        od_scatt = opt.total_tau[i_spectral, l] * opt.total_omega[i_spectral, l]

        # Now see how much of that is due to Rayleigh
        # (which is, for now, always 3 coefficients)
        if (opt.rayleigh_tau[i_spectral, l] / od_scatt) > threshold
            n_coef_array[l] = 3
        end

        # Now check for each aerosol type
        for (aer, aer_tau) in opt.aerosol_tau
            if (aer_tau[i_spectral, l] / od_scatt) > threshold

                # If this scatterer has more coefs than what is already provided,
                # update the value!
                n_coef_array[l] = max(
                    n_coef_array[l],
                    aer.optical_property.max_coefs
                )

            end
        end

    end

    return n_coef_array
end