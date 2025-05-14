#=

    Low-Streams Interpolation (LSI) is a technique by Chris O'Dell to accelerate RT
    calculations for hyperspectral applications. This code is a direct implementation
    of the publication 10.1029/2009JD012803.

    Note that equation A3 has a typo: the x-bin indices in the numerator are switched,
    it should be

    ε(τ_i, x_A) = [ε(τ_i, x_i, 2) - ε(τ_i, x_i, 1)] / [x_(i,2) - x_(i,1)] etc.

=#


"""
Pretty printing for Low-Streams Interpolation RT

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", rt::LSIRTMethod)
    # Just print this for now
    println(io, "Low-Streams Interpolation RT $(rt.monochromatic_RT.model)")

end

"""
Pretty printing for Low-Streams Interpolation RT

$(SIGNATURES)
"""
function show(io::IO, rt::LSIRTMethod)
    # Just print this for now
    println(io, "Low-Streams Interpolation RT $(rt.monochromatic_RT.model)")

end

"""
    Calculates √ξ (see O'Dell 2009) for a spectral index `i` from a given
    `AbstractOpticalProperties` object. Needed for LSI.

    $(TYPEDSIGNATURES)

    NOTE: This current implementation is not very performant, so a specialized
    implementation to calculate √ξ for all spectral points in a band is needed
    to make this faster. However, it seems to run in less than 10ms and allocate
    about ~1MB of memory, which is currently still acceptable.
    The major drawback is the fact that if one runs this function for all spectral
    points in a band, the internal loops here reduce the performance and scaling. It would
    be better to run the spectral loop inside all these layer loops.

"""
function calculate_ξ_sqrt(
    opt::AbstractOpticalProperties,
    i::Integer
    )

    N_layer = size(opt.total_tau, 2)

    # Total gas absorption optical depth, per layer
    τ_gas_layer = opt.tmp_Nlay1
    @views τ_gas_layer[:] .= 0.0

    τ_gas = 0.0
    for gas_tau in values(opt.gas_tau)
        for l in 1:N_layer
            τ_gas_layer[l] += gas_tau[i,l]
            τ_gas += gas_tau[i,l]
        end
    end

    # Total column scattering optical depth
    τ_sca = 0.0
    # Cumulative sum starting from TOA
    τ_sca_running = opt.tmp_Nlay2
    @views τ_sca_running[:] .= 0

    for l in 1:N_layer
        τ_sca_running[l] += opt.total_tau[i,l] * opt.total_omega[i,l]
    end

    # Calculate cumulative sum, so from here on τ_sca_running[l] is the sum
    # of all preceding scattering optical depths, counted from the top of the
    # atmosphere (l=1).
    for l in 2:N_layer
        τ_sca_running[l] += τ_sca_running[l-1]
    end

    # Set critical value, down to which we integrate the gas optical depth
    # to get τ_g'
    if τ_sca <= 2.0
        critical_value = τ_sca / 2.0
    else
        critical_value = 1.0
    end

    τ_gas_prime = 0.0
    # Let this loop not go all the way down to the bottom layer.
    for l in 1:N_layer-1

        τ_gas_prime += τ_gas_layer[l]

        if τ_sca_running[l] >= critical_value
            break
        end


    end

    return sqrt(τ_gas_prime / τ_gas)

end

"""
Generates the `LSIRTMethod` object from an existing `MonochromaticRTMethod`
and some additional parameters.
"""
function LSIRTMethod(
    bin_boundaries::Vector,
    rt::MonochromaticRTMethod,
    high_options::Union{T, Vector{T}}
    #τ_crit_min=0.0,
    #τ_crit_max=Inf
) where {T <: AbstractDict}

    #=
        Calculate the assignment of spectral points to bins. Points outside
        the bin ranges are assigned the lowest or highest bin.
    =#

    # First, calculate the total gas optical depth, and use a tmp array for
    # the calculation
    τ_gas = rt.optical_properties.tmp_Nhi1
    N_hires = rt.optical_properties.spectral_window.N_hires

    # We need to allocate a new array for this since this object will
    # passed to the new instance of LSIRTMethod.
    ξ_sqrt = zeros(N_hires)

    # Sum up all gas contributions for τ_gas
    @views τ_gas[:] .= 0
    for (gas, gas_tau) in rt.optical_properties.gas_tau
        avx_add_along_columns!(τ_gas, gas_tau)
    end

    # Shift τ_gas == 0 to some tiny value, the algorithm might have
    # problems with values that are exactly zero.
    τ_gas[τ_gas .<= 1e-10] .= 1e-10


    #=
        Computation of ξ = √(τ_g' / τ_g).
    =#

    ξ_sqrt[:] .= calculate_ξ_sqrt.(Ref(rt.optical_properties), 1:N_hires)


    # In case the bin boundaries are not sorted, let's sort them here
    sort!(bin_boundaries)
    tau_bin_assignment = searchsortedfirst.(Ref(bin_boundaries), τ_gas) .- 1

    ξ_sqrt_bin_assignment = similar(tau_bin_assignment)
    ξ_sqrt_bin_assignment[:] .= 0

    N_tau_bins = length(bin_boundaries) - 1
    N_ξ_bins = 2

    used_bin = zeros(Bool, N_tau_bins, N_ξ_bins)

    tau_bins = unique(tau_bin_assignment)
    sort!(tau_bins)

    # For every gas-tau bin, we now make a decision about the √ξ bin
    for i_bin in 1:N_tau_bins

        # Collect all points for this τ_gas bin
        idx = tau_bin_assignment .== i_bin

        # If this bin falls outside of the values set by the
        # critical τ parameter, then we don't split up this
        # τ-bin into two √ξ bins.

        #=
        # This might be used later if we implement single-ξ bins
        bin_left = max(1, i_bin)
        bin_right = min(length(bin_boundaries), i_bin + 1)

        if (τ_crit_min > bin_boundaries[bin_left]) |
            (τ_crit_max < bin_boundaries[bin_right])

            # Pul all into a single bin
            ξ_sqrt_bin_assignment[idx] == 1
            continue

        end
        =#

        if sum(idx) < 1
            continue
        end


        # Get the 25% and 75% percentiles of √ξ of the subsets
        ξ_25, ξ_75 = quantile.(Ref(ξ_sqrt[idx]), [0.25, 0.75])

        # "Lower" bin: √ξ < 25% percentile
        idx_25 = idx .& (ξ_sqrt .< ξ_25)
        # "Upper" bin: √ξ < 75% percentile (but not in lower bin)
        idx_75 = idx .& (ξ_sqrt .<= ξ_75) .& (ξ_sqrt .>= ξ_25)

        ξ_sqrt_bin_assignment[idx_25] .= 1
        ξ_sqrt_bin_assignment[idx_75] .= 2

        used_bin[i_bin, :] .= true

    end

    #=
        Create new RT Method objects to perform the binned calculations
    =#
    idx_center = get_scattering_index(rt.optical_properties.spectral_window)
    idx_edge = 1

    RT_bin = create_rt_copy_for_bins(rt, idx_center)
    RT_bin_edge = create_rt_copy_for_bins(rt, idx_edge)

    #=
        Create result arrays
    =#




    bin_rad_lo = Array{typeof(rt.hires_radiance)}(undef, N_tau_bins, N_ξ_bins)
    bin_wf_lo = Array{Vector{typeof(rt.hires_radiance)}}(undef, N_tau_bins, N_ξ_bins)
    bin_rad_hi = Array{typeof(rt.hires_radiance)}(undef, N_tau_bins, N_ξ_bins)
    bin_wf_hi = Array{Vector{typeof(rt.hires_radiance)}}(undef, N_tau_bins, N_ξ_bins)
    bin_edge_rad_lo = Array{typeof(rt.hires_radiance)}(undef, 1, 1)
    bin_edge_rad_hi = Array{typeof(rt.hires_radiance)}(undef, 1, 1)

    # Element type of radiance array (e.g. Float64)
    Trad = eltype(rt.hires_radiance)
    # Short-cut to the type of the radiance, allows us to construct new ones
    RadType = typeof(rt.hires_radiance).name.wrapper

    for ar in [bin_rad_lo, bin_rad_hi, bin_edge_rad_lo, bin_edge_rad_hi]
        for i in eachindex(ar)
            _v = RadType(Trad, 1)
            ar[i] = _v
        end
    end




    #=
        Create the LSI RT method object
    =#

    return LSIRTMethod(
        high_options,
        bin_boundaries,
        τ_gas,
        ξ_sqrt,
        tau_bin_assignment,
        ξ_sqrt_bin_assignment,
        used_bin,
        rt,
        RT_bin,
        RT_bin_edge,
        bin_rad_lo,
        bin_wf_lo,
        bin_rad_hi,
        bin_wf_hi,
        bin_edge_rad_lo,
        bin_edge_rad_hi
    )

end

"""
$(TYPEDSIGNATURES)

Calculates the binned optical properties needed to perform the LSI correction.
"""
function calculate_binned_properties!(
    lsi::LSIRTMethod,
    tau_gas_bin::Integer,
    ξ_bin::Integer,
    do_edge::Bool
    )

    @debug "Calculating binned optical properties for $(tau_gas_bin) / $(ξ_bin)"

    # Create some shortcuts
    old_rt = lsi.monochromatic_RT

    if do_edge
        rt = lsi.RT_bin_edge
    else
        rt = lsi.RT_bin
    end

    this_swin = rt.optical_properties.spectral_window

    @assert this_swin isa BinnedSpectralWindow "Must be a BinnedSpectralWindow!"

    # Note this only works if this_swin is a BinnedSpectralWindow,
    # which should always be the case. The spectral_idx field of a
    # BinnedSpectralWindow always refers to the e.g. band-center of the original
    # spectral window it is derived from.

    spectral_idx = this_swin.spectral_idx

    # Select all points that fall into this τ/√ξ-bin
    idx = (
        (lsi.tau_gas_bin_assignment .== tau_gas_bin) .&
        (lsi.ξ_sqrt_bin_assignment .== ξ_bin)
    )

    #=
        For LSI, the binned optical properties are constructed the following way.
        Optical depth due to gas absorption is the mean of all gas optical depths
        of the spectral points that belong to a certain bin.
    =#

    # Set to zero, just in case
    for (gas, gas_tau) in rt.optical_properties.gas_tau
        @views gas_tau[:,:] .= 0
    end

    for (gas, gas_tau) in old_rt.optical_properties.gas_tau
        @views rt.optical_properties.gas_tau[gas][1,:] .+=
            mean(gas_tau[idx,:], dims=1)[1,:]
    end

    #=
        Optical depth due to scattering (Rayleigh, aerosols) is constructed as a
        profile taken straight from the band center
    =#

    rt.optical_properties.rayleigh_tau[1,:] =
        old_rt.optical_properties.rayleigh_tau[spectral_idx,:]

    for (aer, aer_tau) in old_rt.optical_properties.aerosol_tau
        rt.optical_properties.aerosol_tau[aer][1,:] = aer_tau[spectral_idx,:]
        # We must also grab the aerosol single-scatter albedo
        rt.optical_properties.aerosol_omega[aer][1,:] =
            old_rt.optical_properties.aerosol_omega[aer][spectral_idx,:]
    end

    #=
        For now, we don't have to do anything with coefficients, since those are
        evaluated at the band center anyway..

        Finally, compute all "total" optical properties from the individual
        contributions.
    =#

    opt = rt.optical_properties
    @views opt.total_tau[:,:] .= 0.0

    # Add up total optical depth from contributions of each gas
    for (gas, gas_tau) in opt.gas_tau
        @views opt.total_tau[:,:] += gas_tau[:,:]
    end

    # Add up aerosol contributions
    for (aer, aer_tau) in opt.aerosol_tau
        @views opt.total_tau[:,:] += aer_tau[:,:]
    end

    # Add Rayleigh contributions
    # (this will be zero if no RayleighScattering is present..)
    @views opt.total_tau[:,:] += opt.rayleigh_tau[:,:]

    #=
        Calculate total ω
    =#

    @views opt.total_omega[:,:] .= 0.0

    # Start with Rayleigh
    @views opt.total_omega[:,:] += opt.rayleigh_tau[:,:]
    # Add each aerosol scattering optical depths
    for (aer, aer_tau) in opt.aerosol_tau
        @views opt.total_omega[:,:] += opt.aerosol_omega[aer][:,:] .* aer_tau[:,:]
    end
    # Divide by total extinction tau
    @views opt.total_omega[:,:] ./= opt.total_tau[:,:]

end


function perform_LSI_correction!(lsi::LSIRTMethod)


    #=
        TODO: This is currently hard-coded for XRTM, since there are some performance
        gains to be made if we know the RT is performed by XRTM.
    =#

    # Decide whether we want progress bars for XRTM calculations
    XRTM_PROGRESS = false
    if haskey(ENV, "XRTM_PROGRESS") && (ENV["XRTM_PROGRESS"] == "1")
        XRTM_PROGRESS = true
    end

    #time_bin = @timed begin

        #=
            We have to loop over the number of XRTM configurations known to us; they are
            incompatible, so we must create XRTM instances for each.

            Note that we are only doing this explicitly so we do not have to create two
            new XRTM instances for every bin calculation, which is easily possible, but
            very slow.

            Roughly, the steps are these:

            1) all lsi.low_RT and all lsi.high_RT objects have the same sets of
            configurations. We grab an arbitrary one, and extract the `model_options`
            field. The model options can be either a Vector{Dict} or a Dict, so we must
            first check and create a new one-element Vector{Dict}, if necessary, so we
            can loop over.

            2) The outer-most loop iterates over the known XRTM configurations. For each
            iteration, a new XRTM instance is created.

            3) Inside the main options loop, we call the function to calculate radiances
            and weighting functions, and supply the recently created XRTM instance, so
            that inside of `_calculate_radiances_and_wfs_XRTM!`, the XRTM instance does
            not have to be created again. *THIS IS THE MAIN SPEED-UP*

            4) At the end of the bin calculations, the XRTM instance is destroyed to free
            up memory (explicitly needed since Julia cannot garbage collect it!). The same
            process is done separately for the low-accuracy and high-accuracy
            calculations!


            NOTE
            These calculations cannot easily be threaded, since each binned calculation
            relies on its own `MonochromaticRTMethod` object with a single spectral point
            each, and thus spreading out the calculations on multiple threads requires
            copying those objects, which currently is not implemented (and probably never
            will be).

        =#


        # Make all model options into a list so we can iterate over them
        low_options = Dict[]
        high_options = Dict[]
        if lsi.monochromatic_RT.model_options isa Dict
            push!(low_options, lsi.monochromatic_RT.model_options)
        else
            push!(low_options, lsi.monochromatic_RT.model_options...)
        end

        if lsi.high_options isa Dict
            push!(high_options, lsi.high_options)
        else
            push!(high_options, lsi.high_options...)
        end

        # Low bins
        for this_option in low_options

            xrtm_low = create_XRTM(
                lsi.RT_bin,
                lsi.RT_bin.scene.observer,
                this_option
                )

            # Loop through all tau and xi bins and calcuate radiances!
            for τ_bin in axes(lsi.used_bin, 1), ξ_bin in axes(lsi.used_bin, 2)

                # Skip unused bin
                if !lsi.used_bin[τ_bin, ξ_bin]
                    continue
                end

                calculate_binned_properties!(lsi, τ_bin, ξ_bin, false)


                clear!(lsi.RT_bin)  # Must clear out beforehand
                _calculate_radiances_and_wfs_XRTM!(
                    lsi.RT_bin,
                    lsi.RT_bin.scene.observer,
                    this_option,
                    xrtm_in=xrtm_low
                    )

                # Results are ADDED to the total from the various solvers
                @views lsi.bin_rad_lo[τ_bin, ξ_bin][:] += lsi.RT_bin.hires_radiance[:]

            end

            calculate_binned_properties!(lsi, 1, 1, true)
            _calculate_radiances_and_wfs_XRTM!(
                lsi.RT_bin_edge,
                lsi.RT_bin_edge.scene.observer,
                this_option,
                xrtm_in=xrtm_low
                )


            @views lsi.bin_edge_rad_lo[1, 1][:] += lsi.RT_bin.hires_radiance[:]

            XRTM.destroy(xrtm_low)
        end

        # High bins
        for this_option in high_options

            xrtm_high = create_XRTM(
                lsi.RT_bin,
                lsi.RT_bin.scene.observer,
                this_option
                )

            # Loop through all tau and xi bins and calcuate radiances!
            for τ_bin in axes(lsi.used_bin, 1), ξ_bin in axes(lsi.used_bin, 2)

                # Skip unused bin
                if !lsi.used_bin[τ_bin, ξ_bin]
                    continue
                end

                calculate_binned_properties!(lsi, τ_bin, ξ_bin, false)

                clear!(lsi.RT_bin) # Must clear out beforehand
                _calculate_radiances_and_wfs_XRTM!(
                    lsi.RT_bin,
                    lsi.RT_bin.scene.observer,
                    this_option,
                    xrtm_in=xrtm_high
                    )

                # Results are ADDED to the total from the various solvers
                @views lsi.bin_rad_hi[τ_bin, ξ_bin][:] += lsi.RT_bin.hires_radiance[:]


            end


            calculate_binned_properties!(lsi, 1, 1, true)
            _calculate_radiances_and_wfs_XRTM!(
                lsi.RT_bin_edge,
                lsi.RT_bin_edge.scene.observer,
                this_option,
                xrtm_in=xrtm_high
                )

            @views lsi.bin_edge_rad_hi[1, 1][:] += lsi.RT_bin.hires_radiance[:]

            XRTM.destroy(xrtm_high)
        end


    #=


        _low_rt = first(values(lsi.low_RT)) # Grab some RT object, doesn't matter which
        if _low_rt.model_options isa Dict
            low_options = [_low_rt.model_options]
        else
            low_options = _low_rt.model_options
        end
        _high_rt = first(values(lsi.high_RT)) # Grab some RT object, doesn't matter which
        if _high_rt.model_options isa Dict
            high_options = [_high_rt.model_options]
        else
            high_options = _high_rt.model_options
        end
        # Low bins
        for this_option in low_options

            xrtm_low = create_XRTM(
                _low_rt,
                _low_rt.scene.observer,
                this_option
                )

            desc_str="LSI low bin calculations ($(this_option["solvers"]))"
            prog = Progress(length(keys(lsi.low_RT));
                dt=0.25, barlen=10, desc=desc_str,
                enabled=XRTM_PROGRESS
                )

            for ((tau_bin, ξ_bin), rt) in lsi.low_RT
                # ingore √ξ == 0 for RT
                ξ_bin == 0 && continue
                _calculate_radiances_and_wfs_XRTM!(
                    rt,
                    rt.scene.observer,
                    this_option,
                    xrtm_in=xrtm_low
                    )

                # Let the user know if there are any NaNs etc.
                check_non_finites(rt)
                # Update progress bar
                next!(prog)
            end

            # Low edge bin
            for ((tau_bin, ξ_bin), rt) in lsi.low_RT_edge
                _calculate_radiances_and_wfs_XRTM!(
                    rt,
                    rt.scene.observer,
                    this_option,
                    xrtm_in=xrtm_low
                    )
                # Let the user know if there are any NaNs etc.
                check_non_finites(rt)
            end

            # We don't need this XRTM instance anymore..
            XRTM.destroy(xrtm_low)
        end

        # High bins
        for this_option in high_options

            xrtm_high = create_XRTM(
                _high_rt,
                _high_rt.scene.observer,
                this_option
            )

            desc_str="LSI high bin calculations ($(this_option["solvers"]))"
            prog = Progress(length(keys(lsi.high_RT));
                dt=0.25, barlen=10, desc=desc_str,
                enabled=XRTM_PROGRESS
            )

            for k in collect(eachindex(lsi.high_RT))

                # We have to do it this way because threads does not like Dict keys..
                tau_bin, ξ_bin = k
                rt = lsi.high_RT[k]

                # ingore √ξ == 0 for RT
                ξ_bin == 0 && continue

                _calculate_radiances_and_wfs_XRTM!(
                    rt,
                    rt.scene.observer,
                    this_option,
                    xrtm_in=xrtm_high
                )

                # Let the user know if there are any NaNs etc.
                check_non_finites(rt)
                # Update progress bar
                next!(prog)
            end

            # High edge bin (there should only be one here)
            for ((tau_bin, ξ_bin), rt) in lsi.high_RT_edge
                _calculate_radiances_and_wfs_XRTM!(
                    rt,
                    rt.scene.observer,
                    this_option,
                    xrtm_in=xrtm_high
                )

                # Let the user know if there are any NaNs etc.
                check_non_finites(rt)
            end

            # We don't need this XRTM instance anymore..
            XRTM.destroy(xrtm_high)
        end

    #end # End timed section

    @info "LSI binned calculations done in $(time_bin.time) sec."

    #=
        Construct the error array for each bin, as well as calculate the total gas
        optical depth along with the actual √ξ for the binned optical properties.

        Important NOTE!
        ###############
        Most "containers" here are dictionaries rather than arrays, and the dictionary
        keys tend to be either the bin in tau-gas-space [i], in √ξ-space [j], or both
        [(i,j)]. This also means that, in the case of emtpy bins, the corresponding key
        does not exist in these containers. For example, if bin 4, which lies between
        bin boundaries [4,5] has no spectral points in them, there will be no ε[(4,j)].
        This means, that when traversing bins, one must always iterate over existing bins
        via, e.g. keys(lsi.low_RT) or unique(lsi.tau_gas_bin_assingment).

        Funny thing, also, in Julia you can access a dictionary element d[(i,j)] also
        via omitting the tuple parentheses, i.e.: d[i,j] and d[(i,j)] are the same.

    =#

    # Element type of radiance array (e.g. Float64)
    T = eltype(lsi.monochromatic_RT.hires_radiance)
    # Short-cut to the type of the radiance, allows us to construct new ones
    RadType = typeof(lsi.monochromatic_RT.hires_radiance).name.wrapper

    # Error of Stokes components
    ε = Dict()
    # Partial derivative of ε w.r.t. weighting function ∂ε/∂x
    ∂ε = Dict()
    # Bin-averaged ξ
    ξ = Dict()
    # Bin-averaged τ_gas
    τ_gas = Dict()
    # Which bins should use relative errors? (this is for error catching)
    use_scaled_error = Dict()
    # How many weighting functions do we have?
    k_first = first(keys(lsi.low_RT))
    N_wf = map(length, values(lsi.low_RT[k_first].wfunctions_map)) |> sum
    N_stokes = size(lsi.monochromatic_RT.hires_radiance.S, 2)

    bad_bins = Int[]

    #=
        First, create some handy short cuts that contain the bin errors for both
        Stokes components ε[(i,j)], as well as the weighting functions ε_wf[(i,j)].
    =#
    for k in keys(lsi.low_RT)
        # k = (i_bin, ξ_bin)

        # Subtract all Stokes components
        ε[k] = RadType(T, 1)
        @views ε[k].S[1,:] = (
            lsi.low_RT[k].hires_radiance.S[1,:] -
            lsi.high_RT[k].hires_radiance.S[1,:]
        )

        # Produce relative errors by dividing by the intensity, such that
        # ε[(i,j)] = (low_I[(i,j)] - high_I[(i,j)]) / high_I[(i,j)]
        # and
        # ε[(i,j)] = (low_S[(i,j)] - high_S[(i,j)]) for S = Q,U
        #=
            NOTE !!!

            Unlike in the publication, we only produce relative errors for the intensity
            component, but not for Q, U. This is also the way it is implemented in ACOS.
            Possibly due to the fact that the (1-ε) term in the denominator could grow
            to be very large if ε ≈ 1.

        =#

        fac = abs(ε[k].S[1,1]) / abs(lsi.high_RT[k].hires_radiance.S[1,1])

        if fac > 1

            #@info "Bin error larger than high-bin radiance for $(k)! Ratio: $(fac)"
            #@info "Bin tau is between $(lsi.bin_boundaries[k[1]]) and $(lsi.bin_boundaries[k[1]+1])"
            ε[k].S[1,1] /= lsi.high_RT[k].hires_radiance.S[1,1]
            push!(bad_bins, k[1])
        else
            # Otherwise: OK, we can move on and scale the error to produce
            # ε[(i,j)] = (low_I[(i,j)] - high_I[(i,j)]) / high_I[(i,j)]

            ε[k].S[1,1] /= lsi.high_RT[k].hires_radiance.S[1,1]

        end



        # Calculate bin errors for derivatives
        ∂ε[k] = [RadType(T, 1) for l in 1:N_wf]
        r_high = lsi.high_RT[k].hires_radiance # this is the high-stream radiance result

        for l in 1:N_wf

            wf_low = lsi.low_RT[k].hires_wfunctions[l]
            wf_high = lsi.high_RT[k].hires_wfunctions[l]

            # Again, we produce the scaled error for intensity only, not for Stokes

            ∂ε[k][l].S[1,1] = 1/r_high.S[1,1] * (
                wf_low.S[1,1] - wf_high.S[1,1] * (1 + ε[k].S[1,1])
            )

            if RadType == VectorRadiance
                ∂ε[k][l].S[1,2:end] = wf_low.S[1,2:end] - wf_high.S[1,2:end]
            end

        end

        # Also, calculate the bin-averaged τ_gas and √ξ
        opt = lsi.low_RT[k].optical_properties
        # Remember, binned optical properties only have one spectral point and the optical
        # properties of both low_RT and high_RT ar the same (or should be).
        ξ[k] = calculate_ξ_sqrt(opt, 1)

        # Total column tau gas in this bin, extracted from the optical properties
        # object that is attached to the low/high RT object
        τ_gas[k] = sum(sum(gas_tau) for gas_tau in values(opt.gas_tau))

    end

    # Create the interpolator object, but note that we are interpolating in log(τ) rather
    # than just τ.

    nodes_τ = Float64[]
    nodes_ξ = Float64[]
    nodes_err = Float64[]

    for k in keys(lsi.low_RT)

        push!(nodes_τ, (τ_gas[k]))
        push!(nodes_ξ, ξ[k])
        push!(nodes_err, ε[k].I[1])

    end

    # Equation A1 in O'Dell
    # (center-corrected errors, only a function of τ-bin, not of ξ)
    ε_τ = Dict()
    ∂ε_τ = Dict()
    mean_τ = Dict()

    used_τ_bins = unique(lsi.tau_gas_bin_assignment)
    first_τ_bin = minimum(used_τ_bins)
    sort!(used_τ_bins)

    for i in used_τ_bins

        # Number of spectral points in each sub-bin (ξ_bin == 1 and ξ_bin == 2) that
        # belong to a certain tau bin.
        idx_τ = lsi.tau_gas_bin_assignment .== i
        N1 = sum((idx_τ) .& (lsi.ξ_sqrt_bin_assignment .== 1))
        N2 = sum((idx_τ) .& (lsi.ξ_sqrt_bin_assignment .== 2))

        # This (A1) is needed for the "centering correction", and not well explained in
        # the manuscript ..

        r_low = RadType(T, 1)
        r_high = RadType(T, 1)

        # Some bins might not have valid entries in the ξ_bin == 2 dimension
        @views r_low.S[:,:] = (
            lsi.low_RT[i,1].hires_radiance.S[:,:] * N1 +
            lsi.low_RT[i,2].hires_radiance.S[:,:] * N2
        ) / (N1 + N2)

        @views r_high.S[:,:] = (
            lsi.high_RT[i,1].hires_radiance.S[:,:] * N1 +
            lsi.high_RT[i,2].hires_radiance.S[:,:] * N2
        ) / (N1 + N2)

        ε_τ[i] = RadType(T, 1)
        @views ε_τ[i].S[1,:] = (r_low.S[1,:] - r_high.S[1,:])
        # Again divide by intensity to get relative errors
        @views ε_τ[i].S[1,1] /= r_high.S[1,1]


        # Calculate derivative of the error term w.r.t. any weighting function
        # if ε = (low_bin - high_bin) / high_bin = (low_bin/high_bin) - 1
        # then
        # ∂ε/∂x = 1/high_bin * (∂low_bin/∂x - ∂high_bin/∂x * (1 + ε)).
        # ONLY FOR THE I STOKES COMPONENT
        # For Q and U, the equations are
        # ε = low_bin - high_bin
        # so
        # ∂ε = ∂low_bin/∂x - ∂high_bin/∂x


        ∂ε_τ[i] = [RadType(T, 1) for l in 1:N_wf]
        wf_low = RadType(T, 1)
        wf_high = RadType(T, 1)

        for l in 1:N_wf

            @views wf_low.S[:,:] = (
                lsi.low_RT[i,1].hires_wfunctions[l].S[:,:] * N1 +
                lsi.low_RT[i,2].hires_wfunctions[l].S[:,:] * N2
            ) / (N1 + N2)

            @views wf_high.S[:,:] = (
                lsi.high_RT[i,1].hires_wfunctions[l].S[:,:] * N1 +
                lsi.high_RT[i,2].hires_wfunctions[l].S[:,:] * N2
            ) / (N1 + N2)

            ∂ε_τ[i][l].S[1,1] = 1/r_high.S[1,1] * (
                wf_low.S[1,1] - wf_high.S[1,1] * (1 + ε_τ[i].S[1,1])
            )

            #if i in bad_bins
            #    ∂ε_τ[i][l].S[1,1] = NaN
            #end


            if RadType == VectorRadiance
                ∂ε_τ[i][l].S[1,2:end] = wf_low.S[1,2:end] - wf_high.S[1,2:end]
            end

        end

        # This is the mean of all gas-tau within the full tau-bin
        # (both ξ bins)
        mean_τ[i] = mean(lsi.tau_gas[idx_τ])

    end

    l = 5

    for s in 1:N_stokes

        plot_x = Float64[]
        plot_y1 = Float64[]
        plot_y2 = Float64[]

        for i in used_τ_bins

            push!(plot_x, log(τ_gas[i,1]))
            push!(plot_y1, ∂ε_τ[i][l].S[1,s])
            push!(plot_y2, ε_τ[i].S[1,s])
        end

        #display(scatterplot(plot_x, plot_y1, grid=false, marker=:xcross, title="stokes q=$(s)"))
        #display(scatterplot(plot_x, plot_y2, grid=false, marker=:xcross, title="stokes q=$(s)"))
    end

    #=
        Create an interpolation object, so we can estimate the error ε for any τ. This is
        ONLY used for the "centering correction", where we adjust the errors in a given
        bin i, such that ε[(i,1)] and ε[(i,2)] will correspond to the same τ coordinate.
        It is important here to extrapolate linearly, since xA and xB will very likely be
        outside of the range of observed x (√ξ).

        This has to be created for each Stokes component separately, and we do this
        in an explicit fashion for I,(Q,U) - depending on the type of radiance.
    =#

    # This interpolation object now has the same dimension as the radiances, so this
    # should behave correctly whether we use ScalarRadiance or VectorRadiance
    itp_centering_ε = LinearInterpolation(
        log.([mean_τ[x] for x in used_τ_bins]),
        [ε_τ[x].S[:,:] for x in used_τ_bins],
        extrapolation_bc=Linear()
        )

    # Do the same for weighting functions, but loop over the weighting function index l
    itp_centering_∂ε_I = [
        LinearInterpolation(
            log.([mean_τ[x] for x in used_τ_bins]),
            [∂ε_τ[x][l].S[:,:] for x in used_τ_bins],
            extrapolation_bc=Linear()
            )
            for l in 1:N_wf
    ]

    # Equation A3 in O'Dell

    # xA and xB define the extent of the "interpolation rectangle"; 0 and 1 are
    # the bounds of valid values, so this should, by definition, cover all possible
    # x = √ξ. One could use xA = minimum(x), xB = maximum(x) as well, as Chris has
    # done in his code (as opposed to what's in the paper).
    xA = minimum(lsi.ξ_sqrt) #0.0
    xB = maximum(lsi.ξ_sqrt) #1.0

    ε_τ_A = Dict()
    ε_τ_B = Dict()
    ∂ε_τ_A = Dict()
    ∂ε_τ_B = Dict()

    for i in used_τ_bins

        idx_ξ1 = (lsi.tau_gas_bin_assignment .== i) .& (lsi.ξ_sqrt_bin_assignment .== 1)
        idx_ξ2 = (lsi.tau_gas_bin_assignment .== i) .& (lsi.ξ_sqrt_bin_assignment .== 2)

        ξ1_mean = mean(lsi.ξ_sqrt[idx_ξ1])
        ξ2_mean = mean(lsi.ξ_sqrt[idx_ξ2])

        log_τ_mean = [log(mean(lsi.tau_gas[idx_ξ1])), log(mean(lsi.tau_gas[idx_ξ2]))]

        # Perform the centering correction, this will adjust the computed errors
        # within a bin ε([i,j]) to be representative of common gas optical depth

        for j in [1,2]

            if haskey(ε, (i,j))

                center_adj_I = ε_τ[i].S[:,:] .- itp_centering_ε(log_τ_mean[j])
                ε[i,j].S[:,:] .+= center_adj_I

                for l in 1:N_wf

                    ∂_center_adj_I = ∂ε_τ[i][l].S[:,:] .- itp_centering_∂ε_I[l](log_τ_mean[j])
                    ∂ε[i,j][l].S[:,:] .+= ∂_center_adj_I

                end
            end
        end

        # Extrapolate errors to be at x = xA and x = xB to create the interpolation
        # rectangle.
        ε_τ_A[i] = RadType(T, 1)
        ε_τ_B[i] = RadType(T, 1)

        x_slope_A = (xA - ξ1_mean) / (ξ2_mean - ξ1_mean)
        x_slope_B = (xB - ξ1_mean) / (ξ2_mean - ξ1_mean)

        @views @. ε_τ_A[i].S[:,:] = (ε[(i,2)].S[:,:] - ε[(i,1)].S[:,:]) * x_slope_A +
            ε[(i,1)].S[:,:]
        @views @. ε_τ_B[i].S[:,:] = (ε[(i,2)].S[:,:] - ε[(i,1)].S[:,:]) * x_slope_B +
            ε[(i,1)].S[:,:]

        ∂ε_τ_A[i] = [RadType(T, 1) for l in 1:N_wf]
        ∂ε_τ_B[i] = [RadType(T, 1) for l in 1:N_wf]

        for l in 1:N_wf

            @views @. ∂ε_τ_A[i][l].S[:,:] = (∂ε[(i,2)][l].S[:,:] - ∂ε[(i,1)][l].S[:,:]) *
                x_slope_A + ∂ε[(i,1)][l].S[:,:]
            @views @. ∂ε_τ_B[i][l].S[:,:] = (∂ε[(i,2)][l].S[:,:] - ∂ε[(i,1)][l].S[:,:]) *
                x_slope_B + ∂ε[(i,1)][l].S[:,:]

        end

    end

    center_idx = lsi.low_RT[(first_τ_bin, 1)].optical_properties.spectral_window.spectral_idx
    edge_idx = lsi.low_RT_edge[(first_τ_bin, 1)].optical_properties.spectral_window.spectral_idx

    ww_c = lsi.monochromatic_RT.optical_properties.spectral_window.ww_grid[center_idx]
    ww_edge = lsi.monochromatic_RT.optical_properties.spectral_window.ww_grid[edge_idx]

    #=
        Construct nodes for the interpolant

        NOTE:
        Just like with the more recent version in ACOS, we implement two different
        interpolants: one in log(τ_gas) and one in τ_gas. Once τ_gas is above a certain
        threshold, interpolation in τ_gas is used, otherwise we stick with the original
        log(τ_gas) formulation.

    =#

    ε_ξ_nodes_log = ([log(mean_τ[i]) for i in used_τ_bins], [xA, xB])
    ε_ξ_nodes = ([mean_τ[i] for i in used_τ_bins], [xA, xB])

    # Construct array for Stokes components
    ε_τ_AB = [[ε_τ_A[i], ε_τ_B[i]] for i in used_τ_bins]
    ε_int =[[ε_τ_AB[i][j].S[1,s] for i in 1:length(used_τ_bins), j in [1,2]]
            for s in 1:N_stokes]

    # Construct array for derivative components
    # which you then can access via ∂ε_int[l] with l being the weighting function index..
    ∂ε_int = []
    for l in 1:N_wf
        ∂ε_τ_AB = [[∂ε_τ_A[i][l], ∂ε_τ_B[i][l]] for i in used_τ_bins]
        this_int = [[∂ε_τ_AB[i][j].S[1,s] for i in 1:length(used_τ_bins), j in [1,2]]
            for s in 1:N_stokes]
        push!(∂ε_int, this_int)
    end

    ε_itp_log = [Interpolations.extrapolate(
        Interpolations.interpolate(ε_ξ_nodes_log, ε_int[s], Gridded(Linear())),
        Flat()
    ) for s in 1:N_stokes]

    ε_itp = [Interpolations.extrapolate(
        Interpolations.interpolate(ε_ξ_nodes, ε_int[s], Gridded(Linear())),
        Flat()
    ) for s in 1:N_stokes]

    ∂ε_itp_log = []
    ∂ε_itp = []
    for l in 1:N_wf

        this_itp = [Interpolations.extrapolate(
            Interpolations.interpolate(ε_ξ_nodes_log, ∂ε_int[l][s], Gridded(Linear())),
            Flat()
            ) for s in 1:N_stokes]
        push!(∂ε_itp_log, this_itp)

        this_itp = [Interpolations.extrapolate(
            Interpolations.interpolate(ε_ξ_nodes, ∂ε_int[l][s], Gridded(Linear())),
            Flat()
            ) for s in 1:N_stokes]
        push!(∂ε_itp, this_itp)

    end
    #=
        Slope correction terms (which are not a function of wavelength)
    =#

    error_center = RadType(T, 1)
    error_edge = RadType(T, 1)

    @views error_center.S[1,:] = (
        lsi.low_RT[(first_τ_bin,1)].hires_radiance.S[1,:] -
        lsi.high_RT[(first_τ_bin,1)].hires_radiance.S[1,:]
        )

    @views error_edge.S[1,:] = (
        lsi.low_RT_edge[(first_τ_bin,1)].hires_radiance.S[1,:] -
        lsi.high_RT_edge[(first_τ_bin,1)].hires_radiance.S[1,:]
        )

    # Again divide all Stokes components by intensity (Eq. 3)
    error_center.S[1,1] /= lsi.high_RT[(first_τ_bin,1)].hires_radiance.S[1,1]
    error_edge.S[1,1] /= lsi.high_RT_edge[(first_τ_bin,1)].hires_radiance.S[1,1]

    ∂_error_center = [RadType(T, 1) for l in 1:N_wf]
    ∂_error_edge = [RadType(T, 1) for l in 1:N_wf]

    for l in 1:N_wf
        @views @. ∂_error_center[l][1,:] = (
            lsi.low_RT[(first_τ_bin,1)].hires_wfunctions[l].S[1,:] -
            lsi.high_RT[(first_τ_bin,1)].hires_wfunctions[l].S[1,:]
            )

        ∂_error_center[l][1,1] -= error_center.S[1,1] *
            lsi.high_RT[(first_τ_bin,1)].hires_wfunctions[l].S[1,1]
        ∂_error_center[l][1,1] /=
            lsi.high_RT[(first_τ_bin,1)].hires_radiance.S[1,1]

        @views @. ∂_error_edge[l][1,:] = (
            lsi.low_RT_edge[(first_τ_bin,1)].hires_wfunctions[l].S[1,:] -
            lsi.high_RT_edge[(first_τ_bin,1)].hires_wfunctions[l].S[1,:]
            )

        ∂_error_edge[l][1,1] -= error_edge.S[1,1] *
            lsi.high_RT_edge[(first_τ_bin,1)].hires_wfunctions[l].S[1,1]
        ∂_error_edge[l][1,1] /=
            lsi.high_RT_edge[(first_τ_bin,1)].hires_radiance.S[1,1]
    end

    #=
        Final, in-place correction of Stokes components for every spectral point.
    =#

    _lsi_correct_stokes!(
        lsi.monochromatic_RT.hires_radiance,
        lsi.monochromatic_RT.hires_wfunctions,
        lsi.monochromatic_RT.optical_properties.spectral_window.ww_grid,
        ww_c,
        ww_edge,
        lsi.tau_gas,
        lsi.ξ_sqrt,
        error_edge,
        error_center,
        ∂_error_edge,
        ∂_error_center,
        ε_itp,
        ε_itp_log,
        ∂ε_itp,
        ∂ε_itp_log
    )

    # After the correction of radiances and WFs has taken place, we must produce
    # the new Jacobians, derived from the weighting functions.
    orig_rt = lsi.monochromatic_RT
    for (i, sve) in enumerate(orig_rt.state_vector.state_vector_elements)
        if calculate_jacobian_before_isrf(sve)
            @debug "Re-calculating Jacobian for $(sve)"
            calculate_rt_jacobian!(orig_rt.hires_jacobians[sve], orig_rt, sve)
        end
    end

    =#

end

"""
    Specialized high-performance function to correct Stokes radiances. This is meant to be
    used within the `perform_LSI_correction!` function.
"""
function _lsi_correct_stokes!(
    radiance,
    wfunctions,
    ww_array,
    ww_center,
    ww_edge,
    τ_array,
    x_array,
    error_edge,
    error_center,
    ∂_error_edge,
    ∂_error_center,
    ε_interpolant_linear,
    ε_interpolant_log,
    ∂ε_interpolant_linear,
    ∂ε_interpolant_log,
)

    # Element type of radiance array (e.g. Float64)
    T = eltype(radiance)
    # Short-cut to the type of the radiance, allows us to construct new ones
    RadType = typeof(radiance).name.wrapper

    # Make some minor allocations ahead of time to be used inside the spectral loo.
    # Note that all these radiance-type arrays have one spectral element only, so further
    # below, they have to be
    stokes_error = RadType(T, 1)
    slope_adjust = RadType(T, 1)
    ∂_stokes_error = RadType(T, 1)
    ∂_slope_adjust = RadType(T, 1)

    N_hires = length(ww_array)
    N_wf = length(∂ε_interpolant_linear)
    N_stokes = size(stokes_error, 2)

    first = true

    for idx_spec in 1:N_hires

        # TODO Make this faster!

        τ = τ_array[idx_spec]
        log_τ = log(τ)
        x = x_array[idx_spec]
        ww = ww_array[idx_spec]

        # Performance issue: casting directly causes allocations
        # such as: stokes_error[1,:] = ε_interpolant_linear(τ, x)
        # It works, but is quite slow. Probably due to some needed conversion
        # into a new vector or something like that.

        for s in 1:N_stokes
            if τ > 100.0
                stokes_error.S[1,s] = ε_interpolant_linear[s](τ, x)
            else
                stokes_error.S[1,s] = ε_interpolant_log[s](log_τ, x)
            end
        end

        #=
            Perform the slope correction (Eq. 6 in O'Dell)
        =#

        for s in 1:N_stokes
            slope_adjust.S[1,s] = (error_edge.S[1,s] - error_center.S[1,s])
            slope_adjust.S[1,s] *= (ww_center - ww) / (ww_center - ww_edge)

            stokes_error.S[1,s] += slope_adjust.S[1,s]
        end

        # In-place correction of low-stream radiances
        radiance.S[idx_spec,1] /= (1.0 + stokes_error.S[1,1])

        # Correct Vector radiance if needed
        if RadType == VectorRadiance
            radiance.S[idx_spec,2] -= stokes_error.S[1,2]
            radiance.S[idx_spec,3] -= stokes_error.S[1,3]
        end

        # .. do the same for weighting functions
        for l in 1:N_wf

            # Calculate interpolated error
            for s in 1:N_stokes
                if τ > 100.0
                    ∂_stokes_error.S[1,s] = ∂ε_interpolant_linear[l][s](τ, x)
                else
                    ∂_stokes_error.S[1,s] = ∂ε_interpolant_log[l][s](log_τ, x)
                end
            end

            # Slope-correct
            for s in 1:N_stokes
                ∂_slope_adjust.S[1,s] = (∂_error_edge[l].S[1,s] - ∂_error_center[l].S[1,s])
                ∂_slope_adjust.S[1,s] *= (ww_center - ww) / (ww_center - ww_edge)

                ∂_stokes_error.S[1,s] += ∂_slope_adjust.S[1,s]
            end

            wfunctions[l].S[idx_spec,1] -= radiance.S[idx_spec,1] * ∂_stokes_error.S[1,1]
            wfunctions[l].S[idx_spec,1] /= (1 + stokes_error.S[1,1])

            if RadType == VectorRadiance
                wfunctions[l].S[idx_spec,2] -= ∂_stokes_error.S[1,2]
                wfunctions[l].S[idx_spec,3] -= ∂_stokes_error.S[1,3]
            end

        end


    end

end