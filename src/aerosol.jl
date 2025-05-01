


"""
Reads a pair of MIE and MOM text files to create an appropriate aerosol property
object.

$(TYPEDSIGNATURES)

# Details
We expect the format of the Colorado State University MIE and MOM files. TODO Add more
documentation about Siewert vs. de Rooij conventions here.
"""
function read_mie_mom_file(
    mie_fname::String,
    mom_fname::String,
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits};
    convention::NTuple{6, Symbol}=(:a1, :a2, :a3, :a4, :b1, :b2)
    )

    mie_raw = readlines(mie_fname)
    mom_raw = readlines(mom_fname)

    #=
        Process the MIE file first
    =#

    # Remove lines that start with "#"
    idx_comment = [i for i in 1:length(mie_raw) if lstrip(mie_raw[i])[1] == '#']
    deleteat!(mie_raw, idx_comment)

    N_ww = length(mie_raw)

    ww = zeros(N_ww)
    Q_ext = zeros(N_ww)
    Q_sca = zeros(N_ww)
    omega = zeros(N_ww)

    # Do this slowly (I know a one-liner, but let's not go crazy..)
    # One-liner..
    #mie_numeric = reduce(hcat, [parse.(Ref(Float64), x) for x in split.(mie_raw)])'

    # Split each string into substrings, separated by spaces
    split_mie_str  = [split(x) for x in mie_raw]
    # Parse into Float64s
    split_mie_num = [parse.(Ref(Float64), x) for x in split_mie_str]
    # Turn into a matrix
    mie_mat = reduce(hcat, split_mie_num)'

    #=
        Fill our value arrays
    =#

    @views ww[:] = mie_mat[:, 1]
    @views Q_ext[:] = mie_mat[:, 2]
    @views Q_sca[:] = mie_mat[:, 3]
    @views omega[:] = mie_mat[:, 4]

    @debug "[AEROSOL] $(N_ww) spectral points in $(mie_fname): $(ww)"


    #=
        Process the MOM file second
    =#

    idx_comment = [i for i in 1:length(mom_raw) if lstrip(mom_raw[i])[1] == '#']
    deleteat!(mom_raw, idx_comment)

    # We have to loop through each section of the mom file, representing a phase function
    # expansion at a given wavelength.

    this_mom = []
    all_mom = []
    for (i_line, line) in enumerate(mom_raw)
        #@info "$(i_line) out of $(length(mom_raw))"
        sline = split(line)
        if length(sline) == 2
            # New section
            this_ww = parse(Float64, sline[1])
            this_ncoef = parse(Int, sline[2])

            if length(this_mom) != 0
                # Turn list of vectors into a matrix
                this_mom_matrix = reduce(hcat, this_mom)'
                push!(all_mom, this_mom_matrix)
            end

            # Reset placeholder list
            empty!(this_mom)

        elseif length(sline) == 6
            push!(this_mom, [parse(Float64, x) for x in sline])
        else
            error("[AEROSOL] Unexpected line in MOM file")
        end

        # Take care of the very last section
        if i_line == length(mom_raw)
            # Turn list of vectors into a matrix
            this_mom_matrix = reduce(hcat, this_mom)'
            push!(all_mom, this_mom_matrix)
        end

    end

    max_coefs = maximum([size(x, 1) for x in all_mom])
    coefficients = zeros(max_coefs, 6, N_ww)
    work = similar(coefficients)

    for (i_ww, this_mom) in enumerate(all_mom)
        this_n_coef = size(this_mom, 1)
        @views coefficients[1:this_n_coef, :, i_ww] = this_mom
    end

    # Reverse sign for γ and ε to match the Siewert convention for XRTM
    # THIS IS IMPORTANT!
    @turbo coefficients[:,5,:] *= -1
    @turbo coefficients[:,6,:] *= -1

    return MieMomAerosolProperty(
        ww,
        ww_unit,
        Q_sca,
        Q_ext,
        omega,
        max_coefs,
        coefficients,
        work,
        convention
    )

end

"""
    Calculates the needed Ångstrom scaling to obtain optical depth from
"""
function calculate_angstrom_scaling(
    aer::AbstractAerosolType,
    ww,
    ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}
)

    return calculate_angstrom_scaling(
        aer.optical_property,
        ww,
        ww_unit
    )

end

function calculate_angstrom_scaling(
    prop::MieMomAerosolProperty,
    ww::Number,
    ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}
)

    # Convert input wavelength to the wavelength of the MieMomAerosolProperty
    ww_use = ustrip(prop.ww_unit, ww * ww_unit)

    # Find out which pre-calculated properties are closest in wavelength
    idx1, idx2 = sortperm(abs.(ww_use .- prop.ww))[1:2]
    ww1 = prop.ww[idx1]
    ww2 = prop.ww[idx2]

    # Check how far the requested ww is from the known ww's:
    if abs(ww_use - ww1) > 2 * abs(ww2 - ww1)
        @warn "[AEROSOL] Known aerosol propertes far away from requested spectral point: "
        @warn "[AEROSOL] $(ww1) / $(ww2) vs. $(ww_use)"
        @warn "[AEROSOL] Spurious results possible!"
    end

    # Calculate Ångstrom coefficient for extinction
    α_ext = -(
        log((prop.Q_ext[idx1]) / (prop.Q_ext[idx2])) /
        log((prop.ww[idx1]) / (prop.ww[idx2]))
    )
    # Calculate Ångstrom coefficient for scattering
    α_sca = -(
        log((prop.Q_sca[idx1]) / (prop.Q_sca[idx2])) /
        log((prop.ww[idx1]) / (prop.ww[idx2]))
    )

    return α_ext, α_sca, idx1

end

function calculate_aerosol_tau_at_all_ww!(
    opt::EarthAtmosphereOpticalProperties,
    aer::AbstractAerosolType,
    ref_tau_ext::Vector
    )

    # Dispatch to specific type
    # (can be different, depending on the aerosol property type)
    return calculate_aerosol_tau_at_all_ww!(
        opt, aer, aer.optical_property, ref_tau_ext
        )

end

function calculate_aerosol_tau_at_all_ww!(
    opt::EarthAtmosphereOpticalProperties,
    aer::AbstractAerosolType,
    prop::MieMomAerosolProperty,
    ref_tau_ext::Vector

)

    # For this function, aer.optical_property === prop, we demand this, otherwise
    # the calculations below really don't make sensible.
    # Also, we have to make the `prop` argument explicit so we can dispatch
    # to the correct MieMomAerosolProperty type hier.

    @assert aer.optical_property === prop "Optical properties inside $(aer) are not the" *
        " same as the supplied ones!"

     # We need suitable Ångstrom coefficients for the reference wavelength
     α_ext_ref, α_sca_ref, idx_ref = calculate_angstrom_scaling(
        prop, aer.ww_reference, aer.ww_unit
    )

    # Reference wavelength/wavenumber in units of the optical property
    # wavelength/wavenumber
    ww_ref = ustrip(prop.ww_unit, aer.ww_reference * aer.ww_unit)

    # Obtain the Q_sca and Q_ext at the reference wavelength via Ångstrom ansatz
    Q_sca_ref = prop.Q_sca[idx_ref] /
        (ww_ref / prop.ww[idx_ref])^(-α_sca_ref)
    Q_ext_ref = prop.Q_ext[idx_ref] /
        (ww_ref / prop.ww[idx_ref])^(-α_ext_ref)

    # Calculating the optical depth profiles, we need (potentially) another set of
    # Ångstrom coefficients for the wavelengths of the band.
    α_ext, α_sca, idx = calculate_angstrom_scaling(
        prop, opt.spectral_window.ww_grid[1], opt.spectral_window.ww_unit
        )

    # Determine the per-wavelength scale factors needed to get the correct AOD
    @turbo for i in 1:opt.spectral_window.N_hires

        # Calculate Qsca and Qext for this wavelength
        Q_ext = aer.optical_property.Q_ext[idx] *
            (opt.spectral_window.ww_grid[i] / prop.ww[idx])^(-α_ext)
        Q_sca = aer.optical_property.Q_sca[idx] *
            (opt.spectral_window.ww_grid[i] / prop.ww[idx])^(-α_sca)

        for l in 1:length(ref_tau_ext)
            # Compute the total extinction
            opt.aerosol_tau[aer][i,l] = Q_ext / Q_ext_ref * ref_tau_ext[l]
            # Get the single-scatter albedo
            opt.aerosol_omega[aer][i,l] = (Q_sca / Q_sca_ref) / (Q_ext / Q_ext_ref) *
                prop.omega[idx_ref]
        end

        #=
            NOTE

            Note that the single-scatter albedo opt.aerosol_omega for a given aerosol
            does not change with layer! We keep the (spectral, layer)-shape of the
            aerosol ω arrays to make array calculations easier later on, but there is
            redundancy here.
            Any aerosol-type has only one ω for a given wavelength, since the optical
            depths due to extinction and scattering only vary with wavelength or
            wavenumber. Total (aerosol) ω are expected to change with layer when there
            is more than one aerosol type present, or gases also contribute to
            extninction.

        =#

    end

    if any(@views opt.aerosol_omega[aer][:,:] .> 1)
        @warn "[AEROSOL] Aerosol single-scatter albedo > 1 for $(aer)"
    end

end

"""
    Calculates and sets the linearized inputs for XRTM for aerosol optical depth. Keep in
    mind that this function will be called for every spectral point in the hires
    calculation! Note that this function calculates the inputs needed for ∂I/∂AOD,
    and not ∂I/∂log(AOD). Appropriate transformation of the resulting weighting function
    is performed elsewhere.
"""
function set_XRTM_wf!(
    xrtm,
    wf_idx::Integer,
    spectral_idx::Integer,
    rt::AbstractRTMethod,
    sve::AerosolOpticalDepthSVE,
    first_XRTM_call::Bool
)

    # Cast the current AOD value to a unitless quantity
    # (this can be either AOD or log(AOD))
    this_aod = get_current_value_with_unit(sve) |> NoUnits
    if sve.log
        this_aod = exp(this_aod)
    end

    # Useful short-cuts
    opt = rt.optical_properties
    this_aer_tau = opt.aerosol_tau[sve.aerosol]
    this_aer_ssa = opt.aerosol_omega[sve.aerosol]
    this_ray_tau = opt.rayleigh_tau

    tmp_Nlay1 = opt.tmp_Nlay1
    tmp_Nlay2 = opt.tmp_Nlay2

    if Threads.nthreads() > 1

        # Use a thread-safe variant if there are multiple threads calling this..

        tmp_Nlay1 = similar(opt.tmp_Nlay1)
        tmp_Nlay2 = similar(opt.tmp_Nlay2)

    end



    # Make this calculation only once!
    if first_XRTM_call
        #=
            The number of coefficient elements going into XRTM will depend on whether XRTM
            was initilized in vector mode, or not. We can ask XRTM itself for that
            information. If the number of stokes inside XRTM == 1, use 1, otherwise 6.
            Note that we create a new array here with a different amount of coefficient
            elements, since we want to retain the possibility to have set up the retrieval
            in "Vector" mode, but do some RT calculates for scalar only.
        =#
        lcoef = zeros(
            size(rt.optical_properties.total_coef, 1),
            XRTM.get_n_stokes(xrtm) == 1 ? 1 : 6
        )

        # Placeholder for interpolated aerosol coefficients (need only 2 dimensions)
        this_aer_coef = @view rt.optical_properties.tmp_coef[:,:,1]
        this_aer_coef[:] .= 0

        if Threads.nthreads() > 1

            # Use a thread-safe variant if there are multiple threads calling this..
            this_aer_coef = zeros(
                size(rt.optical_properties.tmp_coef, 1),
                size(rt.optical_properties.tmp_coef, 2),
                1
            )

        end

        # Aerosol coeffs: expansion moments, matrix element, wavelength
        # We need the aerosol coefficients for this mixture, and this wavelength center
        # (this interpolation call takes <0.1 ms generally)
        swin = rt.optical_properties.spectral_window
        idx_aer = get_scattering_index(swin)
        # Interpolate coefficients for this aerosol for spectral point `swin.ww_grid[idx_aer]`
        interpolate_aer_coef!(this_aer_coef, sve.aerosol, swin.ww_grid[idx_aer], swin.ww_unit)

    end



    #=
        Pre-compute the linearized inputs for ∂τ/∂AOD and ∂ω/∂AOD
    =#

    #=
        Part 1: ∂τ/∂AOD
        ∂τ/∂AOD = τ_aerosol / AOD
        (this should hold true for all AbstractAerosols, regardless of shape)
    =#

    # Short-cut to pre-allocated vector with N_layer length
    ∂τ_∂AOD = tmp_Nlay1
    ∂τ_∂AOD[:] .= 0
    @turbo for l in axes(∂τ_∂AOD, 1)
        τ_aer = this_aer_tau[spectral_idx, l]
        ∂τ_∂AOD[l] = τ_aer / this_aod
    end
    # Set ∂τ/∂AOD linearized input for all layers
    XRTM.set_ltau_l_n1(
        xrtm,
        wf_idx-1, # Derivative input (0-based for XRTM)
        ∂τ_∂AOD
    )

    #=
        Part 2: ∂ω/∂AOD
        ∂ω/∂AOD =  τ_aerosol / AOD * (ω_aerosol - ω) / τ
        (this should hold true for all AbstractAerosols, regardless of shape)
    =#

    # Short-cut to pre-allocated vector with N_layer length
    ∂ω_∂AOD = tmp_Nlay2
    ∂ω_∂AOD[:] .= 0
    @turbo for l in axes(∂ω_∂AOD, 1)
        τ = opt.total_tau[spectral_idx,l]
        ω_aer = this_aer_ssa[spectral_idx, l]
        ω = opt.total_omega[spectral_idx,l]

        ∂ω_∂AOD[l] = ∂τ_∂AOD[l] * (ω_aer - ω) / τ
    end
    # Set ∂ω/∂AOD linearized input for all layers
    XRTM.set_omega_l_n1(
        xrtm,
        wf_idx-1,
        ∂ω_∂AOD
    )

    if first_XRTM_call # Set this only once for bin center values

        for l in 1:rt.scene.atmosphere.N_layer

            @views lcoef[:,:] .= 0
            denom = 0.0
            denom += this_ray_tau[idx_aer, l]

            # This must be the sum over all aerosol optical thicknesses multiplied by the
            # aerosol ω, not just the one whose Jacobian we want! The aerosol single-
            # scatter albedo `ssa` must also be the correct one, not the one tied to this
            # state vector element.
            for (aer_key, aer_tau) in opt.aerosol_tau # Loops through all aerosols
                ssa = opt.aerosol_omega[aer_key] # Grabs the SSA from this aerosol
                denom += ssa[idx_aer, l] * aer_tau[idx_aer, l] # Adds SSA * τ_aer
            end

            # (1): ∂τ/∂AOD * ω_aerosol * (- β_total) / (τ_ray + Σ ω_aerosol * τ_aerosol)
            @turbo for c in axes(lcoef, 1) # Loops through pfmoms
                for q in axes(lcoef, 2) # Loops through elements (1 or 6)

                    lcoef[c,q] -= ∂τ_∂AOD[l] * this_aer_ssa[idx_aer, l] * (
                        rt.optical_properties.total_coef[c,q,l]) / denom

                end
            end

            # (2): ∂τ/∂AOD * ω_aerosol * (β_aerosol) / (τ_ray + Σ ω_aerosol * τ_aerosol)
            @turbo for c in axes(this_aer_coef, 1) # Loops through pfmoms
                for q in axes(lcoef, 2) # Loops through elements (1 or 6)

                    lcoef[c,q] += ∂τ_∂AOD[l] * this_aer_ssa[idx_aer, l] * (
                        this_aer_coef[c,q]) / denom

                end
            end

            # Reminder: linearized coefficient inputs for XRTM *must match* the number
            # of coefficients *used* (not necessarily supplied).
            # You cannot have N used coefficients, but then only vary M != N.
            ncoef = XRTM.get_n_coef(xrtm, l-1)

            XRTM.set_coef_l_11(
                xrtm,
                l-1,
                wf_idx-1,
                lcoef[1:ncoef,:]
            )

        end
    end

end


"""
$(TYPEDSIGNATURES)

Calculates and sets the linearized inputs for XRTM for aerosol height. Keep in mind that
this function will be called for every spectral point in the hires calculation! Note that
this function calculates the inputs needed for ∂I/∂height, and not ∂I/∂log(height).
Appropriate transformation of the resulting weighting function is performed elsewhere.
"""
function set_XRTM_wf!(
    xrtm,
    wf_idx::Integer,
    spectral_idx::Integer,
    rt::AbstractRTMethod,
    sve::AerosolHeightSVE,
    first_XRTM_call::Bool
)


    # Useful short-cuts
    opt = rt.optical_properties
    this_aer_tau = opt.aerosol_tau[sve.aerosol]
    this_aer_ssa = opt.aerosol_omega[sve.aerosol]
    this_ray_tau = opt.rayleigh_tau

    tmp_Nlay1 = opt.tmp_Nlay1
    tmp_Nlay2 = opt.tmp_Nlay2

    if Threads.nthreads() > 1

        # Use a thread-safe variant if there are multiple threads calling this..

        tmp_Nlay1 = similar(opt.tmp_Nlay1)
        tmp_Nlay2 = similar(opt.tmp_Nlay2)

    end

    # Cast the current height value to a unitless quantity
    # (this can be either height or log(height))
    this_height = get_current_value_with_unit(sve) |> NoUnits
    if sve.log
        this_height = exp(this_height)
    end

    if first_XRTM_call
        # The number of coefficient elements going into XRTM will depend on whether XRTM
        # was initilized in vector mode, or not. We can ask XRTM itself for that information.
        # If the number of stokes inside XRTM == 1, use 1, otherwise 6.
        n_stokes = XRTM.get_n_stokes(xrtm) == 1 ? 1 : 6
        lcoef = zeros(
            size(rt.optical_properties.total_coef, 1),
            n_stokes
        )

        # Placeholder for interpolated aerosol coefficients
        this_aer_coef = @view rt.optical_properties.tmp_coef[:,:,1]
        @views this_aer_coef[:] .= 0
        # Aerosol coeffs: expansion moments, matrix element, wavelength

        if Threads.nthreads() > 1

            # Use a thread-safe variant if there are multiple threads calling this..
            this_aer_coef = zeros(
                size(rt.optical_properties.tmp_coef, 1),
                size(rt.optical_properties.tmp_coef, 2),
                1
            )

        end


        # We need the aerosol coefficients for this mixture, and this wavelength center
        # (this interpolation call takes <0.1 ms generally)
        swin = rt.optical_properties.spectral_window
        idx_aer = get_scattering_index(swin)
        interpolate_aer_coef!(this_aer_coef, sve.aerosol, swin.ww_grid[idx_aer], swin.ww_unit)
    end


    #=
        Pre-compute the linearized inputs for ∂τ/∂height and ∂ω/∂height
    =#

    #=
        Part 1: ∂τ/∂height
    =#


    # Short-cut to pre-allocated vector with N_layer length
    ∂τ_∂height = tmp_Nlay1
    # Dispatch to aerosol-type specific function!
    for l in axes(∂τ_∂height, 1)
        ∂τ_∂height[l] = calculate_layer_dtau_dheight(
            rt.scene.atmosphere,
            sve.aerosol,
            l
        )
    end

    # Set ∂τ/∂height linearized input for all layers
    XRTM.set_ltau_l_n1(
        xrtm,
        wf_idx-1, # Derivative input (0-based for XRTM)
        ∂τ_∂height
    )

    #=
        Part 2: ∂ω/∂height for layer l
        ∂ω/∂height =  ∂τ/∂height * (ω_aerosol - ω) / τ
        (this should hold true for all AbstractAerosols, regardless of shape)
    =#

    # Short-cut to pre-allocated vector with N_layer length
    ∂ω_∂height = tmp_Nlay2

    @turbo for l in axes(∂ω_∂height, 1)
        τ = rt.optical_properties.total_tau[spectral_idx,l]
        ω_aer = this_aer_ssa[spectral_idx, l]
        ω = rt.optical_properties.total_omega[spectral_idx,l]

        ∂ω_∂height[l] = ∂τ_∂height[l] * (ω_aer - ω) / τ
    end

    # Set ∂ω/∂height linearized input for all layers
    XRTM.set_omega_l_n1(
        xrtm,
        wf_idx-1,
        ∂ω_∂height
    )

    if first_XRTM_call # Set this only once for bin center values
        for l in 1:rt.scene.atmosphere.N_layer

            @views lcoef[:,:] .= 0

            denom = 0.0
            denom += this_ray_tau[idx_aer, l]

            # This must be the sum over all aerosol optical thicknesses multiplied by the
            # aerosol ω, not just the one whose Jacobian we want! The aerosol single-
            # scatter albedo `ssa` must also be the correct one, not the one tied to this
            # state vector element.
            for (aer_key, aer_tau) in opt.aerosol_tau
                ssa = opt.aerosol_omega[aer_key]
                denom += ssa[idx_aer, l] * aer_tau[idx_aer, l]
            end

            # (1): ∂τ/∂height * ω_aerosol * (- β_total) /
            #                                (τ_ray + Σ ω_aerosol * τ_aerosol)
            @turbo for c in axes(lcoef, 1)
                for q in axes(lcoef, 2)

                    lcoef[c,q] -= ∂τ_∂height[l] * this_aer_ssa[idx_aer, l] * (
                        opt.total_coef[c,q,l]) / denom

                end
            end

            # (2): ∂τ/∂height * ω_aerosol * (β_aerosol) /
            #                               (τ_ray + Σ ω_aerosol * τ_aerosol)
            @turbo for c in axes(this_aer_coef, 1)
                for q in axes(lcoef, 2)

                    lcoef[c,q] += ∂τ_∂height[l] * this_aer_ssa[idx_aer, l] * (
                        this_aer_coef[c,q]) / denom

                end
            end

            # Reminder: linearized coefficient inputs for XRTM *must match* the number
            # of coefficients *used* (not necessarily supplied).
            # You cannot have N used coefficients, but then only vary M < N.
            ncoef = XRTM.get_n_coef(xrtm, l-1)

            XRTM.set_coef_l_11(
                xrtm,
                l-1,
                wf_idx-1,
                lcoef[1:ncoef,:]
            )

        end
    end

end




"""
$(TYPEDSIGNATURES)

Calculates and sets the linearized inputs for XRTM for aerosol height. Keep in mind that
this function will be called for every spectral point in the hires calculation! Note that
this function calculates the inputs needed for ∂I/∂width, and not ∂I/∂log(width).
Appropriate transformation of the resulting weighting function is performed elsewhere.
"""
function set_XRTM_wf!(
    xrtm,
    wf_idx::Integer,
    spectral_idx::Integer,
    rt::AbstractRTMethod,
    sve::AerosolWidthSVE,
    first_XRTM_call::Bool
)


    # Useful short-cuts
    opt = rt.optical_properties
    this_aer_tau = opt.aerosol_tau[sve.aerosol]
    this_aer_ssa = opt.aerosol_omega[sve.aerosol]
    this_ray_tau = opt.rayleigh_tau

    tmp_Nlay1 = opt.tmp_Nlay1
    tmp_Nlay2 = opt.tmp_Nlay2

    if Threads.nthreads() > 1

        # Use a thread-safe variant if there are multiple threads calling this..

        tmp_Nlay1 = similar(opt.tmp_Nlay1)
        tmp_Nlay2 = similar(opt.tmp_Nlay2)

    end

    # Cast the current width value to a unitless quantity
    # (this can be either width or log(width))
    this_width = get_current_value_with_unit(sve) |> NoUnits
    if sve.log
        this_width = exp(this_width)
    end

    if first_XRTM_call
        # The number of coefficient elements going into XRTM will depend on whether XRTM
        # was initilized in vector mode, or not. We can ask XRTM itself for that information.
        # If the number of stokes inside XRTM == 1, use 1, otherwise 6.
        n_stokes = XRTM.get_n_stokes(xrtm) == 1 ? 1 : 6
        lcoef = zeros(
            size(rt.optical_properties.total_coef, 1),
            n_stokes
        )

        # Placeholder for interpolated aerosol coefficients
        this_aer_coef = @view rt.optical_properties.tmp_coef[:,:,1]
        @views this_aer_coef[:] .= 0
        # Aerosol coeffs: expansion moments, matrix element, wavelength

        if Threads.nthreads() > 1

            # Use a thread-safe variant if there are multiple threads calling this..
            this_aer_coef = zeros(
                size(rt.optical_properties.tmp_coef, 1),
                size(rt.optical_properties.tmp_coef, 2),
                1
            )

        end


        # We need the aerosol coefficients for this mixture, and this wavelength center
        # (this interpolation call takes <0.1 ms generally)
        swin = rt.optical_properties.spectral_window
        idx_aer = get_scattering_index(swin)
        interpolate_aer_coef!(this_aer_coef, sve.aerosol, swin.ww_grid[idx_aer], swin.ww_unit)
    end


    #=
        Pre-compute the linearized inputs for ∂τ/∂width and ∂ω/∂width
    =#

    #=
        Part 1: ∂τ/∂width
    =#


    # Short-cut to pre-allocated vector with N_layer length
    ∂τ_∂width = tmp_Nlay1
    # Dispatch to aerosol-type specific function!
    for l in axes(∂τ_∂width, 1)
        ∂τ_∂width[l] = calculate_layer_dtau_dwidth(
            rt.scene.atmosphere,
            sve.aerosol,
            l
        )
    end

    # Set ∂τ/∂width linearized input for all layers
    XRTM.set_ltau_l_n1(
        xrtm,
        wf_idx-1, # Derivative input (0-based for XRTM)
        ∂τ_∂width
    )

    #=
        Part 2: ∂ω/∂width for layer l
        ∂ω/∂width =  ∂τ/∂width * (ω_aerosol - ω) / τ
        (this should hold true for all AbstractAerosols, regardless of shape)
    =#

    # Short-cut to pre-allocated vector with N_layer length
    ∂ω_∂width = tmp_Nlay2

    @turbo for l in axes(∂ω_∂width, 1)
        τ = rt.optical_properties.total_tau[spectral_idx,l]
        ω_aer = this_aer_ssa[spectral_idx, l]
        ω = rt.optical_properties.total_omega[spectral_idx,l]

        ∂ω_∂width[l] = ∂τ_∂width[l] * (ω_aer - ω) / τ
    end

    # Set ∂ω/∂width linearized input for all layers
    XRTM.set_omega_l_n1(
        xrtm,
        wf_idx-1,
        ∂ω_∂width
    )

    if first_XRTM_call # Set this only once for band center values
        for l in 1:rt.scene.atmosphere.N_layer

            @views lcoef[:,:] .= 0

            denom = 0.0
            denom += this_ray_tau[idx_aer, l]

            # This must be the sum over all aerosol optical thicknesses multiplied by the
            # aerosol ω, not just the one whose Jacobian we want! The aerosol single-
            # scatter albedo `ssa` must also be the correct one, not the one tied to this
            # state vector element.
            for (aer_key, aer_tau) in opt.aerosol_tau
                ssa = opt.aerosol_omega[aer_key]
                denom += ssa[idx_aer, l] * aer_tau[idx_aer, l]
            end

            # (1): ∂τ/∂width * ω_aerosol * (- β_total) /
            #                                (τ_ray + Σ ω_aerosol * τ_aerosol)
            @turbo for c in axes(lcoef, 1)
                for q in axes(lcoef, 2)

                    lcoef[c,q] -= ∂τ_∂width[l] * this_aer_ssa[idx_aer, l] * (
                        opt.total_coef[c,q,l]) / denom

                end
            end

            # (2): ∂τ/∂width * ω_aerosol * (β_aerosol) /
            #                               (τ_ray + Σ ω_aerosol * τ_aerosol)
            @turbo for c in axes(this_aer_coef, 1)
                for q in axes(lcoef, 2)

                    lcoef[c,q] += ∂τ_∂width[l] * this_aer_ssa[idx_aer, l] * (
                        this_aer_coef[c,q]) / denom

                end
            end

            # Reminder: linearized coefficient inputs for XRTM *must match* the number
            # of coefficients *used* (not necessarily supplied).
            # You cannot have N used coefficients, but then only vary M < N.
            ncoef = XRTM.get_n_coef(xrtm, l-1)

            XRTM.set_coef_l_11(
                xrtm,
                l-1,
                wf_idx-1,
                lcoef[1:ncoef,:]
            )

        end
    end

end


