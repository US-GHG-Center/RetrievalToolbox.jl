"""
$(TYPEDSIGNATURES)

Creates an empty `EarthOpticalProperties` object based on a supplied spectral window
`swin`, a scene `scene`, a state vector `sv` and the type of radiance to be used
(`RadType`). In addition, a float number type `T` must be supplied.
"""
function create_empty_EarthAtmosphereOpticalProperties(
    swin::AbstractSpectralWindow,
    scene::EarthScene,
    sv::AbstractStateVector,
    RadType::Type{<:Radiance},
    T::Type{<:AbstractFloat}
)

    atm_elements = scene.atmosphere.atm_elements # Atmosphere elements short-cut
    N_layer = scene.atmosphere.N_layer

    # Dictionaries to keep track of gas optical depths and derivatives
    gas_tau = Dict{GasAbsorber{T}, Array{T, 2}}()
    gas_derivatives = Dict{GasAbsorber{T}, Dict{String, AbstractArray}}()
    # Dictionaries to keep track of aerosol optical depths and single-scatter albedo
    aerosol_tau = Dict{AbstractAerosolType, Array{T, 2}}()
    aerosol_omega = Dict{AbstractAerosolType, Array{T, 2}}()

    # Rayleigh extinction optical depth
    rayleigh_tau = zeros(T, swin.N_hires, N_layer)
    # Rayleigh extinction derivative
    rayleigh_deriv = zeros(T, swin.N_hires, N_layer)
    # Total optical depth and single-scatter albedo
    total_tau = zeros(T, swin.N_hires, N_layer)
    total_omega = zeros(T, swin.N_hires, N_layer)

    for atm in atm_elements

        # Gases
        if atm isa GasAbsorber
            # per-gas layer-resolved optical depth
            gas_tau[atm] = zeros(T, swin.N_hires, N_layer)

            # related partial derivatives w.r.t.
            # surface pressure, temperature, VMRs
            # -> zero-arrays are created later when
            #    looping through the state vector elements

            # Note: this only creates an empty Dict, and is
            # filled just below if needed.
            gas_derivatives[atm] = Dict{String, AbstractArray}()

            if sv isa RetrievalStateVector
                for _sve in sv.state_vector_elements

                    if _sve isa SurfacePressureSVE
                        gas_derivatives[atm]["dTau_dpsurf"] =
                            zeros(T, swin.N_hires)
                    end

                    if _sve isa GasLevelScalingFactorSVE
                        gas_derivatives[atm]["dTau_dVMR"] =
                            zeros(T, swin.N_hires, N_layer, 2)
                    end

                    if _sve isa GasVMRProfileSVE
                        gas_derivatives[atm]["dTau_dVMR"] =
                            zeros(T, swin.N_hires, N_layer, 2)
                    end

                    if _sve isa TemperatureOffsetSVE
                        gas_derivatives[atm]["dTau_dT"] =
                            zeros(T, swin.N_hires, N_layer)
                    end

                end
            end

        end

        # Aersols
        if atm isa AbstractAerosolType

            aerosol_tau[atm] = zeros(T, swin.N_hires, N_layer)
            aerosol_omega[atm] = zeros(T, swin.N_hires, N_layer)

        end

        # TODO: add other atmospheric elements

    end # End of atm-loop

    # Allocate dry and wet air molecule numbers
    nair_dry = zeros(T, N_layer)
    nair_wet = zeros(T, N_layer)

    # Allocate vector(s) to hold temporary values
    tmp_Nhi1 = zeros(T, swin.N_hires)
    tmp_Nhi2 = zeros(T, swin.N_hires)
    tmp_Nlay1 = zeros(T, N_layer)
    tmp_Nlay2 = zeros(T, N_layer)

    # Determine the maximal number of phasefunction expansion
    # coefficients found in all aerosols of this scene.
    max_coef = 1
    for atm in atm_elements
        if atm isa AbstractAerosolType
            if atm.optical_property.max_coefs > max_coef
                max_coef = atm.optical_property.max_coefs
            end
        end
    end

    # If there are no aerosols, but there is RayleighScattering, adjust here
    if findanytype(atm_elements, RayleighScattering)
        # Bump up to 3 coefficients
        max_coef = max(max_coef, 3)
    end

    # Allocate matrix to hold phasefunction expansion coefficient
    # for this spectral window.
    # NOTE/TODO
    # Currently, we only use one per band (no spectral variation)
    # NOTE/TODO
    # Expand this for polarized RT
    if (length(keys(aerosol_tau)) > 0) | findanytype(atm_elements, RayleighScattering)

        # Check if we are doing polarization or not
        if RadType === ScalarRadiance
            n_elem = 1
        elseif RadType === VectorRadiance
            n_elem = 6
        end

        total_coef = zeros(T, max_coef, n_elem, N_layer)
        tmp_coef = zeros(T, max_coef, n_elem, N_layer)
        tmp_coef_scalar = zeros(T, max_coef, 1, N_layer)
    else
        # No aerosols, no Rayleigh scattering => no need for coefficients
        total_coef = nothing
        tmp_coef = nothing
        tmp_coef_scalar = nothing
    end

    # Allocate optical property object
    return EarthAtmosphereOpticalProperties(
        swin,
        gas_tau,
        gas_derivatives,
        aerosol_tau,
        aerosol_omega,
        rayleigh_tau,
        rayleigh_deriv,
        total_tau,
        total_omega,
        total_coef,
        nair_dry,
        nair_wet,
        tmp_Nhi1,
        tmp_Nhi2,
        tmp_Nlay1,
        tmp_Nlay2,
        tmp_coef,
        tmp_coef_scalar
    )

end


function calculate_earth_optical_properties!(
    buf::EarthAtmosphereBuffer,
    state_vector::AbstractStateVector;
    N_sublayer = 10
    )


    for swin in keys(buf.optical_properties)
        calculate_earth_optical_properties!(
            buf.optical_properties[swin],
            buf.scene,
            state_vector,
            N_sublayer=N_sublayer
        )
    end

end

"""
Creates an **EarthAtmosphereOpticalProperties** object, given some
**EarthScene** and **AbstractSpectralWindow**.

$(SIGNATURES)

# Details

Earth scenes have optical absorption and scattering due to gases,
aerosols and Rayleigh scattering. At the moment, the partial derivative of surface
pressure w.r.t. optical depth at the bottom-most layer (∂psurf / ∂τ) is calculated
via finite differencing, as the analytic calculation proved error-prone. As such, a
finite perturbation parameter is required, which is set to a default value.

"""
function calculate_earth_optical_properties!(
    opt::EarthAtmosphereOpticalProperties,
    scene::EarthScene,
    state_vector::AbstractStateVector;
    N_sublayer = 10,
    psurf_perturb = 1.0u"Pa"
    )

    # Create some convenient short-hand variables
    swin = opt.spectral_window
    atm = scene.atmosphere

    # Determine, whether we want optical depth derivatives
    return_dpsurf = false
    return_dT = false
    return_dVMR = false

    if typeof(state_vector) == RetrievalStateVector
        for _sve in state_vector.state_vector_elements

            if _sve isa SurfacePressureSVE
                @debug "[OPT] Surface pressure retrieval: need dTau_dpsurf"
                return_dpsurf = true
            end

            if _sve isa GasLevelScalingFactorSVE
                @debug "[OPT] Gas scaling factor retrieval: need dTau_dVMR"
                return_dVMR = true
            end

            if _sve isa GasVMRProfileSVE
                @debug "[OPT] Gas profile retrieval: need dTau_dVMR"
                return_dVMR = true
            end

            if _sve isa TemperatureOffsetSVE
                @debug "[OPT] Temperature offset retrieval: need dTau_dT"
                return_dT = true
            end

        end
    end

    # Finite differencing for psurf jacobian.
    # NOTE
    # Analytic calculation of the surface pressure Jacobian was
    # implemented at some point, but taken out. It does not consume
    # much time or memory and makes the underlying function to
    # calculate optical depths a little more legible.
    # This might be changed in the future again.

    if return_dpsurf

        Δp = ustrip(atm.pressure_unit, psurf_perturb)
        # Adjust the lowest pressure level
        atm.pressure_levels[atm.N_level] += Δp

        # Calculate optical depth profiles due to gases
        calculate_gas_optical_depth_profiles!(
            opt, # This will me modified!
            atm,
            swin,
            N_sublayer=N_sublayer,
            return_dT=false,
            return_dVMR=false,
            do_only_last_layer=true,
        )

        if findanytype(atm.atm_elements, AbstractRayleighScattering)
            # Check if we have Rayleigh scattering as
            # an atmospheric element..

            # Calculate OD profiles with perturbed surface pressure
            calculate_rayleigh_optical_depth_profiles!(opt, scene)
        end


        # We store τ(psurf + Δp) for the bottom layer
        for this_gas in keys(opt.gas_derivatives)
            @. @views opt.gas_derivatives[this_gas]["dTau_dpsurf"][:] =
                opt.gas_tau[this_gas][:,atm.N_layer]
        end



        @. @views opt.rayleigh_derivatives[:,:] = opt.rayleigh_tau[:,:]

        # Restore the surface pressure
        atm.pressure_levels[atm.N_level] -= Δp

    end

    # Calculate optical depth profiles due to gases
    calculate_gas_optical_depth_profiles!(
        opt, # This will me modified!
        atm,
        swin,
        N_sublayer=N_sublayer,
        return_dT=return_dT,
        return_dVMR=return_dVMR
    )

    if findanytype(atm.atm_elements, AbstractRayleighScattering)
        # Check if we have Rayleigh scattering as
        # an atmospheric element..

        # Calculate OD profiles
        calculate_rayleigh_optical_depth_profiles!(opt, scene)
    end



    if return_dpsurf
        for this_gas in keys(opt.gas_derivatives)

            # Calculate τ(psurf + Δp) - τ(psurf)
            @views @. opt.gas_derivatives[this_gas]["dTau_dpsurf"][:] -=
                opt.gas_tau[this_gas][:,atm.N_layer]

            # Calculate [τ(psurf + Δp) - τ(psurf)] / Δp to get
            # ≈ ∂τ/∂psurf
            @views @. opt.gas_derivatives[this_gas]["dTau_dpsurf"][:] /= Δp
        end

        # Same for change of Rayleigh scattering due to pressure level shift
        @. @views opt.rayleigh_derivatives[:,:] -= opt.rayleigh_tau[:,:]
        @. @views opt.rayleigh_derivatives[:,:] ./ Δp
    end

    #=
        Calculate aerosol optical properties
    =#

    calculate_aerosol_optical_properties!(opt, scene)


    #=
        Calculate the total layer contributions
    =#

    @views opt.total_tau[:,:] .= 0.0

    # Add up total optical depth from contributions of each gas
    for (gas, gas_tau) in opt.gas_tau
        @views opt.total_tau[:,:] .+= gas_tau[:,:]
    end

    # Add up aerosol contributions
    for (aer, aer_tau) in opt.aerosol_tau
        @views opt.total_tau[:,:] .+= aer_tau[:,:]
    end

    # Add Rayleigh contributions (this will be zero if no
    # RayleighScattering is present..)
    @views opt.total_tau[:,:] .+= opt.rayleigh_tau[:,:]

    #=
        Calculate total ω
    =#

    @views opt.total_omega[:,:] .= 0.0

    # Start with Rayleigh
    @views opt.total_omega[:,:] .+= opt.rayleigh_tau[:,:]
    # Add each aerosol scattering optical depths
    for (aer, aer_tau) in opt.aerosol_tau
        @views opt.total_omega[:,:] .+= opt.aerosol_omega[aer][:,:] .* aer_tau[:,:]
    end
    # Divide by total extinction tau
    @views opt.total_omega[:,:] ./= opt.total_tau[:,:]



    # Calculate the total phase function expansion coefficients
    # (layer resolved)


    # Need these calculations only if we have atmospheric elements that
    # cause scattering..

    if findanytype(scene.atmosphere.atm_elements, AbstractRayleighScattering) |
        findanytype(scene.atmosphere.atm_elements, AbstractAerosolType)

        calculate_total_coef!(opt, scene.atmosphere)
    end

end

"""
Calculates wavelength and layer-resolved optical depths for gaseous
absorbers defined in the atmosphere object **atm**, for wavelengths
given in the spectral window object **swin**. The result is a Dict
where each key corresponds to a gas embedded in **atm**.

$(TYPEDSIGNATURES)

# Details

Gas optical depth calculation is performed by looping from the surface
up, going to the top of the atmosphere.

Each layer is subdivided into sub-layers, whose number can be supplied
via the **N_sublayer** keyword.

Gas concentrations are assumed to change linearly with pressure.
"""
function calculate_gas_optical_depth_profiles!(
    opt::EarthAtmosphereOpticalProperties,
    atm::EarthAtmosphere,
    swin::AbstractSpectralWindow;
    N_sublayer = 10,
    return_dT = false,
    return_dVMR = false,
    do_only_last_layer = false,
    )

    # N-point Gauss rule for the sub-layer integration.
    # Calculate the x_i and w_i only once,
    # and then cheaply scale them to the appropriate integration limits
    # later on, when needed.

    N_gauss = 3
    x_gauss, w_gauss = gauss(N_gauss)

    # Tau gas is zero'd out always!
    for gas in keys(opt.gas_tau)
        @views opt.gas_tau[gas] .= 0.0
    end

    # Air columns zero'd out as well!
    @views opt.nair_dry[:] .= 0.0
    @views opt.nair_wet[:] .= 0.0

    # Manually set the temperature perturbation size
    T_perturb = 10.0 # in [K]

    # Take out gas absorbers only, we don't want aerosols or others
    gases = [c for c in atm.atm_elements if c isa GasAbsorber]

    # We must only consider gas absorbers, which lie within
    # our spectral window. Note that this only impacts which
    # gas properties we calculate, the arrays that hold them
    # might still exist (even though they will be zero always).

    remove_idx = []
    for i in 1:length(gases)
        if (minimum(gases[i].spectroscopy.ww) > maximum(swin.ww_grid)) |
            (maximum(gases[i].spectroscopy.ww) < minimum(swin.ww_grid))

            push!(remove_idx, i)
            @debug "[OPT] Removing gas $(gases[i]) from optical property calculations."
        end
    end

    # Remove the references to the gases not in this window
    deleteat!(gases, remove_idx)


    @assert length(gases) > 0 "Need at least one gas in atmosphere"
    @assert N_sublayer > 0 "Number of sub-layers must be > 0"


    if return_dVMR
        for gas in keys(opt.gas_derivatives)
            @views opt.gas_derivatives[gas]["dTau_dVMR"][:,:,:] .= 0.0
        end
    end

    if return_dT
        for gas in keys(opt.gas_derivatives)
            @views opt.gas_derivatives[gas]["dTau_dT"][:,:] .= 0.0
        end
    end

    # Note that we do not zero out dTau_dpsurf here! (not needed)

    # Rebind for convenience
    p = atm.pressure_levels
    psurf = atm.pressure_levels[end]

    p_met = atm.met_pressure_levels
    T = atm.temperature_levels
    sh = atm.specific_humidity_levels
    grav = atm.gravity_levels
    wl = swin.ww_grid


    # Create interpolation objects to sample met profiles
    # at any pressure.
    # This operation is quick and fast, does not allocate much
    # much memory. By default, extrapolation (linear)
    # is activated, since the met profiles tend to not fully cover
    # the full pressure range used in the retrieval grid.

    T_int = linear_interpolation(p_met, T, extrapolation_bc = Line())
    sh_int = linear_interpolation(p_met, sh, extrapolation_bc = Line())
    grav_int = linear_interpolation(p_met, grav, extrapolation_bc = Line())

    spec_wl_idx_left = Dict()
    is_swin_matched_with_spec = Dict()
    gas_idx_map = Dict()

    # Speed-up possibility:
    # if the wavelength grid is matched with the spectroscopy grid,
    # we do not need to interpolate between spectroscopy grid points.
    for this_gas in gases

        this_spec = this_gas.spectroscopy

        @debug "[OPT] Calculating optical properties for $(this_gas)"

        @debug "[OPT] Checking if gas $(this_gas.gas_name) is matched with $(swin.window_name)."
        #=
        Check if the retrieval wavelength grid is an exact
        subset of the spectroscopy wavelength grid.
        Note that this does not yet account for skips,
        e.g. if the wavelength grid is the same as the spectroscopy one,
        but only skips every n'th point.
        =#

        # The unit stripping is required to ensure that both spectral window
        # and spectroscopy object are using the same wavelength units.

        start_idx = vector_contained_in_vector(
            ustrip.(
                Ref(swin.ww_unit),
                this_spec.ww * this_spec.ww_unit
            ),
            wl
        )

        if start_idx == -1
            # If not, we can find the wavelength indices that are
            # later used for wavelength interpolation
            @debug "[OPT] Gas $(this_gas) is NOT matched with the wavelength grid!"
            @debug "[OPT] Gas optical depth calculations will be slower."
            is_swin_matched_with_spec[this_gas] = false
            spec_wl_idx_left[this_gas] = _find_ww_indices(
                ustrip.(
                    Ref(swin.ww_unit),
                    this_spec.ww * this_spec.ww_unit
                ),
                wl
            )

            if findany(spec_wl_idx_left[this_gas], -1)
                @warn "[OPT] Spectroscopy for $(this_gas) out-of-bounds. " *
                    "Gas optical depth values for some spectral points will be zero."

            end

        else
            # If so, the wavelength indices will correspond
            # truly to the position of wavelengths in the spectroscopy
            # wavelength grid
            @debug "[OPT] Gas $(this_gas) is matched with the wavelength grid!"
            is_swin_matched_with_spec[this_gas] = true
            spec_wl_idx_left[this_gas] = collect(
                start_idx:start_idx+length(wl)-1)
        end

    end

    # Main layer loop
    # This loops from the bottom-most layer
    # up to to to top-most layer of the RT grid.
    if do_only_last_layer
        # Useful when calculating ∂τ/∂psurf: use only the last layer
        layer_iterator = atm.N_layer+1:-1:atm.N_layer
        layer_iterator = [atm.N_layer + 1]
    else
        layer_iterator =  atm.N_layer+1:-1:2
    end
    # "l" refers to a level index, so to speak
    # use "l-1" to index the corresponding layer
    for l in layer_iterator

        # "Lower" here refers to "closer to the surface"
        # "Higher" means closer to the top of the atmosphere
        # .. so p_lower > p_higher

        p_lower = p[l]
        p_higher = p[l-1]

        # Calculate the iterator for the sub-layer loop
        p_sub_iterator = LinRange(p_lower, p_higher, N_sublayer + 1)

        # Sub-layer loop
        # (we go lower -> higher, in altitude)
        for (k, _this_p) in enumerate(p_sub_iterator[2:end])

            # Calculate the x_i and w_i for the specific integration
            # interval required by this sublayer.

            p_gauss_scaled = @. (x_gauss + 1) / 2 * (
                p_sub_iterator[k] - p_sub_iterator[k+1]
            ) + p_sub_iterator[k+1]

            w_gauss_scaled = @. w_gauss / 2 * (
                p_sub_iterator[k] - p_sub_iterator[k+1]
            )

            # Remember, the integration within the sublayer is evaluated
            # as INTEGRAL += (this_w * INTEGRAND)
            for (i_sub, (this_p, this_w)) in enumerate(
                zip(p_gauss_scaled, w_gauss_scaled))

                this_p_fac = (this_p - p_higher) / (p_lower - p_higher)

                this_T = T_int(this_p)
                this_sh = sh_int(this_p)
                this_H2O = this_sh / (1 - this_sh) * MM_AIR_TO_H2O
                this_grav = grav_int(this_p)

                C_tmp = 1.0 / this_grav * ustrip(NA) / ustrip(MM_DRY_AIR)

                # Dry and wet (moist) air molecules per m^2
                # (and correct for units)
                unit_factor = ustrip(u"m^-2",
                                     1.0 * unit(NA) * atm.pressure_unit
                                     / (atm.gravity_unit * unit(MM_DRY_AIR))
                                     )

                opt.nair_dry[l-1] += (1 - this_sh) * C_tmp * this_w * unit_factor
                opt.nair_wet[l-1] += (1 + this_sh * (MM_AIR_TO_H2O - 1)) *
                    C_tmp * this_w * unit_factor

                C_tmp *= this_w #dp_sub

                # Loop over gases
                for this_gas in gases

                    # Calculate the VMR at this point in the sublayer
                    # Account for VMR units here, this is the ONLY part of this
                    # function where the GasAbsorber is queried for its vmr_levels.
                    VMR_lower = this_gas.vmr_levels[l] * this_gas.vmr_unit
                    VMR_higher = this_gas.vmr_levels[l-1] * this_gas.vmr_unit

                    # Here, the VMR is automatically cast to a unitless quantity
                    this_VMR = (1.0 - this_p_fac) * VMR_lower + this_p_fac * VMR_higher

                    # Rebind for convenience
                    spec = this_gas.spectroscopy

                    # Read absorption coefficient at all wavelengths "wl", for
                    # this partiular p, T, H2O

                    # zero out
                    @views @. opt.tmp_Nhi1[:] = 0.0

                    get_absorption_coefficient_value_at!(
                        opt.tmp_Nhi1, # Output!
                        spec,
                        wl,
                        spec_wl_idx_left[this_gas],
                        is_swin_matched_with_spec[this_gas],
                        this_p,
                        this_T,
                        this_H2O,
                    )

                    # At this point, opt.tmp_Nhi1 contains σ(p, T, q)

                    # Calculate the correct unit conversion factor, which all
                    # depends on the user-supplied units to the various quantities.
                    # This will raise an error if units do not match!

                    # NOTE
                    # We do not account for the gas VMR units here, since that is
                    # used directly in the lines below.
                    # Also, VMRs are ALWAYS assumed to by dry-air mixing ratios
                    unit_fac = 1.0 * (
                        spec.cross_section_unit * unit(NA) * atm.pressure_unit
                        / (atm.gravity_unit * unit(MM_DRY_AIR))
                    ) |> NoUnits


                    # Optical depth contribution is added
                    @views @. opt.gas_tau[this_gas][:, l-1] += (
                        opt.tmp_Nhi1[:] * this_VMR * C_tmp * (1 - this_sh) * unit_fac
                    )

                    ######################################
                    # Derivatives w.r.t. level VMR changes
                    ######################################
                    if return_dVMR
                        # "Higher" (i.e. closer to TOA)
                        @views @. opt.gas_derivatives[this_gas]["dTau_dVMR"][:,l-1,1] +=
                            opt.tmp_Nhi1[:] * this_p_fac * C_tmp *
                            (1 - this_sh) * unit_fac
                        # "Lower" (i.e. closer to surface)
                        @views @. opt.gas_derivatives[this_gas]["dTau_dVMR"][:,l-1,2] +=
                            opt.tmp_Nhi1[:] * (1 - this_p_fac) * C_tmp *
                            (1 - this_sh) * unit_fac
                    end


                    if return_dT

                        # zero out
                        @views @. opt.tmp_Nhi2[:] = 0.0

                        # Perturb temperature
                        get_absorption_coefficient_value_at!(
                            opt.tmp_Nhi2, # Output!
                            spec,
                            wl,
                            spec_wl_idx_left[this_gas],
                            is_swin_matched_with_spec[this_gas],
                            this_p,
                            this_T + T_perturb, # Add some small T perturbation
                            this_H2O,
                        )

                        @views @. opt.gas_derivatives[this_gas]["dTau_dT"][:,l-1] +=
                            opt.tmp_Nhi2[:] * this_VMR * C_tmp * (1 - this_sh) * unit_fac
                        end

                end # End gas loop

            end # End Gauss integration loop

        end # End sub-layer loop

    end # End layer loop

    if return_dT
        for this_gas in gases
            for l in 1:atm.N_layer
                # Calculate τ(T + ΔT) - τ(T)
                @views @. opt.gas_derivatives[this_gas]["dTau_dT"][:,l] -=
                    opt.gas_tau[this_gas][:, l]
                # Calculate [τ(T + ΔT) - τ(T)] / ΔT
                @views @. opt.gas_derivatives[this_gas]["dTau_dT"][:,l] /=
                    T_perturb
            end
        end
    end

end

"""
Lower-level implementation for Rayleigh optical
depth calculation. Improves speed a little bit
because Julia can infer types here.

Note that this calculation could further be
improved by using sub-layer integration for the
optical depth and making use of (usually) finer
vertically resolved meteorological profiles.
"""
function _calculate_rayleigh_optical_depth_profiles!(
    ray_tau,
    p_levels,
    p_layers,
    met_p_levels,
    T_levels,
    grav_levels,
    wavelength,
    )

    N_hires, N_layer = size(ray_tau)

    # Create interpolation objects for met profiles
    grav_int(x) = linear_interpolation(met_p_levels, grav_levels,
                                       extrapolation_bc = Line())(x)

    for j in 1:N_layer

        this_p = p_layers[j]
        delta_p = p_levels[j+1] - p_levels[j] # This has + sign, always

        this_grav = grav_int(this_p)

        @inbounds for i in 1:N_hires

            # Both wavelength and Rayleigh cross section have units
            wl = wavelength[i]

            # Rayleigh cross section is in cm^2
            ray_sigma = calculate_rayleigh_sigma(wl)

            # This quantity then comes out as unit-less and can
            # thus be inserted into the ray_tau array which does not
            # carry units.

            ray_tau[i,j] = 1.0 * (
                NA *
                ray_sigma * delta_p /
                (this_grav * MM_DRY_AIR)
            ) |> upreferred


        end
    end

end

"""
Calculates wavelength and layer-resolved optical depths due to
Rayleigh scattering for a full **EarthAtmosphereBuffer**, such
that all spectral windows inside are process accordingly.

$(SIGNATURES)

# Details

Add more details
"""
function calculate_rayleigh_optical_depth_profiles!(
    opt::EarthAtmosphereOpticalProperties,
    scene::EarthScene;
    )

    atm = scene.atmosphere

    # Create unit-ful quantities here so that the lower-lying function
    # can assess units and conversions correctly.

    T = atm.temperature_levels * atm.temperature_unit
    grav = atm.gravity_levels * atm.gravity_unit
    plev = atm.pressure_levels * atm.pressure_unit
    play = atm.pressure_layers * atm.pressure_unit
    met_plev = atm.met_pressure_levels * atm.met_pressure_unit

    # Dispatch to high-performance function
    _calculate_rayleigh_optical_depth_profiles!(
            opt.rayleigh_tau,
            plev,
            play,
            met_plev,
            T,
            grav,
            opt.spectral_window.ww_grid * opt.spectral_window.ww_unit
        )

end

"""
    Calculates the refractive index of air for a given wavelength `λ` according
    to Peck and Reeder,  J. Opt. Soc. Am. 62 (1972)
"""
function refractive_index_peck_reeder(λ::Unitful.Length)

    if λ < 0.185u"µm" | λ > 1.69u"µm"
        @debug "[OPT] Warning: Peck-Reeder formula for refractive index not in valid range. "
    end

    # Convert to microns, as required by the formula
    λ_mu = ustrip(u"µm", λ)

    λ_m_2 = 1.0 / (λ_mu * λ_mu)

    res1 = 8060.51 + 2480990.0 / (132.274 - λ_m_2) +
        17455.7 / (39.32957 - λ_m_2)

    res = res1 / 1e8 + 1.0

    return res

end

function create_refracted_sza(
    opt::EarthAtmosphereOpticalProperties,
    scene::EarthScene;
)

    atm = scene.atmosphere
    sza_per_layer = zeros(atm.N_layer)

    T_int = linear_interpolation(
        atm.met_pressure_levels,
        atm.temperature_levels,
        extrapolation_bc = Line()
        )

    for i in atm.N_layer:-1:1

        if i == atm.N_layer
            elev_deg = 90.0 - scene.solar_zenith
            this_sza = scene.solar_zenith
        else
            elev_deg = 90.0 - sza_per_layer[i+1]
            this_sza = sza_per_layer[i+1]
        end

        p_pa = ustrip(u"Pa", atm.pressure_levels[i] * atm.pressure_unit)
        T_K = T_int(atm.pressure_levels[i])

        Δ_sza = rad2deg(8.15e-7 / tan(deg2rad(elev_deg + 10.3 / (elev_deg + 5.11)))
            * p_pa / T_K)


        sza_per_layer[i] = this_sza - Δ_sza

    end

    return sza_per_layer

end

function create_sphericity_factors(
    opt::EarthAtmosphereOpticalProperties,
    scene::EarthScene;
    solar_only=true
)

    # See Dahlback and Stamnes, https://doi.org/10.1016/0032-0633(91)90061-E
    # (appendix B)

    factors_solar = ones(scene.atmosphere.N_layer)

    atm = scene.atmosphere

    # Calculate the per-layer SZA for a refracting atmosphere
    sza_per_layer = create_refracted_sza(opt, scene)

    # Create interpolation object for altitude
    # (we need altitudes for RT grid)
    altitude_int = linear_interpolation(
        ustrip.(Ref(atm.pressure_unit), atm.met_pressure_levels * atm.met_pressure_unit),
        atm.altitude_levels,
        extrapolation_bc = Line()
        )

    # Let's do this without allocations just by looping over all layers, do the
    # calculations and then loop over all gases to multiply the newly calculated
    # factors.

    r_p = ustrip(atm.altitude_unit, EARTH_RADIUS)

    for j in 1:atm.N_layer

        #=
            Solar path contributions
        =#

        alt_up = altitude_int(atm.pressure_levels[j])
        alt_lo = altitude_int(atm.pressure_levels[j+1])

        r_p_term_solar = (r_p + alt_lo)^2 * (sind(sza_per_layer[j]))^2

        # Altitude of the layer boundary above (closer to TOA)
        r_up = r_p + alt_up
        # Altitude of the layer boundary below (closer to ground)
        r_lo = r_p + alt_lo

        Δs = sqrt(r_up^2 - r_p_term_solar) - sqrt(r_lo^2 - r_p_term_solar)
        Δh = alt_up - alt_lo

        # Store sphericity factor
        factors_solar[j] = (Δs / Δh)

    end

    # Override .. for testing only
    @views factors_solar[:] .= 1 / cosd(scene.solar_zenith)

    if solar_only
        return factors_solar
    end

end

"""
Calculates the derivative ∂τ/∂height for a GaussAerosol at layer `l`.

$(TYPEDSIGNATURES)
"""
function calculate_layer_dtau_dheight!(
    result::AbstractArray,
    atm::EarthAtmosphere,
    aer::GaussAerosol
)

    #=
        Recall: let the shape be S(p) = exp(-(p-p0)^2 / (2σ^2)),
        and thus τ_aer = AOD * S(p) / ∑S(p'), note the normalization factor by which
        we must divide. Also note that ∂/∂p0 S(p) = (p-p0)/σ^2 * S(p)

        ∂τ_aer(p)/∂p0 = AOD * ∂/∂p0 [S(p) / ∑S(p')] = ... (some algebra) =
        =

    =#

    if aer.relative_pressure
        # Result is: (x - x0) * τ_aer / σ^2, where
        # x[l] = p[l] / psurf and x0 = p0 / psurf,
        # p0 being the layer height `aer.pressure`

        psurf = atm.pressure_levels[end]
        p0 = aer.pressure
        σ = aer.width
        pprime_fac = 1.0 / psurf

    else

        p0 = ustrip(atm.pressure_unit, aer.pressure * aer.pressure_unit)
        σ = ustrip(atm.pressure_unit, aer.width * aer_pressure_unit)
        pprime_fac = 1.0

    end

    _calculate_layer_dtau_dheight!(
        result,
        atm.pressure_layers,
        p0,
        σ,
        pprime_fac
    )

    @views result[:] .*= aer.total_optical_depth

end

function _calculate_layer_dtau_dheight!(
    result,
    play,
    p0,
    σ,
    pprime_fac
)

    Ssum = 0
    Ssum2 = 0

    for i in eachindex(play)
        pp = pprime_fac * play[i]
        Ssum += exp(-(pp - p0)^2 / (2*σ^2))
        Ssum2 += (pp - p0) / σ^2 * exp(-(pp - p0)^2 / (2*σ^2))
    end

    for i in eachindex(play)

        p = pprime_fac * play[i]
        S = exp(-(p-p0)^2 / (2*σ^2))

        result[i] = ((p - p0) / σ^2 * Ssum - Ssum2) * S
        result[i] /= Ssum^2

    end

end


"""
Calculates the derivative ∂τ/∂width for a GaussAerosol at layer `l`.

$(TYPEDSIGNATURES)
"""
function calculate_layer_dtau_dwidth!(
    result::AbstractArray,
    atm::EarthAtmosphere,
    aer::GaussAerosol
)

    #=
        Recall: let the shape be S(p) = exp(-(p-p0)^2 / (2σ^2)),
        and thus τ_aer = AOD * S(p) / ∑S(p'), note the normalization factor by which
        we must divide. Also note that ∂/∂σ S(p) = (p-p0)^2/σ^3 * S(p)

        ∂τ_aer(p)/∂σ = AOD * ∂/∂σ [S(p) / ∑S(p')] = ... (some algebra) =
        =

    =#

    if aer.relative_pressure
        # Result is: (x - x0) * τ_aer / σ^2, where
        # x[l] = p[l] / psurf and x0 = p0 / psurf,
        # p0 being the layer height `aer.pressure`

        psurf = atm.pressure_levels[end]
        p0 = aer.pressure
        σ = aer.width
        pprime_fac = 1.0 / psurf

    else

        p0 = ustrip(atm.pressure_unit, aer.pressure * aer.pressure_unit)
        σ = ustrip(atm.pressure_unit, aer.width * aer_pressure_unit)
        pprime_fac = 1.0

    end

    _calculate_layer_dtau_dwidth!(
        result,
        atm.pressure_layers,
        p0,
        σ,
        pprime_fac
    )

    @views result[:] .*= aer.total_optical_depth

end

function _calculate_layer_dtau_dwidth!(
    result,
    play,
    p0,
    σ,
    pprime_fac
)

    Ssum = 0
    Ssum2 = 0

    for i in eachindex(play)
        pp = pprime_fac * play[i]
        Ssum += exp(-(pp - p0)^2 / (2*σ^2))
        Ssum2 += (pp - p0)^2 / σ^3 * exp(-(pp - p0)^2 / (2*σ^2))
    end

    for i in eachindex(play)

        p = pprime_fac * play[i]
        S = exp(-(p-p0)^2 / (2*σ^2))

        result[i] = ((p - p0)^2 / σ^3 * Ssum - Ssum2) * S
        result[i] /= Ssum^2

    end

end


function calculate_aerosol_optical_properties!(
    opt::EarthAtmosphereOpticalProperties,
    scene::EarthScene
)

    for (aer, aer_tau) in opt.aerosol_tau
        @debug "[OPT] Calculating aerosol optical propertes for $(aer)"
        # Dispatches to specific aerosol type
        calculate_aerosol_optical_depth_profile!(opt, scene, aer)
    end

end


"""
    Calculates optical depth profile for a gaussian aerosol
"""
function calculate_aerosol_optical_depth_profile!(
    opt::EarthAtmosphereOpticalProperties,
    scene::EarthScene,
    aer::GaussAerosol
)


    if aer.relative_pressure

        # Height and width are given in fractions of surface pressure

        p0 = aer.pressure * scene.atmosphere.pressure_levels[end]
        w0 = aer.width * scene.atmosphere.pressure_levels[end]

    else

        # Height and width are in absolute pressure units

        # Convert aerosol pressure-related quantities to the same units as the
        # retrieval pressure grid.
        p0 = ustrip(scene.atmosphere.pressure_unit, aer.pressure * aer.pressure_unit)
        w0 = ustrip(scene.atmosphere.pressure_unit, aer.width * aer.pressure_unit)

    end


    # At the reference wavelength, the aerosol extinction profile is calculated with the
    # user-defined AOD.
    ref_tau_ext = zeros(scene.atmosphere.N_layer)
    p = scene.atmosphere.pressure_layers

    @turbo for l in 1:scene.atmosphere.N_layer

        ref_tau_ext[l] = exp(-(p[l] - p0)^2 / (2*w0*w0))

    end

    @views ref_tau_ext[:] ./= sum(ref_tau_ext)
    @views ref_tau_ext[:] .*= aer.total_optical_depth

    # Now dispatch to specific function to calculate the extinction and ω at all
    # wavelengths, starting from the reference profile at the reference wavelength!

    calculate_aerosol_tau_at_all_ww!(
        opt,
        aer,
        ref_tau_ext,
        )

    return nothing

end

"""
    Calculates the total per-layer phase function coefficients
"""
function calculate_total_coef!(
    opt::EarthAtmosphereOpticalProperties,
    atm::EarthAtmosphere
)
    # Skip all of this in case there is no appropriate matrix inside the
    # optical property container `opt`
    if isnothing(opt.total_coef)
        return nothing
    end

    swin = opt.spectral_window

    # Set to zero
    @views opt.total_coef[:,:,:] .= 0

    # Number of layers in this
    N_layer = size(opt.total_coef, 3)
    N_elem = size(opt.total_coef, 2)
    # Separate array for denominator
    coef_denom = similar(opt.total_coef)


    # Temp buffer for interpolated aerosol coefficients
    if findanytype(atm.atm_elements, AbstractAerosolType)
        aer_coef_tmp = @view opt.tmp_coef[:,:,1]
    end

    # Figure out if we have a RayleighScattering object that requires us to
    # calculate the Rayleigh phase function coefficients
    have_rayleigh = findanytype(atm.atm_elements, AbstractRayleighScattering)

    ray_coef = zeros(3, N_elem)

    for i in [get_scattering_index(swin)] # For now, only band center

        @views coef_denom[:,:,:] .= 0

        #=
            First step: calculate the aerosol terms of the numerator only, for each layer
        =#

        for (aer, aer_tau) in opt.aerosol_tau

            #=
                Interpolate the aerosol phase function coefficients for this wavelength
                Note that we are doing this calculation *inside* the aerosol type loop,
                since we cannot expect that all aerosols have optical properties with the
                same wavelength/wavenumber definitions.
            =#
            @views aer_coef_tmp[:,:] .= 0 # Re-set temp array
            interpolate_aer_coef!(aer_coef_tmp, aer, swin.ww_grid[i], swin.ww_unit)

            @turbo for l in 1:N_layer

                 # Loop over moments
                 for i1 in axes(aer_coef_tmp, 1)
                    # Loop over elements (1 for scalar, 6 for vector)
                    for i2 in axes(aer_coef_tmp, 2)

                        # Sum contributions
                        opt.total_coef[i1, i2, l] += (
                            aer_coef_tmp[i1, i2]
                            * opt.aerosol_omega[aer][i,l]
                            * opt.aerosol_tau[aer][i,l]
                        )

                    end # End element loop
                end # End moment loop

            end # End layer loop

        end # End aerosol type loop

        #=
            Second step: add Rayleigh components (if needed) to the numerator
        =#

        if have_rayleigh
            # Compute the Rayleigh scattering coefficients for this spectral point,
            # note that this function wants wavelength/wavenumber in physical units.
            create_rayleigh_coefs!(ray_coef, swin.ww_grid[i] * swin.ww_unit, N_elem == 6)

            @turbo for l in 1:N_layer
                for i1 in axes(ray_coef, 1) # Moments
                    for i2 in axes(ray_coef, 2) # Elements

                        opt.total_coef[i1, i2, l] +=
                            opt.rayleigh_tau[i,l] * ray_coef[i1, i2]

                    end # End element loop
                end # End moment loop
            end # End layer loop

        end # End Rayleigh branch


        #=
            Third step: calculate the denominator
        =#

        # 3a; aerosol component
        for (aer, aer_tau) in opt.aerosol_tau
            @turbo for l in 1:N_layer
                for i1 in axes(coef_denom, 1)
                    for i2 in axes(coef_denom, 2)
                            coef_denom[i1,i2,l] += (
                                opt.aerosol_omega[aer][i,l]
                                * opt.aerosol_tau[aer][i,l]
                            )
                    end
                end
            end # End layer loop
        end

        # 3b; Rayleigh component
        if have_rayleigh
            @turbo for l in 1:N_layer
                for i1 in axes(coef_denom, 1)
                    for i2 in axes(coef_denom, 2)
                        coef_denom[i1,i2,l] += opt.rayleigh_tau[i,l]
                    end
                end
            end # End layer loop
        end # End Rayleigh branch
    end

    # Step 4, divide coef / denom
    @views opt.total_coef[:,:,:] ./= coef_denom[:,:,:]


    # Sometimes rounding errors will cause the first element to be slightly different
    # from 1.0, so we are correcting this.
    for l in 1:N_layer
        # .. this is considered "close enough"
        if abs(opt.total_coef[1,1,l] - 1) < 1e-8
            opt.total_coef[1,1,l] = 1
        end
    end


    # Sometimes rounding errors will cause the first element to be slightly different
    # from 1.0, so we are correcting this.
    for l in 1:N_layer
        # .. this is considered "close enough"
        if abs(opt.total_coef[1,1,l] - 1) < 1e-8
            opt.total_coef[1,1,l] = 1
        end
    end

end



function interpolate_aer_coef!(
    coef::AbstractArray,
    aer::AbstractAerosolType,
    ww::Number,
    ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}
)

    prop = aer.optical_property
    # Convert input wavelength to the wavelength of the aerosol optical property
    ww_use = ustrip(prop.ww_unit, ww * ww_unit)

    # Find out which pre-calculated properties are closest in wavelength
    idx1, idx2, _ = sortperm(abs.(ww_use .- prop.ww))

    # Calculation of the interpolation factor `f`
    f = (ww_use - prop.ww[idx1]) / (prop.ww[idx2] - prop.ww[idx1])

    @turbo for i1 in 1:prop.max_coefs
        for i2 in axes(coef, 2)
            coef[i1,i2] = (1 - f) * prop.coefficients[i1,i2,idx1] +
                f * prop.coefficients[i1,i2,idx2]
        end
    end

end