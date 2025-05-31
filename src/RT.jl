"""
$(TYPEDSIGNATURES)

Returns the number of Stokes components used in this RT method object.
"""
function get_nstokes(rt::AbstractRTMethod)

    return size(rt.hires_radiance, 2)

end


function clear!(rt::MonochromaticRTMethod)

    # Radiances can be accessed like 2D arrays
    rt.hires_radiance[:,:] .= 0

    # Weighting functions are vectors of radiances
    for v in rt.hires_wfunctions
        v[:,:] .= 0
    end

    # Jacobians are dictionaries StateVectorElement -> Radiance
    for (k,v) in rt.hires_jacobians
        v[:,:] .= 0
    end

end

"""
$(TYPEDSIGNATURES)

Checks every entry of `rt.hires_radiance` and `rt.hires_jacobians` for non-finite values
and reports if there are any, and where. Finally, `true` is returned if none are found,
and `false` is returned if at least one entry has a non-finite value.
"""
function check_non_finites(rt::MonochromaticRTMethod)

    N_stokes = get_nstokes(rt)
    all_finite = true

    for s in 1:N_stokes # Loop through Stokes components

        radiance_bad_pos = Int[]
        jacobian_bad_pos = Dict{AbstractStateVectorElement, Vector{Int}}()

        for i in axes(rt.hires_radiance, 1) # Loop through spectral index
            if !isfinite(rt.hires_radiance[i,s])
                push!(radiance_bad_pos, i)
                all_finite = false
            end
        end

        # Loop through Jacobian attached to state vector element
        for (sve, jac) in rt.hires_jacobians
            for i in axes(jac, 1) # Loop through spectral index
                if !isfinite(jac[i,s])
                    if !haskey(jacobian_bad_pos, sve)
                        jacobian_bad_pos[sve] = Int[]
                    end
                    push!(jacobian_bad_pos[sve], i)
                    all_finite = false
                end
            end
        end

        #=
        if !(all_finite)
            @warn "Non-finite(s) in hi-res radiance for Stokes component $(s): " *
                "$(radiance_bad_pos)"
            @warn "Non-finite(s) in hi-res Jacobian $(sve) for Stokes component $(s): " *
                "$(jacobian_bad_pos)"
        end
        =#
    end


    return all_finite

end


"""
$(TYPEDSIGNATURES)

Update the solar scaler variable in the RT object `rt`, given some state vector `sv`.
"""
function solar_scaler_statevector_update!(
    rt::AbstractRTMethod,
    sv::AbstractStateVector
    )

    # First - establish if we have any SolarScalerSVEs
    if !any_SVE_is_type(sv, SolarScalerPolynomialSVE)
        return
    end


    idx_half = Int(floor(rt.optical_properties.spectral_window.N_hires // 2))
    # Bind to temp array
    delta_wl = rt.optical_properties.tmp_Nhi1
    @views delta_wl[:] .= 0.0

    @views @. delta_wl[:] = rt.optical_properties.spectral_window.ww_grid[:] -
            rt.optical_properties.spectral_window.ww_grid[idx_half]

    @views rt.solar_scaler[:] .= 0.0

    for (idx, sve) in StateVectorIterator(sv, SolarScalerPolynomialSVE)

        c = sve.iterations[end]
        @turbo for i in eachindex(rt.solar_scaler)
            rt.solar_scaler[i] += c * delta_wl[i] ^ sve.coefficient_order
        end

    end

end


#=
    Generic functions that further dispatch to specific RT and observers
=#


function calculate_dI_dVMR(
    rt::AbstractRTMethod,
    gas::GasAbsorber
)
    @warn "No specific implementation for $(rt) available!"
    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function calculate_dI_dTau(
    rt::AbstractRTMethod
)

    # Call specific implementation
    return calculate_dI_dTau(rt, rt.scene.observer)

end

# Insert functions for the Beer-Lambert RT type
include("RT_BeerLambert.jl")
# Insert functions for XRTM RT module
include("RT_XRTM.jl")
# Insert functions for the Monochromatic RT type
include("RT_Monochromatic.jl")
# Low-Streams Interpolation (RT acceleration method)
include("RT_LSI.jl")

"""
Creates a deepcopy-like copy of `rt`, but only creates new arrays for e.g.
radiances, optical depths, etc.; with the notable exception that these RT
objects have only one spectral point in them, as to represent binned optical
properties (or optical states).

# Details

TODO: add explanation and motivation
"""
function create_rt_copy_for_bins(
    rt::MonochromaticRTMethod,
    spectral_idx::Integer;
    options::Union{Nothing, U, Vector{U}}=nothing
    ) where {U <: AbstractDict}

    # Create some handy shortcuts and useful values
    opt = rt.optical_properties
    old_swin = opt.spectral_window
    N_layer = rt.scene.atmosphere.N_layer
    T = eltype(opt.total_tau)

    #=
        First, create a new spectral window object with exactly
        one spectral point in it (taken from the band center, for example).

        Note, this is a different `AbstractSpectralWindow` type, since it is not
        really a retrieval window in the same sense. We only use this for binned
        RT calculations, and they do not have the same meaning apart from being a
        required object for e.g. `EarthAtmosphereOpticalProperties`.
    =#

    new_ww = old_swin.ww_grid[spectral_idx]

    new_swin = BinnedSpectralWindow(
        old_swin.window_name,
        [new_ww],
        old_swin.ww_unit,
        old_swin.ww_reference,
        1,
        old_swin, # Store reference to the original spectral window
        spectral_idx
    )

    #=
        Second step - produce the internal contents of the optical properties
    =#

    gas_tau = Dict(gas => zeros(1, N_layer) for gas in keys(opt.gas_tau))

    gas_derivatives = Dict{GasAbsorber{T}, Dict{String, AbstractArray}}()
    for gas in keys(opt.gas_derivatives)
        d = Dict{String, AbstractArray}()
        for (s,v) in opt.gas_derivatives[gas]

            # It's a little more complicated here, since we might have
            # various array shapes for the different partial derivatives.
            # Let us just infer the shape and replace the N_hires with 1.

            # This just copies over the old array size, but replaces a
            # potential wavelength axis with 1.
            new_size = [x == old_swin.N_hires ? 1 : x for x in size(v)]
            d[s] = zeros(T, new_size...)
        end
        gas_derivatives[gas] = d
    end

    aerosol_tau = Dict{AbstractAerosolType, Matrix{T}}()
    aerosol_omega = Dict{AbstractAerosolType, Matrix{T}}()
    for aer in keys(opt.aerosol_tau)
        aerosol_tau[aer] = zeros(T, 1, N_layer)
        aerosol_omega[aer] = zeros(T, 1, N_layer)
    end

    rayleigh_tau = zeros(T, 1, N_layer)
    rayleigh_derivatives = zeros(T, 1, N_layer)
    total_tau = zeros(T, 1, N_layer)
    total_omega = zeros(T, 1, N_layer)

    # Total coefficients can be a `Nothing` or an array
    if isnothing(opt.total_coef)
        total_coef = nothing
        tmp_coef = nothing
        tmp_coef_scalar = nothing
    else
        total_coef = copy(opt.total_coef)
        tmp_coef = copy(opt.total_coef)
        tmp_coef_scalar = copy(opt.tmp_coef_scalar)
    end

    nair_dry = copy(opt.nair_dry)
    nair_wet = copy(opt.nair_wet)

    tmp_Nhi1 = zeros(T, 1)
    tmp_Nhi2 = zeros(T, 1)
    tmp_Nlay1 = zeros(T, N_layer)
    tmp_Nlay2 = zeros(T, N_layer)

    # This generates a new optical property object (no copy or reference)
    new_opt = EarthAtmosphereOpticalProperties(
        new_swin,
        gas_tau,
        gas_derivatives,
        aerosol_tau,
        aerosol_omega,
        rayleigh_tau,
        rayleigh_derivatives,
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


    #=
        Third step, we create the remaining inputs for the new RT object
        TODO (for the time being, this works only with `MonochromaticRTMethod`)
    =#

    new_model = rt.model
    if isnothing(options)
        @debug "No RT model options supplied, inheriting from `rt`"
        new_model_options = rt.model_options
    else
        @debug "User-supplied model options for new RT object!"
        new_model_options = options
    end

    # The following quantities need to be referenced only, there is no need
    # to cast a copy of those items.
    new_scene = rt.scene
    new_solar_model = rt.solar_model
    new_state_vector = rt.state_vector
    new_wfunctions_map = rt.wfunctions_map

    # These, however, need to be created anew, since they must match up with the
    # number of spectral points in the spectral window object, which now only has
    # one single wavelength/wavenumber in it..

    # This provides us with, e.g. `ScalarRadiance` if `rt.hires_radiance`
    # is a `ScalarRadaiance{Float64}`.
    RadType = typeof(rt.hires_radiance).name.wrapper

    new_hires_radiance = RadType(T, 1)

    if isnothing(rt.hires_jacobians)
        new_hires_jacobians = nothing
    else
        new_hires_jacobians = Dict(
            sve => RadType(T, 1)
            for sve in rt.state_vector.state_vector_elements
        )
    end

    if isnothing(rt.hires_wfunctions)
        new_hires_wfunctions = nothing
    else
        N_wf = size(rt.hires_wfunctions, 1)
        new_hires_wfunctions = [
            RadType(T, 1) for i in 1:N_wf ]
    end

    new_solar_scaler = ones(T, 1)
    new_hires_solar = RadType(T, 1)
    new_hires_solar.I[1] = rt.hires_solar.I[spectral_idx]

    new_rt = MonochromaticRTMethod(
        new_model,
        new_model_options,
        new_scene,
        new_opt,
        new_solar_model,
        new_state_vector,
        new_hires_solar,
        new_hires_radiance,
        new_hires_jacobians,
        new_hires_wfunctions,
        new_wfunctions_map,
        rt.radiance_unit,
        new_solar_scaler
    )

    return new_rt
end