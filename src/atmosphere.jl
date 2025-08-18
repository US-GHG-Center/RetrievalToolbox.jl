"""
$(TYPEDSIGNATURES)

Creates an empty `EarthAtmosphere` object with the specified number of retrieval grid
levels, meteorological grid levels, array type `T` as well as units for the profiles.
"""
function create_empty_EarthAtmosphere(
    Nlev::Integer,
    Nlev_met::Integer,
    T::Type{<:Real};
    pressure_unit::Unitful.PressureUnits=u"Pa",
    met_pressure_unit::Unitful.PressureUnits=u"Pa",
    temperature_unit::Unitful.Units{U, Unitful.𝚯, nothing} where U=u"K",
    specific_humidity_unit::Unitful.DimensionlessUnits=u"kg/kg",
    altitude_unit::Unitful.LengthUnits=u"km",
    gravity_unit::Unitful.AccelerationUnits=u"m/s^2"
    )

    # Create an empty vector for atmosphere elements
    atm_elements = AbstractAtmosphereElement[]

    return EarthAtmosphere(
        atm_elements,
        Nlev,
        Nlev - 1,
        zeros(T, Nlev),
        zeros(T, Nlev - 1),
        pressure_unit,
        Nlev_met,
        Nlev_met - 1,
        zeros(T, Nlev_met),
        zeros(T, Nlev_met - 1),
        met_pressure_unit,
        zeros(T, Nlev_met),
        zeros(T, Nlev_met - 1),
        temperature_unit,
        zeros(T, Nlev_met),
        zeros(T, Nlev_met - 1),
        specific_humidity_unit,
        zeros(T, Nlev_met),
        zeros(T, Nlev_met - 1),
        altitude_unit,
        zeros(T, Nlev_met),
        zeros(T, Nlev_met - 1),
        gravity_unit
    )

end

"""
$(TYPEDSIGNATURES)

In-place calculation of mid-layer values for all relevant profiles in an `EarthAtmosphere`
object (p, p MET, q, T, z, g).
"""
function calculate_layers!(atm::EarthAtmosphere)

    levels_to_layers!(atm.pressure_layers, atm.pressure_levels)
    levels_to_layers!(atm.met_pressure_layers, atm.met_pressure_levels)
    levels_to_layers!(atm.specific_humidity_layers, atm.specific_humidity_levels)
    levels_to_layers!(atm.temperature_layers, atm.temperature_levels)
    levels_to_layers!(atm.altitude_layers, atm.altitude_levels)
    levels_to_layers!(atm.gravity_layers, atm.gravity_levels)

end


"""
Creates the 'classic' University of Leicester-type pressure grid

$(SIGNATURES)

# Details

A pressure grid is generated with fixed 5 stratospheric levels. The sixth level
is halfway between tropopause pressure and level 5, and the seventh level is the
tropopause pressure itself. Finally, the remaining levels are evenly spaced
between the tropopause and the surface pressure. The total number of pressure
levels must be greater than 8 (default: 20). If the tropopause pressure is larger
than 8000 Pa (i.e. higher in the atmosphere), levels 6 and 7 are then set to
predefined values and the tropopause is then considered to be level 7.
"""
function create_UoL_pressure_grid(
    p_surf::Unitful.Pressure,
    p_tropo::Unitful.Pressure;
    N_total=20)

    @assert N_total >= 8 "Sorry, need at least 8 pressure levels for this scheme."
    @assert p_surf > p_tropo "Surface pressure ($(p_surf)) must be > " *
        "tropopause pressure ($(p_tropo))."

    # Convert to Pa here, but only if it has units,
    # this way we can either supply this function with
    # arbitrary pressure units, or directly with just
    # `p_surf` and `p_tropo` in Pa.

    p_levels = zeros(typeof(p_surf), N_total)

    p_levels[1] = 10.0u"Pa"
    p_levels[2] = 100.0u"Pa"
    p_levels[3] = 1000.0u"Pa"
    p_levels[4] = 5000.0u"Pa"
    p_levels[5] = 8000.0u"Pa"

    if p_tropo <= 8000u"Pa"
        @debug "[ATMOS] Tropopause pressure < 8000 Pa"
        p_levels[6] = 11500.0u"Pa"
        p_levels[7] = 15000.0u"Pa"
    else
        p_levels[6] = (p_tropo - 8000.0u"Pa") / 2.0 + 8000.0u"Pa"
        p_levels[7] = p_tropo
    end

    for i in 1:N_total-7
        p_levels[7 + i] = p_levels[7] + (p_surf - p_levels[7]) / (N_total - 7) * (i)
    end

    return p_levels

end

"""
Creates the ACOS-type pressure grid on 20 levels

$(SIGNATURES)

"""
function create_ACOS_pressure_grid(
    psurf::Unitful.Pressure
    )

    p = zeros(typeof(psurf), 20)

    p[1] = 0.0001 * psurf
    for i in 2:20
        p[i] = psurf * (i - 1) / 19
    end

    return p

end


"""
$(TYPEDSIGNATURES)

Calculates local gravity based on altitude levels. The optional argument `g` can be
supplied if users want a latitude=dependent surface-level gravity. Otherwise the standard
gravity is used. Note that `g` must be in compatible units of acceleration.
"""
function calculate_gravity_from_z!(atm::EarthAtmosphere; g=nothing)

    if isnothing(g)
        @debug "[ATMOS] Using standard gravity for calculations."
        g_used = g0 # Use the standard g (not latitude-corrected)
    else
        g_used = g
    end

    for l in 1:atm.N_met_level

        z = atm.altitude_levels[l] * atm.altitude_unit
        glevel = g_used * (EARTH_RADIUS / (EARTH_RADIUS + z)) |> atm.gravity_unit
        atm.gravity_levels[l] = glevel |> ustrip

    end

end

"""
$(TYPEDSIGNATURES)

Calculates local gravity given some latitude and altitude

# Details

This code was taken from MS3 / the CSU simulator suite, with heritage related to the OCO-1
retrieval algorithm written at JPL. See https://github.com/nasa/RtRetrievalFramework at
`./lib/implementation/altitude_hydrostatic.cc`.

"""
function JPL_gravity(
    latitude::Number,
    altitude::Number;
    earth_equatorial_radius=6378178.0u"m"
    )

    if !(altitude isa Unitful.Quantity)
        altitude = altitude * u"m"
    end


    # Earth mass times gravitational constant
    gm = 3.9862216*10^14 * u"m^3/s^2"
    # Earth's angular velocity
    omega = 7.292116*10^-5 * u"rad/s"
    # (a/b)**2-1 where a & b are equatorial & polar radii
    con = 0.006738
    # 2nd harmonic coefficient of Earth's gravity field
    shc = 1.6235*10^-3

    # Convert from geodetic latitude (GDLAT) to geocentric latitude (GCLAT).
    gclat = atan(tand(latitude) / (1.0 + con))
    radius = altitude + earth_equatorial_radius

    ff = (radius / earth_equatorial_radius)^2
    hh = radius * omega^2
    ge = gm / earth_equatorial_radius^2

    gravity = (ge*(1-shc*(3.0*sin(gclat)^2-1.0)/ff)/ff-hh*cos(gclat)^2)*
        (1.0+0.5*(sin(gclat)*cos(gclat)*(hh/ge+2.0*shc/ff^2))^2)

    # This returns a quantity with a physical unit!
    return gravity

end


"""
$(TYPEDSIGNATURES)

Calculates altitude and gravity for an `EarthScene`, assuming that all other
needed quantities have been inserted accordingly, namely: pressure levels,
temperatures layers, specific humidity layers and location. Note! Layer-based values are
calculated in this function call and overwrite existing values!
"""
function calculate_altitude_and_gravity!(scene::EarthScene)

    # Rebind for convenience
    atm = scene.atmosphere
    loc = scene.location

    # Need to calculate layer-based values for T and q!
    calculate_layers!(atm)

    # Calculate z and g levels
    calculate_altitude_and_gravity_levels!(
        atm.altitude_levels,
        atm.gravity_levels,
        atm.met_pressure_levels,
        atm.temperature_layers,
        atm.specific_humidity_layers,
        atm.met_pressure_levels[end],
        loc
    )

    # At this point, atm.altitude_levels is in `m`, and atm.gravity_levels in `m/s^2`, so
    # must potentially convert back to whatever units the `atm` object demands.
    atm.altitude_levels[:] .*= (1.0u"m" / atm.altitude_unit)
    atm.gravity_levels[:] .*= (1.0u"m/s^2" / atm.gravity_unit)

    # Must calculate layer-based values for gravity and altitude.
    calculate_layers!(atm)

end


"""
$(TYPEDSIGNATURES)

Calculates altitude and gravity levels for Earth-type atmospheres, in-place. For now, this
function over-writes `altitude_levels` and `gravity_levels` in units of `m` and `m/s^2`
respectively!

# Details

Given some atmospheric inputs (p, T, q) and the scene latitude and altitude, this function
calculates the altitude and gravity profiles (on levels) corresponding to the pressure
levels. These outputs should then be used to construct atmosphere objects
(EarthAtmosphere).

At this point, `p_levels` must be in [Pa], `T_layers` in [K], `SH_layers` in [1], and the
`location.altitude` in [m]. Alternatively, Unitful arrays with units can be used.
"""
function calculate_altitude_and_gravity_levels!(
    z_levels::AbstractVector,
    g_levels::AbstractVector,
    p_levels::AbstractVector,
    T_layers::AbstractVector,
    SH_layers::AbstractVector,
    p_surf::Number,
    location::EarthLocation
    )

    Rd = GAS_CONSTANT / MM_DRY_AIR

    N_layers = length(T_layers)

    z_levels[end] = location.altitude * location.altitude_unit |> u"m" |> ustrip
    g_levels[end] = JPL_gravity(location.latitude, location.altitude) |> u"m/s^2" |> ustrip

    for i in N_layers:-1:1

        if i == N_layers
            dP = p_surf - p_levels[i]
            logratio = log(p_surf / p_levels[i])
        else
            dP = p_levels[i+1] - p_levels[i]
            logratio = log(p_levels[i+1] / p_levels[i])
        end

        # This value has a unit
        this_g_layer = JPL_gravity(location.latitude, z_levels[i+1])

        # Apply temperature unit here
        Tv = T_layers[i]u"K" * (1.0 + SH_layers[i] * (1.0 - MM_H2O_TO_AIR) / MM_H2O_TO_AIR)

        # This value also has a unit
        dz = logratio * Tv * Rd / this_g_layer
        # .. so does this
        this_g_layer = JPL_gravity(location.latitude, z_levels[i+1]u"m" + 0.5 * dz)
        # .. and this
        dz = logratio * Tv * Rd / this_g_layer
        #constant = dP / (MM_DRY_AIR * this_g_layer)

        # Here we must cast back to m and m/s^2

        z_levels[i] = z_levels[i+1]u"m" + dz |> u"m" |> ustrip
        g_levels[i] = JPL_gravity(location.latitude, z_levels[i]) |> u"m/s^2" |> ustrip

    end

end

#=

    Updating and rolling back the model atmosphere via `_update!` and `_rollback!`
    ==============================================================================

    In most retrieval applications, we want to be able to retrieve/adjust certain aspects
    of, e.g., an `EarthAtmosphere` object, or, an `AbstractAtmosphereElement` one. For
    example, when using the `TemperatureOffsetSVE` in a state vector, the forward model
    must be able to take the current value of this SVE and adjust the temperature profile
    inside the `EarthAtmosphere`. At the same time, we also want to be able to roll back
    the the atmosphere to represent the state at any given iteration, but most importantly
    the state it was in when it was first evaluated (the first-guess state).

    Some SVEs alter the state in a way that makes it difficult (or impossible) to re-set
    the atmosphere at the level of the inverse solver. This is due to the fact that the
    inversion method might not have access to many of the underlying objects of the
    forward model. Often, the forward model is a function that only takes the state vector
    as its single argument to produce radiances and Jacobians.

    The suggested way of making sure that the atmosphere is updated when needed, and can
    be rolled back (e.g. non-linear schemes in which iterations can be rejected), is the
    following. Forward models *must* use the `atmosphere_element_statevector_update!`
    function, before optical properties, radiances and Jacobians are calculated. This
    function will in-place alter the atmosphere object according to the contents of the
    state vector. *After* those calculates have taken place, users then *must* call the
    roll-back functions `atmosphere_element_statevector_update!` that revert the
    atmosphere back to its first-guess state.

    Note that the two functions are designed to work together; the `.._update!` assumes
    that the atmosphere is in its first-guess state. Therefore, applying this function in
    an inappropriate manner, such as not having rolled back the atmosphere prior to
    calling it, will have unintended consequences.

    It is also important to note that this feature is truly only needed for state vector
    elements which modify the atmospheric state in a relative fashion, e.g. scaling a
    first-guess value, or adding/subtracting to or from the first-guess value.

=#


"""
$(TYPEDSIGNATURES)

Default function to update the atmosphere element `atm_element` according to the
current value of the state vector element `sve`. This does nothing!
"""
function atmosphere_element_statevector_update!(
    atm_element::AbstractAtmosphereElement,
    sve::AbstractStateVectorElement)

    # By default, do nothing - many state vector elements
    # have no impact on the model atmosphere.
    return nothing

end

"""
$(TYPEDSIGNATURES)

Function to update all atmosphere elements in the vector `atm_elements` according to the
current value of all state vector elements in the state vector `SV`.
"""
function atmosphere_element_statevector_update!(
    atm_elements::Vector{T},
    SV::AbstractStateVector
    ) where {T<:AbstractAtmosphereElement}

    for atm in atm_elements
        for sve in SV.state_vector_elements
            atmosphere_element_statevector_update!(atm, sve)
        end
    end

end


"""
$(TYPEDSIGNATURES)

Default behavior: do nothing
"""
function atmosphere_element_statevector_rollback!(
    atm_element::AbstractAtmosphereElement,
    sve::AbstractStateVectorElement
)

    return nothing

end

function atmosphere_element_statevector_rollback!(
    atm_elements::Vector{T},
    SV::AbstractStateVector
    ) where {T<:AbstractAtmosphereElement}

    for atm in atm_elements
        for sve in SV.state_vector_elements
            atmosphere_element_statevector_rollback!(atm, sve)
        end
    end

end



"""
$(TYPEDSIGNATURES)

Updates atmospheric element `atm_element` if a `GasVMRProfileSVE` is
present.
"""
function atmosphere_element_statevector_update!(
    atm_element::GasAbsorber,
    sve::GasVMRProfileSVE
)

    if atm_element === sve.gas

        # Set gas VMR to current state vector value, accounting for units
        atm_element.vmr_levels[sve.level] = ustrip(
            atm_element.vmr_unit, sve.iterations[end] * sve.unit
        )

    end

end


"""
$(TYPEDSIGNATURES)

Updates the `GasAbsorber` volume mixing ratio profile according to `sve`.
"""
function atmosphere_element_statevector_update!(
    atm_element::GasAbsorber,
    sve::GasLevelScalingFactorSVE
    )

    #=
    We check if the `atm_element` is a `GasAbsorber` element,
    then make sure the names correspond to one another,
    and finally scale the VMR according to where the
    start and end levels tell us to.
    =#

    if atm_element === sve.gas

        idx1 = sve.start_level
        idx2 = sve.end_level

        #=
            We scale the VMR levels of the corresponding gas by the scale
            factor carried by this state vector element; within the range defined, which
            is from `idx1` through `idx2`.
        =#

        scale_factor = sve.iterations[end] * sve.unit |> NoUnits
        @views atm_element.vmr_levels[idx1:idx2] .*= scale_factor

    end

end

"""
$(TYPEDSIGNATURES)

Rolls back the `GasAbsorber` volume mixing ratio profile to its first-guess state
as defined in `sve`.
"""
function atmosphere_element_statevector_rollback!(
    atm_element::GasAbsorber,
    sve::GasLevelScalingFactorSVE
)

    if atm_element === sve.gas

        idx1 = sve.start_level
        idx2 = sve.end_level

        scale_factor_current = sve.iterations[end] * sve.unit |> NoUnits
        @views atm_element.vmr_levels[idx1:idx2] ./= scale_factor_current

    end

end


function atmosphere_element_statevector_update!(
    atm_element::GaussAerosol,
    sve::AerosolOpticalDepthSVE
    )

    # Only process if the SVE points to this aerosol!
    if sve.aerosol === atm_element

        # Need to update the AOD of the GaussAerosol according to the state vector element

        # Grab the appropriate value and turn into unitless quantity
        if sve.log
            this_value = exp(get_current_value_with_unit(sve) |> NoUnits)
        else
            this_value = get_current_value_with_unit(sve) |> NoUnits
        end

        @debug "[ATMOS] Updating $(atm_element) AOD to $(this_value)"
        atm_element.total_optical_depth = this_value

    end

end


function atmosphere_element_statevector_update!(
    atm_element::GaussAerosol,
    sve::AerosolHeightSVE
    )

    # Only process if the SVE points to this aerosol!
    if sve.aerosol === atm_element

        # Need to update the AOD of the GaussAerosol according to the state vector element

        # Grab the appropriate value and turn into unitless quantity
        if sve.log
            this_value = exp(get_current_value_with_unit(sve) |> NoUnits)
        else
            this_value = get_current_value_with_unit(sve) |> NoUnits
        end

        @debug "[ATMOS] Updating $(atm_element) height to $(this_value)"
        atm_element.pressure = this_value

    end

end

function atmosphere_element_statevector_update!(
    atm_element::GaussAerosol,
    sve::AerosolWidthSVE
    )

    # Only process if the SVE points to this aerosol!
    if sve.aerosol === atm_element

        # Need to update the AOD of the GaussAerosol according to the state vector element

        # Grab the appropriate value and turn into unitless quantity
        if sve.log
            this_value = exp(get_current_value_with_unit(sve) |> NoUnits)
        else
            this_value = get_current_value_with_unit(sve) |> NoUnits
        end

        @debug "[ATMOS] Updating $(atm_element) width to $(this_value)"
        atm_element.width = this_value

    end

end

"""
Calculates the pressure weights according to O'Dell et al. in
10.5194/amt-5-99-2012. Note that this assumes that the gas concentrations
vary linearly with pressure, as implemented in the function
`calculate_gas_optical_depth_profiles`.

$(TYPEDSIGNATURES)
"""
function create_pressure_weights(
    atm::AbstractAtmosphere;
    N_sub::Integer=10)

    # Construct pressure weight array

    # "c" is calculated at the center of the retrieval pressure layers,
    # so we must grab the (potential) high-resolution MET profiles and
    # sample them at the retrieval pressure layer values.

    q_int = linear_interpolation(
        atm.met_pressure_levels,
        atm.specific_humidity_levels,
        extrapolation_bc = Line()
        )

    g_int = linear_interpolation(
        atm.met_pressure_levels,
        atm.gravity_levels,
        extrapolation_bc = Line()
        )

    # Initialize the c and c_sub arrays with the correct units
    c = zeros(atm.N_layer) / atm.gravity_unit / unit(MM_DRY_AIR)


    for i in 1:atm.N_layer

        p_vals = LinRange(atm.pressure_levels[i], atm.pressure_levels[i+1], N_sub)

        # Equation A4 in O'Dell 2012 asks for layer-averaged quantities,
        # so we calculate c = (1-q)/(g * Mdry) for some number of points
        # within the layer, and then average it out

        c_sub = 0.0 / atm.gravity_unit / unit(MM_DRY_AIR)

        # Sub-layer loop to integrate "(1-q) / (g * Mdry)" for a layer
        # using the (potentially) finely resolved `q` (specific humidity ) and `g`
        # (gravity) profiles.
        for (j, this_p) in enumerate(p_vals)

            c_sub += (1.0 - q_int(this_p)*atm.specific_humidity_unit) /
            (g_int(this_p) * atm.gravity_unit * MM_DRY_AIR)

        end

        c[i] = c_sub / N_sub

    end

    # h' is on layers
    csum = sum(c .* diff(atm.pressure_levels) * atm.pressure_unit)
    hprime = c .* diff(atm.pressure_levels) * atm.pressure_unit / csum

    # h is on levels
    h = similar(atm.pressure_levels)
    @views h[:] .= 0

    # f_i = 0.5 for all i
    # In accordance with O'Dell 2012, this means that the gas concentrations
    # within any layer is linearly varying between the surrounding levels, and
    # we grab the mid-point.

    f = one(eltype(h)) / 2
    # fs depends on surface pressure, but in our model surface pressure is equal
    # to the final pressure level, thus fs = 1.0
    fs = one(eltype(h))

    h[1] = (1 - f) * hprime[1]
    h[end-1] = f * hprime[end-1] + (1 - fs * f) * hprime[end]
    h[end] = fs * f * hprime[end]

    for i in 2:atm.N_level - 2
        h[i] = f * hprime[i-1] + (1 - f) * hprime[i]
    end

    return h

end


"""
$(TYPEDSIGNATURES)
Calculates the XGAS from an `AbstractAtmosphere` **atm**.

The column-averaged dry-air mole fraction is calculated according to O'Dell et al.,
10.5194/amt-5-99-2012. Note that the returned quantity comes with the same unit as defined
in the `GasAbsorber` object, which defines the VMR levels of that particular gas.
"""
function calculate_xgas(atm::AbstractAtmosphere; gas_name=nothing)

    # First, calculate the pressure weighting function
    # (this is the same for all gases)
    h = create_pressure_weights(atm)

    # if no gas name is supplied, return all
    if isnothing(gas_name)

        return Dict(
            x.gas_name => sum(h .* x.vmr_levels * x.vmr_unit) for x in atm.atm_elements
            if x isa GasAbsorber
        )

    # Otherwise, just find whichever gas the user wants
    else

        gas = get_gas_from_name(atm, gas_name)
        return sum(h .* gas.vmr_levels * gas.vmr_unit)

    end
end


"""
$(TYPEDSIGNATURES)

Creates a level-based profile from one defined on middle-of-the-layer. Assumes that the
mid-layer value is well-approximated by half the value of the layer-boundary values.
"""
function create_level_from_midlayer(q)

    p = zeros(eltype(q), length(q) + 1)

    for i in 1:length(p)

        if i == 1
            p[i] = q[1]
        elseif i == length(p)
            p[i] = q[end]
        else
            p[i] = 0.5 * (q[i] + q[i-1])
        end

    end

    return p

end

function atmosphere_statevector_update!(
    atm::AbstractAtmosphere,
    SV::RetrievalStateVector
    )

    # Loop through state vector elements and make an update
    # to the atmosphere, as appropriate
    for sve in SV.state_vector_elements
        atmosphere_statevector_update!(atm, sve)
    end

end

function atmosphere_statevector_rollback!(
    atm::AbstractAtmosphere,
    SV::RetrievalStateVector
    )

    # Loop through state vector elements and make a rollback
    # of the atmosphere, as appropriate
    for sve in SV.state_vector_elements
        atmosphere_statevector_rollback!(atm, sve)
    end

end

function atmosphere_statevector_update!(
    atm::AbstractAtmosphere,
    sve::AbstractStateVectorElement
)
    # Do nothing as default!
end

function atmosphere_statevector_rollback!(
    atm::AbstractAtmosphere,
    sve::AbstractStateVectorElement
)
    # Do nothing as default!
end

"""
    Updates the atmosphere `atm` for a `TemperatureOffsetSVE` state vector element.
"""
function atmosphere_statevector_update!(
    atm::AbstractAtmosphere,
    sve::TemperatureOffsetSVE
)

    #=
        Representing temperature differences consistently betweend different units, such
        as K and °C (and possibly °F), is difficult. So for now, we require the SVE
        temperature unit to be the same as the temperature profile in the atmosphere
        object.
    =#

    if sve.unit != atm.temperature_unit
        error("Unit for SVE $(sve) must be the same as atmosphere temperature profile " *
              "($(sve.unit) vs. $(atm.temperature_unit))")
    end

    #=
        The temperature offset SVE adjusts the full T profile of the atmosphere, and it
        represents an offset value with respect to the initial atmosphere state.
    =#

    # If this is the first iteration, we simply add the current value (first guess)
    ΔT = ustrip(atm.temperature_unit, get_current_value_with_unit(sve))
    @views atm.temperature_levels[:] .+= ΔT

end

function atmosphere_statevector_rollback!(
    atm::AbstractAtmosphere,
    sve::TemperatureOffsetSVE
)

    #=
        Representing temperature differences consistently betweend different units, such
        as K and °C (and possibly °F), is difficult. So for now, we require the SVE
        temperature unit to be the same as the temperature profile in the atmosphere
        object.
    =#

    if sve.unit != atm.temperature_unit
        error("[ATMOS] Unit for SVE $(sve) must be the same as atmosphere temperature profile " *
              "($(sve.unit) vs. $(atm.temperature_unit))")
    end

    #=

        At a given iteration (i), the state vector value corresponds to some temperature
        offset ΔT such that the temperature profile for layer l is the first-guess profile
        T_l plus the offset ΔT. Restoring the original temperature profile is easily
        achieved by just subtracting the current state vector value.

    =#

    ΔT = ustrip(atm.temperature_unit, get_current_value_with_unit(sve))
    @views atm.temperature_levels[:] .-= ΔT

end



"""
$(TYPEDSIGNATURES)

Updates the specific humidity profile of the atmosphere `atm`, if an H2O profile
is available as an atmospheric element. Does *not* update gravity and altitude profiles!
"""
function update_specific_humidity_from_H2O!(atm::EarthAtmosphere)

    # Check if we have H2O as a gas inside our atmosphere
    for a in atm.atm_elements
        if a isa GasAbsorber
            # Condition .. water vapor must ALWAYS be declared as "H2O", otherwise
            # this function will not be able to find the appropriate profile.
            if lowercase(a.gas_name) == "h2o"

                @debug "[ATMOS] Water vapor gas identified: $(a)"

                # Calculate SH from H2O (at RT pressure levels)
                sh_from_h2o = H2O_VMR_to_specific_humidity.(a.vmr_levels)
                # Calculate SH at the MET pressure levels
                sh_met_levels = atmospheric_profile_interpolator_linear(
                    atm.pressure_levels,
                    sh_from_h2o,
                    atm.met_pressure_levels
                )

                # Copy values over to atmosphere object
                @views atm.specific_humidity_levels[:] .= sh_met_levels
                # Calculate the mid-layer values
                @views atm.specific_humidity_layers[:] = levels_to_layers(sh_met_levels)

                #= NOTE
                    This does NOT update gravity and altitude!
                    Users must make that update manually by invoking
                    `calculate_altitude_and_gravity_levels!`
                =#

                return true

            end
        end
    end

    return false

end


"""
$(TYPEDSIGNATURES)

Given an atmosphere object, this function returns a reference to the gas object whose name
is the same as some string `name`.
"""
function get_gas_from_name(atm::EarthAtmosphere, name::String)

    glist = filter(x -> x isa GasAbsorber && x.gas_name == name, atm.atm_elements)

    if length(glist) > 1
        @warn "[ATMOS] More than one gas found with name $(name)!"
        return nothing
    end

    if length(glist) == 0
        @warn "[ATMOS] No gas found with name $(name)!"
        return nothing
    end

    return glist[1]

end


function list_example_atmospheres()

    fdir = joinpath(@__DIR__, "..", "data", "atmospheres")
    @info fdir
    fnames = filter(x -> endswith(x, ".csv"), readdir(fdir))
    anames = [replace(x, ".csv" => "") for x in fnames]

    return anames

end


"""
$(TYPEDSIGNATURES)

Creates an `EarthAtmosphere` object based on some representative atmospheres that were
extracted from NASA GMAO's MERRA2 reanalysis and Sourish Basu's CO2/CH4. Must provide the
name of the example atmosphere `name` as well as the intended number of pressure levels
for the retrieval grid, `Nlev`, to be filled out be the user later.
"""
function create_example_atmosphere(
    name::String,
    Nlev::Integer;
    T::Type{<:Real}=Float64,
    surface_pressure::Union{Float64, Nothing}=nothing,
    altitude::Union{Float64, Nothing}=nothing
    )

    if !isnothing(surface_pressure) & !isnothing(altitude)
        throw(ArgumentError(
            "[ATMOS] Must not provide BOTH surface pressure AND altitude. Choose one."
            ))
    end

    fname = joinpath(@__DIR__, "..", "data", "atmospheres", "$(name).csv")

    if !isfile(fname)
        throw(ErrorException(
            "[ATMOS] Example atmosphere $(name) not found!"
        ))
    end

    # Read the content of the example atmosphere file
    @debug "[ATMOS] Opening atmosphere file at $(fname)"

    # First read the header to obtain info on units etc.
    units_list = nothing

    open(fname, "r") do file
        while !eof(file)
            line = readline(file)
            line[1] != '#' && break  # Get out if line does not start with #
            if startswith(line, "# units:")
                # Grab the units
                units_list = uparse.(split(split(line, "units:")[end], ","))
            end
        end
    end

    # These are now our profiles, defined in layers
    csv = CSV.File(fname, comment="#")
    # Note that the units are in the same order as the csv columns, so we can make
    # a helpful dict here that keeps the units easily accessible.
    csv_units = Dict(csv.names[i] => units_list[i] for i in 1:csv.cols)

    Nlev_met = length(csv)

    # Create the empty atmosphere
    atm = create_empty_EarthAtmosphere(
        Nlev, # Number of retrieval pressure levels
        Nlev_met, # Number of met pressure levels
        T; # Number type
        pressure_unit=csv_units[:pressure], # Retrieval grid pressure unit
        met_pressure_unit=csv_units[:pressure], # Met grid pressure unit
        temperature_unit=csv_units[:temperature], # Temperature unit
        specific_humidity_unit=csv_units[:specific_humidity], # Specific humidity unit
        altitude_unit=csv_units[:altitude], # Altitude unit
        gravity_unit=u"m/s^2" # Gravity unit
    )

    ingest!(atm, :met_pressure_levels,
        csv.pressure * csv_units[:pressure])
    ingest!(atm, :temperature_levels,
        csv.temperature * csv_units[:temperature])
    ingest!(atm, :specific_humidity_levels,
        csv.specific_humidity * csv_units[:specific_humidity])
    ingest!(atm, :altitude_levels,
        csv.altitude * csv_units[:altitude])

    # For this example atmosphere, use a simple approach to calculate gravity based
    # on altitude (latitude is unknown yet, and probably not important)
    calculate_gravity_from_z!(atm)

    # Calculate layer quantities
    calculate_layers!(atm)

    return atm
end