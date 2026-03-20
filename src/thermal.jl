"""
$(TYPEDSIGNATURES)

Calculates the isotropic radiance emitted given some temperature `T` at wavelength `λ`.
The result is forced into units of W m⁻² sr⁻¹ µm⁻¹. Users must make sure they then
convert the result into the units of radiance they need. FOR INTERNAL USE MOSTLY.
"""
function _Planck_radiance(
    λ::Unitful.Length,
    T::Unitful.Temperature,
    rad_dim::TYPE_POWER_PER_LENGTH
    )

    # Take physical constants from constants.jl file, rebind to short vars
    kB = BOLTZMANN
    h = PLANCK
    c = SPEED_OF_LIGHT

    return (2 * h * c^2) / (λ^5) / (exp((h * c) / (λ * kB * T)) - 1) |> u"W/m^2/sr/µm"

end

"""
$(TYPEDSIGNATURES)

First derivative of `_Planck_radiance` with respect to temperature `T`, at some
wavelength `λ`.
"""
function _dPlanck_radiance_dT(
    λ::Unitful.Length,
    T::Unitful.Temperature,
    rad_dim::TYPE_POWER_PER_LENGTH
    )

    # Take physical constants from constants.jl file, rebind to short vars
    kB = BOLTZMANN
    h = PLANCK
    c = SPEED_OF_LIGHT

    x = (h * c) / (λ * kB * T)
    return _Planck_radiance(λ, T, rad_dim) * x / T * exp(x) / (exp(x) - 1) |> u"W/m^2/sr/µm/K"

end

"""
$(TYPEDSIGNATURES)

Calculates the isotropic radiance emitted given some temperature `T` at wavenumber `ν`.
The result is forced into units of W m⁻² sr⁻¹ (cm⁻¹)⁻¹. Users must make sure they then
convert the result into the units of radiance they need. FOR INTERNAL USE MOSTLY.
"""
function _Planck_radiance(
    ν::Unitful.Wavenumber,
    T::Unitful.Temperature,
    rad_dim::TYPE_POWER_PER_WAVENUMBER
    )

    # Take physical constants from constants.jl file, rebind to short vars
    kB = BOLTZMANN
    h = PLANCK
    c = SPEED_OF_LIGHT

    return (2 * h * c^2 * ν^3) / (exp((h * c * ν) / (kB * T)) - 1) |>
        u"W/m^2/sr/cm^-1"

end

"""
$(TYPEDSIGNATURES)

First derivative of `_Planck_radiance` with respect to temperature `T`, at some
wavenumber `ν`.
"""
function _dPlanck_radiance_dT(
    ν::Unitful.Wavenumber,
    T::Unitful.Temperature,
    rad_dim::TYPE_POWER_PER_WAVENUMBER
    )

    # Take physical constants from constants.jl file, rebind to short vars
    kB = BOLTZMANN
    h = PLANCK
    c = SPEED_OF_LIGHT

    x = (h * c * ν) / (kB * T)
    return _Planck_radiance(ν, T, rad_dim) * x / T * exp(x) / (exp(x) - 1) |>
        u"W/m^2/sr/cm^-1/K"

end


"""
$(TYPEDSIGNATURES)

Calculates the isotropic radiance emitted given some temperature `T` at wavelength `λ`.
The result is forced into units of ph s⁻¹ m⁻² sr⁻¹ µm⁻¹. Users must make sure they then
convert the result into the units of radiance they need. FOR INTERNAL USE MOSTLY.
"""
function _Planck_radiance(
    λ::Unitful.Length,
    T::Unitful.Temperature,
    rad_dim::TYPE_PHOTON_PER_LENGTH
    )

    # First calculate in Watts
    p = _Planck_radiance(λ, T, DIM_POWER_PER_LENGTH)
    # Then convert to ph/s
    return W_to_ph(p, λ)

end

"""
$(TYPEDSIGNATURES)

Calculates the isotropic radiance emitted given some temperature `T` at wavenumber `ν`.
The result is forced into units of ph s⁻¹ m⁻² sr⁻¹ (cm⁻¹)⁻¹. Users must make sure they
then convert the result into the units of radiance they need. FOR INTERNAL USE MOSTLY.
"""
function _Planck_radiance(
    ν::Unitful.Wavenumber,
    T::Unitful.Temperature,
    rad_dim::TYPE_PHOTON_PER_WAVENUMBER
    )

    # First calculate in Watts
    p = _Planck_radiance(ν, T, DIM_POWER_PER_WAVENUMBER)
    # Then convert to ph/s
    return W_to_ph(p, ν)

end

"""
    Planck_radiance(
        λ::Unitful.Length,
        T::Unitful.Temperature,
        rad_unit::Unitful.Units
    ) -> Unitful quantity converted to `rad_unit`

Top-level function to calculate thermal (Planckian blackbody) radiance emission for some
wavelength `λ` and temperature `T`. Input parameters must be Unitful quantities with
associated and correct units.

Only radiance units compatible with either $(rad_W_wl) or $(rad_ph_wl) will work, since
the lower-lying function will dispatch to another function that must be compatible with
those units.
"""
function Planck_radiance(
    λ::Unitful.Length,
    T::Unitful.Temperature,
    rad_unit::Unitful.Units
    )

    return _Planck_radiance(λ, T, dimension(rad_unit)) |> rad_unit

end

"""
    Planck_radiance(
        λ::Unitful.Length,
        T::Unitful.Temperature,
        rad_unit::Unitful.Units
    ) -> Unitful quantity converted to `rad_unit`

Top-level function to calculate thermal (Planckian blackbody) radiance emission for some
wavenumber `ν` and temperature `T`. Input parameters must be Unitful quantities with
associated and correct units.

Only radiance units compatible with either $(rad_W_wn) or $(rad_ph_wn) will work, since
the lower-lying function will dispatch to another function that must be compatible with
those units.
"""
function Planck_radiance(
    ν::Unitful.Wavenumber,
    T::Unitful.Temperature,
    rad_unit::Unitful.Units
    )

    return _Planck_radiance(ν, T, dimension(rad_unit)) |> rad_unit

end
