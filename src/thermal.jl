#=
    Routines to calculate radiances due to thermal emission
    =======================================================

    The functions in this file compute Planck radiance given some temperature and a
    spectral coordinate, wavelength or wavenumber. Since radiance may be in units of power
    per wavenumber|wavelength or photons/s per wavenumber|wavelength, there is a total of
    four different ways to calculate thermal emission.

    To make this calculation fast enough, i.e. to not encase those four possibilities into
    if/elseif statements, we write a dedicated function for each.

    At the top level, there are only two functions `Planck_radiance` that take either
    length and temperature, or wavelength and temperature, along with the requested
    radiance unit.

    The top-level functions then dispatch further down to more specialized functions that
    calculate the Planck radiance for some combination of

        [length or wavenumber], [power or photons]

    The same strategy also applies for the partial derivatives w.r.t. temperature.
    Unfortunately, Unitful does not play nice with auto diff libraries, so we have to
    write out the derivatives manually into separate functions.

    At the moment, a single call to `Planck_radiance` or `dPlanck_radiance_dT` Takes
    ~200ns, which is still significant. For a spectral window with 30,000 points and
    20 levels, this adds up to 0.12 seconds.

    TODO: make smarter/faster vectorized functions to process multiple spectral points at
    once!

=#

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

    return @inline _Planck_radiance(λ, T, dimension(rad_unit)) |> rad_unit

end

function dPlanck_radiance_dT(
    λ::Unitful.Length,
    T::Unitful.Temperature,
    rad_unit::Unitful.Units
    )

    # Note regarding the units: the unit going into `_dPlanck_radiance_dT` is the radiance
    # unit **without** the per-temperature unit. We back-attach that missing unit when we
    # produce the output.
    # This is a slightly unfortunate hack needed keep the unit dispatch mechanism for
    # Planck radiances here working..
    return @inline _dPlanck_radiance_dT(λ, T, dimension(rad_unit)) |> rad_unit / unit(T)

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

    return @inline _Planck_radiance(ν, T, dimension(rad_unit)) |> rad_unit

end



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
wavelength `λ`. FOR INTERNAL USE MOSTLY.
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
    return _Planck_radiance(λ, T, rad_dim) * x / T * exp(x) / (exp(x) - 1) |>
        u"W/m^2/sr/µm" / unit(T)

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
wavenumber `ν`. FOR INTERNAL USE MOSTLY.
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
        u"W/m^2/sr/cm^-1" / unit(T)

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

First derivative of `_Planck_radiance` with respect to temperature `T`, at some
wavelength `λ`. FOR INTERNAL USE MOSTLY.
"""
function _dPlanck_radiance_dT(
    λ::Unitful.Length,
    T::Unitful.Temperature,
    rad_dim::TYPE_PHOTON_PER_LENGTH
    )

    # First calculate in Watts
    p = _dPlanck_radiance_dT(λ, T, DIM_POWER_PER_LENGTH)
    # Then convert to ph/s. Note that we do a silly trick here since `W_to_ph` does not
    # accept W/m^2/sr/µm/K
    return W_to_ph(p * u"K", λ) / unit(T)

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


