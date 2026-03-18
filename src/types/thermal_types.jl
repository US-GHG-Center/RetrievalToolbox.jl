"""
Defines an element that represents isotropic thermal emission from the atmosphere itself.
"""
struct ThermalAtmosphereIsotropic <: AbstractThermalAtmosphere

    # Nothing needed yet

end


"""
Pretty printing for ThermalAtmosphereIsotropic types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", th::ThermalAtmosphereIsotropic)

    println(io, "ThermalAtmosphereIsotropic")

end

"""
Brief pretty printing for gas absorber types

$(SIGNATURES)
"""
function show(io::IO, th::ThermalAtmosphereIsotropic)

    print(io, "ThermalAtmosphereIsotropic")

end


"""
$(TYPEDFIELDS)

Defines an element that represents a surface with a finite temperature, defined only by
the surface temperature itself. The radiance emission from the surface is isotropoic. Only
Kelvin are supported for the temperature unit, using any other temperature unit will throw
an error!

For now, this surface temperature is applied scene-wide, so each surface that is coupled
to spectral windows, will be considered to have this surface temperature.
"""
struct ThermalSurfaceIsotropic{T} <: AbstractThermalSurface

    temperature :: T
    temperature_unit :: Unitful.TemperatureUnits

    function ThermalSurfaceIsotropic(temp::T, unit::Unitful.TemperatureUnits) where {T}
        if unit !== u"K"
            throw(ArgumentError("temperature_unit must be Kelvin (u\"K\"). \
                Received $unit instead."))
        end
        new{T}(temp, unit)
    end

end

"""
Pretty printing for ThermalSurfaceIsotropic types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", th::ThermalSurfaceIsotropic)

    println(io, "ThermalSurfaceIsotropic ($(th.temperature * th.temperature_unit))")

end

"""
Brief pretty printing for gas absorber types

$(SIGNATURES)
"""
function show(io::IO, th::ThermalSurfaceIsotropic)

    print(io, "ThermalSurfaceIsotropic ($(th.temperature * th.temperature_unit))")

end