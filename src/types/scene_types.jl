struct NoPhysicalScene <: NoAtmosphereScene

end

struct SpaceSunStaringScene <: NoAtmosphereScene

end

mutable struct EarthScene <: AtmosphereScene

    # The EarthAtmosphere belonging to this scene
    atmosphere::EarthAtmosphere
    # The surface dictionary which maps surfaces to their `AbstractSpectralWindow`
    surfaces::Dict{<:AbstractSpectralWindow, <:AbstractSurface}
    # The observer for which we calculate radiances
    observer::AbstractObserver
    # A location object that describes where on Earth this measurement was taken
    location::AbstractLocation
    # The solar zenith angle, to be entered in Degrees!
    solar_zenith::Number
    # The solar azimuth angle, to be entered in Degrees!
    solar_azimuth::Number
    # The time of measurement
    time::DateTime

end

"""
Pretty printing for EarthScene types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", scene::EarthScene)

    println(io, "EarthScene at $(scene.location)")

end

"""
Brief pretty printing for EarthScene types

$(SIGNATURES)
"""
function show(io::IO, scene::EarthScene)

    println(io, "EarthScene at $(scene.location)")

end