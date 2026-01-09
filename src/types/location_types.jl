struct EarthLocation{T1, T2} <: AbstractLocation

    longitude::T1
    latitude::T1
    altitude::T2
    altitude_unit::Unitful.LengthUnits

end

"""
Pretty printing for `EarthLocation` types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", el::EarthLocation)

    println(io, "EarthLocation")
    println(io, "Longitude: $(el.longitude)")
    println(io, "Latitude:  $(el.latitude)")
    println(io, "Altitude:  $(el.altitude * el.altitude_unit)")

end

"""
Brief pretty printing for `EarthLocation` types

$(SIGNATURES)
"""
function show(io::IO, el::EarthLocation)

    print(io, "EarthLocation: $(el.longitude)/$(el.latitude)/$(el.altitude)")

end