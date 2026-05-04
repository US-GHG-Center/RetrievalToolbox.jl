struct EarthLocation{T1, T2} <: AbstractLocation

    longitude::T1
    latitude::T1
    elevation::T2
    elevation_unit::Unitful.LengthUnits

end

"""
Pretty printing for `EarthLocation` types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", el::EarthLocation)

    println(io, "EarthLocation")
    println(io, "Longitude: $(el.longitude)")
    println(io, "Latitude:  $(el.latitude)")
    println(io, "Elevation:  $(el.elevation * el.elevation_unit)")

end

"""
Brief pretty printing for `EarthLocation` types

$(SIGNATURES)
"""
function show(io::IO, el::EarthLocation)

    print(io, "EarthLocation: $(el.longitude)/$(el.latitude)/$(el.elevation * el.elevation_unit)")

end