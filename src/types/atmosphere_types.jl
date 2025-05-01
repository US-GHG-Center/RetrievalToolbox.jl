"""
$(TYPEDFIELDS)
"""
struct EarthAtmosphere{T<:AbstractFloat} <: AbstractAtmosphere

    "Vector of atmosphere elements present in this atmosphere"
    atm_elements::Vector{<:AbstractAtmosphereElement}
    "Number of retrieval levels in this atmosphere"
    N_level::Int
    "Number of retrieval layers in this atmosphere"
    N_layer::Int
    "Pressure level locations"
    pressure_levels::Vector{T}
    "Mid-layer pressure locations"
    pressure_layers::Vector{T}
    "Pressure unit"
    pressure_unit::Unitful.PressureUnits

    "Number of meteorological levels in this atmosphere"
    N_met_level::Int
    "Number of meteorological layers in this atmosphere"
    N_met_layer::Int
    "Meteorological pressure level locations"
    met_pressure_levels::Vector{T}
    "Mid-layer pressure locations for meteorology"
    met_pressure_layers::Vector{T}
    "Pressure units for meteorology"
    met_pressure_unit::Unitful.PressureUnits
    "Temperatures at pressure levels"
    temperature_levels::Vector{T}
    "Temperatures at mid-layer pressures"
    temperature_layers::Vector{T}
    "Temperature unit"
    temperature_unit::Unitful.Units{U, Unitful.𝚯, nothing} where U
    "Specific humidity at pressure levels"
    specific_humidity_levels::Vector{T}
    "Specific humidity at mid-layer pressures"
    specific_humidity_layers::Vector{T}
    "Specific humidity unit"
    specific_humidity_unit::Unitful.DimensionlessUnits
    "Altitude at pressure levels"
    altitude_levels::Vector{T}
    "Altitude at mid-layer pressures"
    altitude_layers::Vector{T}
    "Altitude units"
    altitude_unit::Unitful.LengthUnits
    "Gravity at pressure levels"
    gravity_levels::Vector{T}
    "Gravity at mid-layer pressures"
    gravity_layers::Vector{T}
    "Gravity unit"
    gravity_unit::Unitful.AccelerationUnits

end

"""
Pretty printing for EarthAtmosphere types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", atm::EarthAtmosphere)

    println(io, "EarthAtmosphere (Nlay=$(atm.N_layer))")

end

"""
Brief pretty printing for EarthAtmosphere types

$(SIGNATURES)
"""
function show(io::IO, atm::EarthAtmosphere)

    println(io, "EarthAtmosphere (Nlay=$(atm.N_layer))")

end