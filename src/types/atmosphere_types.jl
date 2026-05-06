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

    "Number of meteorological model levels in this atmosphere"
    N_met_level::Int
    "Meteorological vertical pressure model coordinate"
    met_pressure::Vector{T}
    "Pressure units for meteorology"
    met_pressure_unit::Unitful.PressureUnits
    "Temperatures at meteorology model pressure levels"
    temperature::Vector{T}
    "Temperature unit"
    temperature_unit::Unitful.Units{U, Unitful.𝚯, nothing} where U
    "Specific humidity at meteorology model pressure levels"
    specific_humidity::Vector{T}
    "Specific humidity unit"
    specific_humidity_unit::Unitful.DimensionlessUnits
    "Altitude at meteorology model pressure levels"
    altitude::Vector{T}
    "Altitude units"
    altitude_unit::Unitful.LengthUnits
    "Gravity at meteorology model pressure levels"
    gravity::Vector{T}
    "Gravity unit"
    gravity_unit::Unitful.AccelerationUnits

end

"""
Pretty printing for EarthAtmosphere types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", atm::EarthAtmosphere)

    println(io, "EarthAtmosphere (Nlev=$(atm.N_level))")

end

"""
Brief pretty printing for EarthAtmosphere types

$(SIGNATURES)
"""
function show(io::IO, atm::EarthAtmosphere)

    println(io, "EarthAtmosphere (Nlev=$(atm.N_level))")

end

#=
    We changed the strucutre of EarthAtmosphere. Let's give users a change to change their
    codes before we remove this deprecation warning.
    There is no need to overload `setproperty!`, since EarthAtmosphere is not a mutable
    sturct, and thus users cannot change the fields themselves, only the values of the
    vectors/arrays.
=#

function Base.getproperty(obj::EarthAtmosphere, sym::Symbol)

    if sym === :met_pressure_levels
        @warn "Field `met_pressure_levels` no longer exists. Use `met_pressure`!"
        return getfield(obj, :met_pressure)
    elseif sym === :temperature_levels
        @warn "Field `temperature_levels` no longer exists. Use `temperature`!"
        return getfield(obj, :temperature)
    elseif sym === :temperature_layers
        @error "Field `temperature_layers` no longer exists!"
        throw(FieldError(EarthAtmosphere, sym))
    elseif sym === :specific_humidity_levels
        @warn "Field `specific_humidity_levels` no longer exists. Use `specific_humidity`!"
        return getfield(obj, :specific_humidity)
    elseif sym === :specific_humidity_layers
        @error "Field `specific_humidity_layers` no longer exists!"
        throw(FieldError(EarthAtmosphere, sym))
    elseif sym === :altitude_levels
        @warn "Field `altitude_levels` no longer exists. Use `altitude`!"
        return getfield(obj, :altitude)
    elseif sym === :altitude_layers
        @error "Field `altitude_layers` no longer exists!"
        throw(FieldError(EarthAtmosphere, sym))
    elseif sym === :gravity_levels
        @warn "Field `gravity_levels` no longer exists. Use `gravity`!"
        return getfield(obj, :gravity)
    elseif sym === :gravity_layers
        @error "Field `gravity_layers` no longer exists!"
        throw(FieldError(EarthAtmosphere, sym))
    else
        return getfield(obj, sym)
    end

end