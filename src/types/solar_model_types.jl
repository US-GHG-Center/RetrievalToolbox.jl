"""
An empty `AbstractSolarModel` type for use in applications that
do not require any actual implementation. For example, if one wishes
to calculate only optical properties, but not apply any radiative
transfer model.
"""
struct EmptySolarModel{T} <: AbstractSolarModel end

mutable struct OCOHDFSolarModel{T} <: AbstractSolarModel

    file_name::String
    band_number::Int

    ww::Vector{T}
    transmittance::Vector{T}
    continuum::Vector{T}

    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    irradiance_unit::Unitful.Units

end

mutable struct TSISSolarModel{T} <: AbstractSolarModel

    file_name::String

    ww::Vector{T}
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    irradiance::Vector{T}
    irradiance_unit::Unitful.Units
end

mutable struct ListSolarModel{T} <: AbstractSolarModel

    ww::Vector{T}
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    irradiance::Vector{T}
    irradiance_unit::Unitful.Units

end