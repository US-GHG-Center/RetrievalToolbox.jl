"""
An empty `AbstractSolarModel` type for use in applications that do not require any actual
implementation. For example, if one wishes to calculate only optical properties, but not
apply any radiative transfer model. Another example would be RT in thermal spectral
ranges, where solar irradiance can be omitted.
"""
struct NoSolarModel <: AbstractSolarModel end

"""
$(TYPEDFIELDS)

A type to describe a data-driven solar model that corresponds to the HDF5 file used by the
NASA RtRetrievalFramework code. See `l2_solar_model.h5` found in
https://github.com/nasa/RtRetrievalFramework/tree/master/input/common/input.
"""
mutable struct OCOHDFSolarModel{T} <: AbstractSolarModel

    file_name::String
    band_number::Int

    ww::Vector{T}
    transmittance::Vector{T}
    continuum::Vector{T}

    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    irradiance_unit::Unitful.Units

end

"""
$(TYPEDFIELDS)

A type to describe the data-driven solar model that corresponds to the LASP TSIS solar
irradiance data that can be downloaded at
https://lasp.colorado.edu/lisird/data/tsis1_hsrs_p1nm
"""
mutable struct TSISSolarModel{T} <: AbstractSolarModel

    file_name::String

    ww::Vector{T}
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    irradiance::Vector{T}
    irradiance_unit::Unitful.Units
end

"""
$(TYPEDFIELDS)

A type to hold a general data-driven solar model. The irradiance data must be entered by
users into the `.irradiance` field manually. This can be used if some line list is
available.
"""
mutable struct ListSolarModel{T} <: AbstractSolarModel

    ww::Vector{T}
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    irradiance::Vector{T}
    irradiance_unit::Unitful.Units

end

"""
This solar model type should be used when users want to work exclusively with
sun-normalized radiance and do not need any spectrally varying continuum.
"""
struct UnitSolarModel <: AbstractSolarModel

    irradiance_unit::Unitful.Units

    function UnitSolarModel()
        return new(Unitful.NoUnits)
    end

end