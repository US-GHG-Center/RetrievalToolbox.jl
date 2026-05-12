"""
    MieMomAerosolProperty{T <: Real} <: AbstractAerosolProperty

# Fields

- `ww::Vector{T}`: Wavelength or wavenumber at which the phase function is given
- `ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}`: Unit of wavelength or
  wavenumber
- `Q_sca::Vector{T}`: Scattering efficiency per wavelength or wavenumber
- `Q_ext::Vector{T}`: Extinction efficiency per wavelength or wavenumber
- `omega::Vector{T}`: Single scatter albedo per wavelength or wavenumber
- `max_coefs::Int`: Number of phase function expansion coefficients
- `coefficients::Array{T, 3}`: Phase function expansion coefficients (coefficient, P
  element, ww)
- `work::Array{T, 3}`: Workspace for coefficient calculations (e.g. interpolation to some
  wavelength)
- `convention::NTuple{6, Symbol}`: Convention: which element in the phase matrix do
  entries belong to?

# Details

For now MOM files are expected to have all 6 phase function expansion coefficients,
regardless if further computations require the full phase matrix (polarization) or
not.

The `convention` entry relates the 6 elements of the coefficient array to the
more meaningful elements of the phase matrix expansion coefficients. Expected values
are `(:a1, :a2, :a3, :a4, :b1, :b2)` or some permutation of those.

"""
struct MieMomAerosolProperty{T<:Real} <: AbstractAerosolProperty

    ww::Vector{T}
    ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}

    Q_sca::Vector{T}
    Q_ext::Vector{T}
    omega::Vector{T}
    max_coefs::Int
    coefficients::Array{T, 3}
    work::Array{T, 3}
    convention::NTuple{6, Symbol}

end


"""
    GaussAerosol{T <: Real} <: AbstractAerosolType

Represents an aerosol that is vertically distributed along the (vertical) pressure axis
according to a Gaussian distribution around a central pressure and some width. Functions
down the line interpret the `total_optical_depth` as the full integrated optical depth
from the surface pressure up to the top of the atmosphere.

# Fields

- `aerosol_name::String`: Aerosol name
- `optical_property::AbstractAerosolProperty`: Optical properties of this aerosol
- `relative_pressure::Bool`: If this aerosol layer pressure and width is defined as a
  fraction of surface pressure
- `pressure::T`: Center of the distribution in pressure units
- `width::T`: Width in pressure units (equal to σ of a Gaussian)
- `pressure_unit::Union{Unitful.PressureUnits, Unitful.DimensionlessUnits}`:
- `ww_reference::T`: Pressure unit, or dimensionless unit if coordinates are relative to
  surface pressure
- `ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}`: Wavelength/wavenumber unit
  wavelength/wavenumber at which the total optical depth is specified
- `total_optical_depth::T`: Total-column optical depth of the distribution at valid
  pressure coordinates, and the reference wavelength

"""
mutable struct GaussAerosol{T <: Real} <: AbstractAerosolType

    aerosol_name::String
    optical_property::AbstractAerosolProperty
    relative_pressure::Bool
    pressure::T
    width::T
    pressure_unit::Union{Unitful.PressureUnits, Unitful.DimensionlessUnits}
    ww_reference::T
    ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}
    total_optical_depth::T

end

"""
$(TYPEDSIGNATURES)

Pretty printing for spectral window types
"""
function show(io::IO, ::MIME"text/plain", aer::GaussAerosol)

    println(io, "Gaussian aerosol: $(aer.aerosol_name)")

end

"""
$(TYPEDSIGNATURES)

Brief pretty printing for spectral window types
"""
function show(io::IO, aer::GaussAerosol)

    print(io, "GaussAerosol: $(aer.aerosol_name)")

end

"""
Not implemented yet.
"""
struct TriangleAerosol{T} <: AbstractAerosolType

    optical_property::AbstractAerosolProperty
    height_level::T
    total_optical_depth::T

end

"""
Not implemented yet.
"""
struct HGAerosolProperty{T} <: AbstractAerosolProperty

    g_value::T

end