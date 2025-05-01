"""

    $(TYPEDFIELDS)
    # Details

    For now MOM files are expected to have all 6 phase function expansion coefficients,
    regardless if further computations require the full phase matrix (polarization) or
    not.

    The `convention` entry relates the 6 elements of the coefficient array to the
    more meaningful elements of the phase matrix expansion coefficients. Expected values
    are `(:a1, :a2, :a3, :a4, :b1, :b2)` or some permutation of those.

"""
struct MieMomAerosolProperty{T} <: AbstractAerosolProperty
    
    "Wavelength or wavenumber at which the phase function is given"
    ww::Vector{T}
    "Unit of wavelength or wavenumber"
    ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}
    
    "Scattering efficiency per wavelength or wavenumber"
    Q_sca::Vector{T}
    "Extinction efficiency per wavelength or wavenumber"
    Q_ext::Vector{T}
    "Single scatter albedo per wavelength or wavenumber"
    omega::Vector{T}
    "Number of phase function expansion coefficients"
    max_coefs::Int
    "Phase function expansion coefficients (coefficient, P element, ww)"
    coefficients::Array{T, 3}
    "Workspace for coefficient calculations (e.g. interpolation to some wavelength)"
    work::Array{T, 3}
    "Convention: which element in the phase matrix do entries belong to?"
    convention::NTuple{6, Symbol}

end

struct HGAerosolProperty{T} <: AbstractAerosolProperty

    g_value::T

end


mutable struct GaussAerosol{T} <: AbstractAerosolType

    "Aerosol name"
    aerosol_name::String
    "Optical properties of this aerosol"
    optical_property::AbstractAerosolProperty
    "If this aerosol layer pressure and width is defined as a fraction of surface pressure"
    relative_pressure::Bool
    "Center of the distribution in pressure units"
    pressure::T
    "Width in pressure units (equal to σ of a Gaussian)"
    width::T
    "Pressure unit, or dimensionless unit if coordinates are relative to psurf"
    pressure_unit::Union{Unitful.PressureUnits, Unitful.DimensionlessUnits}
    # Reference wavelength/wavenumber at which the total optical depth is specified
    ww_reference::T
    # Wavelength/wavenumber unit
    ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}
    "Total-column optical depth of the distribution at valid pressure coordinates, 
     and the reference wavelength"
    total_optical_depth::T

end

"""
Pretty printing for spectral window types

$(TYPEDSIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", aer::GaussAerosol)

    println(io, "Gaussian aerosol: $(aer.aerosol_name)")

end

"""
Brief pretty printing for spectral window types

$(TYPEDSIGNATURES)
"""
function show(io::IO, aer::GaussAerosol)

    print(io, "GaussAerosol: $(aer.aerosol_name)")

end


struct TriangleAerosol{T} <: AbstractAerosolType

    optical_property::AbstractAerosolProperty
    height_level::T
    total_optical_depth::T

end
