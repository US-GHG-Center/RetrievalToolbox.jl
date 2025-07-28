"""
$(TYPEDEF)

Type to hold optical property data that contain all necessary information to produce the
RT inputs for a spectral window `spectral_window`.

$(TYPEDFIELDS)
"""
struct EarthAtmosphereOpticalProperties{T} <: AbstractOpticalProperties where T

    # No need for this to be mutable, since all of the fields are Dicts or arrays

    "Spectral window for which to calculate optical properties"
    spectral_window::AbstractSpectralWindow
    "Dictionary to map `GasAbsorber` -> (wavelength, layer) array"
    gas_tau::Dict{GasAbsorber{T}, Array{T, 2}}
    "Dictionary to map `GasAbsorber` -> Dict of string -> (wavelength, layer) array"
    gas_derivatives::Dict{GasAbsorber{T}, Dict{String, AbstractArray}}
    "Dictionary to map `AbstractAerosol` -> (wavelength, layer) array for optical depth"
    aerosol_tau::Dict{<:AbstractAerosolType, Array{T, 2}}
    "Dictionary to map `AbstractAerosol` -> (wavelength, layer) array for single scatter albedo"
    aerosol_omega::Dict{<:AbstractAerosolType, Array{T, 2}}
    "Rayleigh optical extinction array (wavelength, layer)"
    rayleigh_tau::Array{T, 2}
    "Rayleigh optical extinction derivatives ∂τray/∂psurf (wavelength, layer)"
    rayleigh_derivatives::Array{T, 2}
    "Total extinction array (wavelength, layer)"
    total_tau::Array{T, 2}
    "Total single scatter albedo (wavelength, layer)"
    total_omega::Array{T, 2}
    "Total phase function (matrix) expansion coefficient array
     (moment, element(s) - 1 for scalar, 6 for polarized, layer)"
    total_coef::Union{Array{T, 3}, Nothing}
    "Number of dry air molecules"
    nair_dry::Vector{T}
    "Number of wet air molecules"
    nair_wet::Vector{T}
    # Temp array - shape of [spectral_window.N_hires]
    tmp_Nhi1::Vector{T}
    # Temp array - shape of [spectral_window.N_hires]
    tmp_Nhi2::Vector{T}
    # Temp array - shape of [N_layer]
    tmp_Nlay1::Vector{T}
    # Temp array - shape of [N_layer]
    tmp_Nlay2::Vector{T}
    # Temp array - same shape as total_coef
    tmp_coef::Union{Array{T, 3}, Nothing}
    # Temp array - same shape of total_coef but only keep the β coefficents (moment, 1, layer)
    tmp_coef_scalar::Union{Array{T, 3}, Nothing}

end
