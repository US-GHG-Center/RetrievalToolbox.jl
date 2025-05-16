struct NoRTMethod <: AbstractRTMethod

    scene::NoAtmosphereScene
    state_vector::AbstractStateVector

end

struct MonochromaticRTMethod{
    U <: AbstractOpticalProperties,
    V <: AbstractSolarModel,
    R1 <: Radiance,
    S <: AbstractStateVector,
    SS <: AbstractStateVectorElement,
    } <: AbstractRTMethod

    "Which RT model to use, at the moment, only :XRTM is supported"
    model::Symbol
    "Option dictionary to set parameters for the underlying RT code"
    model_options::Union{Vector{T}, T} where {T <: AbstractDict}
    "`AtmosphereScene` for which the RT is to be computed"
    scene::AtmosphereScene
    "`AbstractOpticalProperties` that provide the basic inputs for RT"
    optical_properties::U
    "The solar model"
    solar_model::V
    "The state vector"
    state_vector::S
    "Buffer to hold the high-resolution solar irradiance"
    hires_solar::R1
    "Buffer to hold the high-resolution at-instrument radiance"
    hires_radiance::R1
    "Buffer to hold the high-resolution radiance Jacobians"
    hires_jacobians::Union{
        Nothing,
        Dict{SS, R1}
        }
    "Buffer to hold the RT-computed weighting functions"
    hires_wfunctions::Union{
        Nothing,
        Vector{R1}
    }
    "Mapper dictionary to assign RT weighting functions to required Jacobians"
    wfunctions_map::Dict{Any, Vector{Int}}
    "Radiance unit"
    radiance_unit::Unitful.Units
    "Solar scale array, to be multiplied into the solar irradiance"
    solar_scaler::AbstractVector

end

struct BeerLambertRTMethod <: AbstractRTMethod

    "`AtmosphereScene` for which the RT is to be computed"
    scene::AtmosphereScene
    "`AbstractOpticalProperties` that provide the basic inputs for RT"
    optical_properties::AbstractOpticalProperties
    "The solar model"
    solar_model::AbstractSolarModel
    "The state vector"
    state_vector::AbstractStateVector
    "Buffer to hold the high-resolution solar irradiance"
    hires_solar::Radiance
    "Buffer to hold the high-resolution at-instrument radiance"
    hires_radiance::Radiance
    "Buffer to hold the high-resolution radiance Jacobians"
    hires_jacobians::Union{
        Nothing,
        Dict{<:AbstractStateVectorElement,<:Radiance}
        }
    "Radiance unit"
    radiance_unit::Unitful.Units
    "Solar scale array, to be multiplied into the solar irradiance"
    solar_scaler::AbstractVector

end

struct LSIRTMethod <: AbstractRTMethod

    "Options dictonary to control hi-res RT settings"
    high_options::Union{<:AbstractDict, Vector{<:AbstractDict}}
    "Boundaries of gas optical depth bins"
    bin_boundaries::Vector{Float64}
    "Total (gas) tau per spectral point"
    tau_gas::Vector{Float64}
    "√ξ coordinate per spectral point"
    ξ_sqrt::Vector{Float64}
    "Assigment of spectral points of to tau bins"
    tau_gas_bin_assignment::Vector{Integer}
    "Assigment of spectral points of to ξ bins"
    ξ_sqrt_bin_assignment::Vector{Integer}
    "Array to quickly check if this tau/ξ bin is used or not"
    used_bin::Array{Bool}
    "`MonochromaticRTMethod` object to contain the full low-streams results"
    monochromatic_RT::MonochromaticRTMethod
    "`MonochromaticRTMethod` objects to use for binned calculations at window center"
    RT_bin::MonochromaticRTMethod
    "`MonochromaticRTMethod` objects to use for binned calculations at window edge"
    RT_bin_edge::MonochromaticRTMethod
    "Low-res binned radiances"
    bin_rad_lo::Array{Radiance}
    "Low-res binned weighting functions"
    bin_wf_lo::Array{Vector{Radiance}}
    "High-res binned results"
    bin_rad_hi::Array{Radiance}
    "High-res binned weighting functions"
    bin_wf_hi::Array{Vector{Radiance}}
    "Low-res binned results for edge bin"
    bin_edge_rad_lo::Array{Radiance}
    "Low-res binned weighting functions for edge bin"
    bin_edge_wf_lo::Array{Vector{Radiance}}
    "High-res binned results for edge bin"
    bin_edge_rad_hi::Array{Radiance}
    "High-res binned weighting functions for edge bin"
    bin_edge_wf_hi::Array{Vector{Radiance}}
end