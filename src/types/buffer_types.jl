"""
$(TYPEDFIELDS)

An instrument buffer type to hold pre-allocated arrays for use in instrument-related
functions such as the application of ISRFs.
"""
struct InstrumentBuffer{T<:AbstractFloat}

    tmp1::Vector{T}
    tmp2::Vector{T}
    low_res_output::Vector{T}

end

"""
$(TYPEDFIELDS)

An RT buffer type for scalar radiances, representing, for example, an instrument that is
capable of detecting total intensity only, which is what most remote sensing instruments
do. Note that this is *different* from the instrument having polarization sensitivity. The
role of this buffer is to provide pre-allocated vectors to hold the various radiance and
Jacobian values *after* potential instrument models are applied.

# Notes

The `radiance_unit` field is a `Union` of some arbitrary unit and a generic number type.
This is to allow e.g. DOAS-type retrievals, where the measured radiance is considered as a
ratio of some sort, and thus does not have a physical unit. Users are meant to then supply
`u"1"` as the `radiance_unit`. **Warning:** if users supply any number, such as
`radiance_unit=4.5`, expect that to act as a scaling factor by which state vector elements
will be scaled, as to bring a radiance-valued state vector element (for example,
`ZeroLevelOffsetPolynomialSVE`) to the same unit as the RT Buffer. It is advised to not do
so.

## Jacobians

At-instrument Jacobians are stored inside a Dict that map every state vector element to
a `ScalarRadiance` object. Users can also instantiate a `ScalarRTBuffer` object without
Jacobians, by simply passing `nothing`.

## Indices

Indices allow users to connect spectral windows and positions inside the `radiance` vector
to ensure the correct relationship between them. For example, if `swin` is the spectral
window, and the `ScalarRTBuffer` is `rt_buf`, the at-detector radiances should be stored
at `rt_buf.radiance.I[rt_buf.indices[swin]]`. Note that radiances are **not automatically
copied there, but must be copied by users themselves**.
"""
struct ScalarRTBuffer{T1<:AbstractFloat, T2<:Integer} <: AbstractRTBuffer

    "Dict to map spectral windows to dispersion objects"
    dispersion::Dict{<:AbstractSpectralWindow,<:AbstractDispersion}
    "Scalar radiance vector representing the at-instrument radiance"
    radiance::ScalarRadiance{T1, Vector{T1}}
    "At-instrument level jacobians, either a Dict that maps state vector elements to a
    `ScalarRadiance`, or nothing."
    jacobians::Union{
        Nothing,
        Dict{<:AbstractStateVectorElement, ScalarRadiance{T1, Vector{T1}}}
    }
    "Dict to map spectral windows to positions in `radiance`."
    indices::Dict{<:AbstractSpectralWindow, Vector{T2}}
    "Radiance unit of `radiance`."
    radiance_unit::Union{Unitful.Units, Unitful.Number}

end

"""
$(TYPEDFIELDS)

An RT buffer type for vector radiances, representing an instrument that is capable of
detecting the polarization state of light. Chances are *this is not what you need*, since
almost all known remote sensing devices measure intensity only - even if they have some
polarization sensitivity.
"""
struct VectorRTBuffer{T1<:AbstractFloat, T2<:Integer} <: AbstractRTBuffer

    dispersion::Dict{<:AbstractSpectralWindow,<:AbstractDispersion}
    radiance::VectorRadiance{T1, Matrix{T1}}
    jacobians::Union{
        Nothing,
        Dict{<:AbstractStateVectorElement, VectorRadiance{T1, Matrix{T1}}}
    }
    indices::Dict{<:AbstractSpectralWindow, Vector{T2}}
    radiance_unit::Union{Unitful.Units, Unitful.Number}

end

"""
$(TYPEDFIELDS)

A buffer for use in `EarthAtmosphere` type simulations and/or retrievals, containing other
buffers as well as the `EarthScene` to describe the location and the atmospheric state.

# Details

This buffer is needed to represent a full single-scene modeled measurement at one time and
location. It is highly recommended to not instantiate this buffer object manually, but to
use the corresponding `EarthAtmosphereBuffer` function that helps produce this. Manually
creating this object carries the risk of some of the dictionaries having incorrect keys.
"""
struct EarthAtmosphereBuffer <: AbstractAtmosphereBuffer

    """A vector of type `AbstractSpectralWindow` to hold all spectral windows used in this
       buffer"""
    spectral_window::Vector{<:AbstractSpectralWindow}
    "An `EarthScene`"
    scene::EarthScene
    "A Dict to hold the optical properties attached to an `AbstractSpectralWindow`"
    optical_properties::Dict{<:AbstractSpectralWindow, EarthAtmosphereOpticalProperties}
    "A Dict to hold the RT buffer attached to an `AbstractSpectralWindow`"
    rt::Dict{<:AbstractSpectralWindow, <:AbstractRTMethod}
    "An RT buffer object to hold retrieval-wide radiances and Jacobians"
    rt_buf::AbstractRTBuffer
    "A convolution buffer object to help with the ISRF application"
    inst_buf::InstrumentBuffer
end