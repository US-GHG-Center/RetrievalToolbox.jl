struct InstrumentBuffer{T<:AbstractFloat}

    tmp1::Vector{T}
    tmp2::Vector{T}
    low_res_output::Vector{T}

end

"""

    $(TYPEDFIELDS)

    # Notes
    Note that the `radiance_unit` field is a `Union` of some arbitrary unit and a generic
    number type. This is to allow e.g. DOAS-type retrievals, where the measured radiance 
    is considered as a ratio of some sort, and thus does not have a physical unit. Users
    are meant to then supply `u"1"` as the `radiance_unit`. **Warning:** if users supply
    any number, such as `radiance_unit=4.5`, expect that to act as a scaling factor by
    which state vector elements will be scaled, as to bring a radiance-valued state vector
    element (for example, `ZeroLevelOffsetPolynomialSVE`) to the same unit as the RT
    Buffer. 

"""
struct ScalarRTBuffer{T1<:AbstractFloat, T2<:Integer} <: AbstractRTBuffer

    dispersion::Dict{<:AbstractSpectralWindow,<:AbstractDispersion}
    radiance::ScalarRadiance{T1, Vector{T1}}
    jacobians::Union{
        Nothing,
        Dict{<:AbstractStateVectorElement, ScalarRadiance{T1, Vector{T1}}}
    }
    indices::Dict{<:AbstractSpectralWindow, Vector{T2}}
    radiance_unit::Union{Unitful.Units, Unitful.Number}

end

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
A buffer for use in EarthAtmosphere type simulations and/or retrievals.

$(TYPEDFIELDS)

# Details

"""
struct EarthAtmosphereBuffer <: AbstractAtmosphereBuffer

    """A vector of type `AbstractSpectralWindow` to hold all
       spectral windows used in this buffer"""
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

"""
we haven't really used this one much so far..
"""
struct OEBuffer{T<:AbstractFloat}

    K::Matrix{T} # Jacobian
    Se::Diagonal{T} # Error covariance
    Sa::Matrix{T} # Prior covariance
    Shat::Matrix{T} # Posterior covariance
    G::Matrix{T} # Gain matrix
    AK::Matrix{T} # Averaging kernel matrix

    function OEBuffer(N_spec, N_sv, ::Type{T}) where T

        # For now, we only support diagonal error covaraince matrices
        Se = Diagonal{T}(undef, N_spec)
        @views Se[:,:] .= 0

        return new{T}(
            zeros(T, N_spec, N_sv), # K
            Se, # Se
            zeros(T, N_sv, N_sv), # Sa
            zeros(T, N_sv, N_sv), # Shat
            zeros(T, N_sv, N_spec), # G
            zeros(T, N_sv, N_sv), # AK
        )

    end


end
