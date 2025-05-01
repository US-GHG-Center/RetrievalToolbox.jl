struct NoSurface <: AbstractSurface end


struct LambertianPolynomialSurface{T} <: AbstractSurface

    swin::AbstractSpectralWindow
    order::Int
    coefficients::Vector{T}

    function LambertianPolynomialSurface(
        swin::AbstractSpectralWindow,
        coefficients::Vector{T}) where {T}

        return new{T}(swin, length(coefficients) - 1, coefficients)

    end

end

"""
Type to hold a combination of BRDF Kernels

$(TYPEDFIELDS)
"""
struct BRDFSurface{T <: BRDFKernel} <: AbstractSurface

    kernels::Vector{T}

end

"""
Type that implements the Lambertian BRDF kernel for which the BRDF amplitude
can be a spectrally varying polynomial of arbitrary order.

$(TYPEDFIELDS)
"""
struct LambertianPolynomialKernel{T} <: BRDFKernel

    swin::AbstractSpectralWindow
    order::Int
    coefficients::Vector{T}
    asymmetry::T
    anisotropy::T

    function LambertianPolynomialKernel(
        swin::AbstractSpectralWindow,
        coefficients::Vector{T},
        ) where T

        return new{T}(swin, length(coefficients) - 1, coefficients)

    end

end

get_short_name(k::Type{LambertianPolynomialKernel}) = "LambertianKernel"

"""
$(TYPEDSIGNATURES)
Returns the string needed to register this surface type with XRTM
"""
get_XRTM_name(k::LambertianPolynomialKernel) = "lambertian"

"""
Type that implements the Rahman-Pinty-Verstraete BRDF kernel for which the BRDF amplitude
can be a spectrally varying polynomial of arbitrary order. The other two parameters,
asymmetry and anisotropy, are kept fixed throughout the spectral window which this surface
is attached to.

$(TYPEDFIELDS)
"""
struct RPVPolynomialKernel{T} <: BRDFKernel

    swin::AbstractSpectralWindow
    order::Int
    coefficients::Vector{T}
    hotspot::T
    asymmetry::T
    anisotropy::T

    function RPVPolynomialKernel(
        swin::AbstractSpectralWindow,
        coefficients::Vector{T},
        hotspot::T,
        asymmetry::T,
        anisotropy::T) where T

        return new{T}(swin, length(coefficients) - 1, coefficients,
            hotspot, asymmetry, anisotropy)

    end

end

get_short_name(k::Type{RPVPolynomialKernel}) = "RPVKernel"
"""
$(TYPEDSIGNATURES)
Returns the string needed to register this surface type with XRTM
"""
get_XRTM_name(k::RPVPolynomialKernel) = "rahman"
