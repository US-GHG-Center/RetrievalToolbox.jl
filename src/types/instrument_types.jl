struct TableISRF{T} <: AbstractISRF

    "Δλ array on which the ISRF is defined (Δλ or Δν, spectral sample)"
    ww_delta::Array{T, 2}
    "Units for the `ww_delta` field. The ISRF then has units 1 / `ww_delta`"
    ww_delta_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    "ISRF relative response (Δλ, spectral sample)"
    relative_response::Array{T, 2} # ISRF delta lambda, spectral index

end


"""
Pretty printing for TableISRF types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", isrf::TableISRF)

    println(io, "TableISRF")
    println(io, "N Δλ, N spectral: $(size(isrf.ww_delta, 1)), $(size(isrf.ww_delta, 2))")

end

"""
Brief pretty printing for TableISRF types

$(SIGNATURES)
"""
function show(io::IO, isrf::TableISRF)

    print(io, "TableISRF: ($(size(isrf.ww_delta, 1)), $(size(isrf.ww_delta, 2)))")

end


"""
A type to store a Gaussian ISRF that is defined only by its FWHM value
"""
struct GaussISRF{T} <: AbstractISRF
    "Full width at half of the maximum"
    FWHM::T
    "Unit belonging to FWHM"
    FWHM_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
end

"""
Pretty printing for GaussISRF types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", isrf::GaussISRF)

    σ = FWHM_to_sigma(isrf.FWHM)
    println(io, "GaussISRF")
    println(io, "FWHM = $(isrf.FWHM), σ = $(σ)")

end

"""
Brief pretty printing for TableISRF types

$(SIGNATURES)
"""
function show(io::IO, isrf::GaussISRF)

    σ = FWHM_to_sigma(isrf.FWHM)
    print(io, "GaussISRF: FWHM = $(isrf.FWHM), σ = $(σ)")

end