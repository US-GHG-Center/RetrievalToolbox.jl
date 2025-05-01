struct OCOInstrument <: AbstractInstrument

end

struct GeoCarbInstrument <: AbstractInstrument

end

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
function show(io::IO, ::MIME"text/plain", ils::TableISRF)

    println(io, "TableISRF")
    println(io, "N Δλ, N spectral: $(size(ils.ww_delta, 1)), $(size(ils.ww_delta, 2))")

end

"""
Brief pretty printing for TableISRF types

$(SIGNATURES)
"""
function show(io::IO, isrf::TableISRF)

    print(io, "TableISRF: ($(size(isrf.ww_delta, 1)), $(size(isrf.ww_delta, 2)))")

end
