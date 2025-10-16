"""
$(TYPEDFIELDS)

Type to hold surface-leaving SIF radiance
"""
mutable struct SIFradiance{T}

    radiance_at_reference :: T
    radiance_unit

    ww_reference :: T
    ww_unit :: Union{Unitful.WavenumberUnits, Unitful.LengthUnits}

    ww_waveform :: Vector
    waveform :: Vector

end

"""
$(TYPEDSIGNATURES)

Pretty printing for SIF radiance types
"""
function show(io::IO, ::MIME"text/plain", sif::SIFradiance)

    println(io, "SIF radiance at $(sif.ww_reference * sif.ww_unit): " *
     "$(sif.radiance_at_reference * sif.radiance_unit)")

end

"""
$(TYPEDSIGNATURES)

Pretty printing for SIF radiance types
"""
function show(io::IO, sif::SIFradiance)

    println(io, "SIF radiance at $(sif.ww_reference * sif.ww_unit): " *
     "$(sif.radiance_at_reference * sif.radiance_unit)")

end



function SIFradiance(
    radiance_at_reference,
    radiance_unit,
    ww_reference,
    ww_unit
)

    # Load the reference waveform
    fname = joinpath(@__DIR__, "..", "..", "data", "SIF", "SIF_waveform.csv")
    csv = CSV.File(fname, header=false)

    ww_waveform = csv.Column1 * u"nm" # nanometer is hardcoded
    waveform = csv.Column2 # radiance waveform is relative, no units


    return SIFradiance(
        radiance_at_reference,
        radiance_unit,
        ww_reference,
        ww_unit,
        ww_waveform,
        waveform
    )


end