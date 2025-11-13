"""
$(TYPEDFIELDS)

Type to hold surface-leaving SIF radiance
"""
mutable struct SIFRadiance{T} <: AbstractSIFRadiance

    radiance_at_reference :: T
    radiance_unit :: Unitful.Units

    ww_reference :: T
    ww_unit :: Union{Unitful.WavenumberUnits, Unitful.LengthUnits}

    ww_waveform :: Vector
    waveform :: Vector

    itp :: Interpolations.AbstractInterpolation

end

"""
$(TYPEDSIGNATURES)

Pretty printing for SIF radiance types
"""
function show(io::IO, ::MIME"text/plain", sif::SIFRadiance)

    println(io, "SIF radiance at $(sif.ww_reference * sif.ww_unit): " *
     "$(sif.radiance_at_reference * sif.radiance_unit)")

end

"""
$(TYPEDSIGNATURES)

Pretty printing for SIF radiance types
"""
function show(io::IO, sif::SIFRadiance)

    println(io, "SIF radiance at $(sif.ww_reference * sif.ww_unit): " *
     "$(sif.radiance_at_reference * sif.radiance_unit)")

end


"""
$(TYPEDSIGNATURES)

Constructor to produce a new SIFradiance object.
"""
function SIFRadiance(
    radiance_at_reference,
    radiance_unit,
    ww_reference,
    ww_unit::Union{Unitful.WavenumberUnits, Unitful.LengthUnits}
)

    # Load the reference waveform
    fname = joinpath(@__DIR__, "..", "..", "data", "SIF", "SIF_waveform.csv")
    csv = CSV.File(fname, header=false)

    ww_waveform = csv.Column1 * u"nm" # nanometer is hardcoded
    waveform = csv.Column2 # radiance waveform is relative, no units

    itp = linear_interpolation(ww_waveform, waveform, extrapolation_bc = Throw())

    return SIFRadiance(
        radiance_at_reference,
        radiance_unit,
        ww_reference,
        ww_unit,
        ww_waveform,
        waveform,
        itp
    )

end