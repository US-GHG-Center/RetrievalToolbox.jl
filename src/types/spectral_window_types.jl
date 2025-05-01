"""
$(TYPEDFIELDS)

A type to hold a spectral window
"""
struct SpectralWindow{
    T<:AbstractFloat
    } <: AbstractSpectralWindow

    "Label for this spectral window"
    window_name::String
    "This is the 'what user wants' lower limit"
    ww_min::T
    "This is the 'what user wants' upper limit"
    ww_max::T

    "Wavelength or wavenumber grid (high-resolution) at the instrument"
    ww_grid::Vector{T}
    "Wavelength or wavenumber unit"
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    "Reference λ or ν from SV elements which have spectral dependence"
    ww_reference::T
    "Number of high-resolution spectral elements"
    N_hires::Int

    function SpectralWindow(
        window_name::String,
        ww_min::T,
        ww_max::T,
        ww_grid::Vector{T},
        ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
        ww_reference::T
    ) where {T}

        #=
            If this object is manually created, let's run some quick checks to make sure
            the values are sensible.
        =#

        # Check ww_min < ww_max
        if ww_unit isa Unitful.LengthUnits
            @assert ww_min <= ww_max "Minimum wavelength must be smaller than " *
                "maximum wavelength!"
        elseif ww_unit isa Unitful.WavenumberUnits
            @assert ww_min <= ww_max "Minimum wavenumber must be smaller than " *
                "maximum wavenumber!"
        end

        # Check ww array to be covering [ww_min, ww_max]
        @assert ww_min > ww_grid[1] "First element of hires array must be > min! " *
            "($(ww_min) must be > $(ww_grid[1]))"
        @assert ww_max < ww_grid[end] "Last element of hires array must be > max!"

        # ww should be sorted
        @assert issorted(ww_grid) "Hires array must be sorted!"

        return new{T}(
            window_name,
            ww_min,
            ww_max,
            ww_grid,
            ww_unit,
            ww_reference,
            length(ww_grid)
        )

    end


end

"""
Pretty printing for spectral window types

$(TYPEDSIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sw::SpectralWindow)

    println(io, "Spectral window: $(sw.window_name)")
    println(io, "Boundaries: $(sw.ww_min) to $(sw.ww_max), in $(sw.ww_unit).")
    println(io, "N = $(sw.N_hires)")

end

"""
Brief pretty printing for spectral window types

$(TYPEDSIGNATURES)
"""
function show(io::IO, sw::SpectralWindow)

    print(io, "SpectralWindow: $(sw.window_name)")

end


struct BinnedSpectralWindow{
    T1<:AbstractFloat,
    T2<:AbstractFloat,
    U<:AbstractSpectralWindow
    } <: AbstractSpectralWindow

    window_name::String
    "Wavelength or wavenumber grid (high-resolution) at the instrument"
    ww_grid::Vector{T1}
    "Wavelength or wavenumber unit"
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
    "Reference λ or ν from SV elements which have spectral dependence"
    ww_reference::T2
    "Number of high-resolution spectral elements"
    N_hires::Int
    "Reference to the original spectral window"
    original_window::U
    "A -special- spectral index, from which to e.g. draw scattering properties from."
    spectral_idx::Integer

end

"""
Brief pretty printing for binned spectral window types

$(TYPEDSIGNATURES)
"""
function show(io::IO, sw::BinnedSpectralWindow)

    print(io, "BinnedSpectralWindow $(sw.window_name) [$(sw.original_window.window_name)]")

end