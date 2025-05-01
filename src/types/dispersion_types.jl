"""


$(TYPEDFIELDS)

"""
struct SimplePolynomialDispersion{
    T1<:AbstractSpectralWindow,
    T2<:AbstractFloat,
    T3<:AbstractRange
    } <: AbstractDispersion

    "Spectral window"
    spectral_window::T1
    "Polynomial coefficients, ordered low to high"
    coefficients::Vector{T2}
    "Spectral samples within the current spectral window!"
    index::Vector{Int}
    "Wavelengths corresponding to the `index`"
    ww::Vector{T2}
    "The full detector grid corresponding to the L1B spectral samples"
    detector_samples::T3
    "Unitful unit of all wavelength and coefficient quantities"
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}
end


function SimplePolynomialDispersion(
    coeffs,
    detector_samples,
    spectral_window
    )

    # We demand the user supplies coefficients
    # with some length unit attached, which we then
    # strip.

    @assert coeffs isa Vector{<:Unitful.AbstractQuantity} "Coefficients must have units!"

    coeff_unit = unit(coeffs[1])
    sw_unit = spectral_window.ww_unit

    this_unit = eltype(coeffs)
    # Unit conversion happens here; so the dispersion coefficients can be
    # supplied in any length unit
    c = ustrip.(sw_unit, coeffs)
    full_grid = Polynomial(c).(detector_samples)

    # Decreasing or increasing?
    # (maybe the wavelength grid comes in decreasing order,
    #  due to the way the instrument optics work..)
    if full_grid[2] > full_grid[1]
        # Increasing ww
        idx_start = searchsortedfirst(full_grid,
                                      spectral_window.ww_min)
        idx_stop = searchsortedfirst(full_grid,
                                     spectral_window.ww_max) + 1
    else
        # Decreasing ww
        idx_start = searchsortedlast(full_grid,
                                     spectral_window.ww_max, rev=true)
        idx_stop = searchsortedlast(full_grid,
                                    spectral_window.ww_min, rev=true)
    end

    idx_start = max(idx_start, 1)
    idx_stop = min(idx_stop, length(detector_samples))

    return SimplePolynomialDispersion(
        spectral_window,
        c,
        collect(idx_start:idx_stop),
        full_grid[idx_start:idx_stop],
        detector_samples,
        coeff_unit
    )

end
