function get_scattering_index(swin::SpectralWindow)

    return Int(round(swin.N_hires // 2))

end

function get_scattering_index(swin::BinnedSpectralWindow)

    return 1

end


"""
$(TYPEDSIGNATURES)

Procudes a `SpectralWindow` object given a name, min and max window edges, as well as a
spectral reference, a window buffer length, spectral units and a reference ABSCO object.
The resuling spectral grid will match the underlying ABSCO grid, such that optical
property calculations can be performed without having to do any spectral interpolation of
the spectroscopy.

Note: This function currently extends the high-res grid by a few points to make the total
high-res grid divisible by 4 to allow for some low-level optimization.
"""
function spectralwindow_from_ABSCO(
    name,
    ww_min::AbstractFloat,
    ww_max::AbstractFloat,
    ww_reference::AbstractFloat,
    buffer::AbstractFloat,
    absco::ABSCOSpectroscopy,
    ww_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits};
    skip=1
    )


    # First, look if the requested buffer puts us beyond the ABSCO length


    # Check which indices we need from the ABSCO table
    ww_start, ww_stop = searchsortedfirst.(
        Ref(absco.ww * absco.ww_unit),
        [ww_min - buffer, ww_max + buffer] .* ww_unit
    )

    # Add points so that the total number of high-res points
    # is divisible by 4.
    # Many loop optimizations might work a little faster if we
    # can divide the number of spectral points by e.g. 4 or 8.
    Nround = 4

    remainder = mod(length(ww_start:skip:ww_stop), Nround)
    ww_stop += Nround - remainder

    # Make sure the wavelength boundary does not exceed the coverage
    # of the ABSCO itself.
    ww_stop = min(ww_stop, length(absco.ww))

    # We have to perform a unit conversion here!
    ww_unit_fac = 1.0 * absco.ww_unit |> ww_unit |> ustrip

    ww = absco.ww[ww_start:skip:ww_stop] * ww_unit_fac

    return SpectralWindow(
        name,
        ww_min,
        ww_max,
        ww,
        ww_unit,
        ww_reference,
    )

end