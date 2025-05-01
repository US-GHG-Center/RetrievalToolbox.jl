
function apply_radiance_correction!(
    rt_buf::AbstractRTBuffer,
    sve::AbstractStateVectorElement
    )

    # Default implementation, which does nothing
    return nothing

end



"""

$(TYPEDSIGNATURES)

"""
function apply_radiance_correction!(
    rt_buf::AbstractRTBuffer,
    sve::ZeroLevelOffsetPolynomialSVE
    )


    # This is the spectral window which the SVE references
    swin = sve.swin
    # This is the dispersion which the spectral window is assigned to
    dispersion = rt_buf.dispersion[swin]
    # These are the dispersion/detector indices for this spectral window
    indices = rt_buf.indices[swin]

    # Convert the ZLO to units specified in the
    # radiance buffer.

    # Figure out the physical unit before the loop
    # NOTE that doing this calculation within the loop is very performance-draining!
    this_unit = get_unit(sve) * dispersion.ww_unit^sve.coefficient_order
    # Find out the conversion factor needed to bring (sve)*Δww^coefficient_order
    # to the same units as the radiance buffer
    this_unit_conversion = ustrip(rt_buf.radiance_unit, 1.0 * this_unit)

    # Δww is in units of the dispersion^coefficient_order now
    ww_delta = (
        dispersion.ww .-
        ustrip(dispersion.ww_unit, swin.ww_reference * swin.ww_unit)
    ) .^ sve.coefficient_order

    @turbo for (i, idx) in enumerate(indices)

        # Calculate c_i * (ww - ww_reference)^o (with unit)
        this_value = get_current_value(sve) * ww_delta[i]

        # Back convert to required units, and add to intensity
        rt_buf.radiance.I[idx] += this_value * this_unit_conversion

    end

    # This should be short-hand!
    #rad = @views rt_buf.radiance.I[:][indices]
    #calculate_zlo!(rad, sve, rt_buf.dispersion[sve.swin])

end