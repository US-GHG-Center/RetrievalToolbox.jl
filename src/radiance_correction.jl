"""
$(TYPEDSIGNATURES)

Generic function to apply a radiance correction to an at-instrument spectrum, for example
stray light or a zero-level offset. The generic implementation does nothing, and users
should implement their own in accordance with the state vector element `sve`.
"""
function apply_radiance_correction!(
    rt_buf::AbstractRTBuffer,
    sve::AbstractStateVectorElement
    )

    # Default implementation, which does nothing
    return nothing

end



"""
$(TYPEDSIGNATURES)

Radiance "correction" due to a zero-level offset as given by the
`ZeroLevelOffsetPolynomialSVE` object. The spectrally-dependent radiance is added to the
total radiance, and the intensity component only.
"""
function apply_radiance_correction!(
    rt_buf::AbstractRTBuffer,
    sve::ZeroLevelOffsetPolynomialSVE
    )

    # Perform the ZLO calculation
    calculate_zlo!(
        rt_buf.radiance.I,
        rt_buf.radiance_unit,
        sve,
        rt_buf.dispersion[sve.swin],
        rt_buf.indices[sve.swin]
    )

end

"""
$(TYPEDSIGNATURES)

Calculates the per-sample zero-level offset as given by the state vector element
`sve` and indexed by the dispersion `dispersion`. This function is called by
`apply_radiance_correction!` and likely does not need to be done by the user.
"""
function calculate_zlo!(
    output::AbstractVector,
    radiance_unit,
    sve::ZeroLevelOffsetPolynomialSVE,
    dispersion::AbstractDispersion,
    indices::AbstractVector
)

    # This is the spectral window which the SVE references
    swin = sve.swin

    # Convert the ZLO to units specified in the radiance buffer.

    # Figure out the physical unit before the loop
    # NOTE that doing this calculation within the loop is very performance-draining!
    this_unit = get_unit(sve) * dispersion.ww_unit^sve.coefficient_order

    # Find out the conversion factor needed to bring (sve)*Δww^coefficient_order
    # to the same units as the radiance buffer. The branch below is needed because
    # `ustrip` does not work if the first argument is a `Number`.

    if radiance_unit isa Number
        this_unit_conversion = radiance_unit * this_unit
    else
        this_unit_conversion = ustrip(radiance_unit, 1.0 * this_unit)
    end

    # Δww is in units of the dispersion^coefficient_order now
    ww_delta = (
        dispersion.ww .-
        ustrip(dispersion.ww_unit, swin.ww_reference * swin.ww_unit)
    ) .^ sve.coefficient_order

    sve_value = get_current_value(sve)

    @turbo for (i, idx) in enumerate(indices)

        # Calculate c_i * (ww - ww_reference)^o (with unit)
        this_value = sve_value * ww_delta[i]

        # Back convert to required units, and add to radiance
        output[idx] += this_value * this_unit_conversion

    end
end


"""
$(TYPEDSIGNATURES)

Calculates the Jacobian for a `ZeroLevelOffsetPolynomialSVE`

# Details

The ZLO radiance offset at wavelength of spectral sample `s` (``λ_s``) is generally
calculated as (and analogously for wavenumbers, ``ν_s``)
```math
\\sum_{i=0}^{O-1} c_i \\cdot (λ_s - λ_\\text{ref})^i
```
where i runs from 0 up to O - 1, where O is the order of the polynomial. See also
`apply_radiance_correction!` for this SVE type.

The partial derivative of the radiance ``I`` with respect to the polynomial coefficient
``c_i`` is then
```math
∂I/∂c_i = (λ_s - λ_\\text{ref})^i.
```

The reference spectral point (``λ_\\text{ref}`` or ``ν_\\text{ref}``) is given by the
spectral window object that is attached to the `sve` as `sve.swin.ww_reference`. Unit
differences between the `rt_buf` RT buffer and the unit of this `sve` are explicitly
taken into account, so this `sve` may have any compatible radiance units as long as the
converstion to the `rt_buf.radiance_unit` is valid.
"""
function calculate_jacobian!(
    rt_buf::AbstractRTBuffer,
    sve::ZeroLevelOffsetPolynomialSVE
    )

    # Grab the dispersion attached to the spectral window
    # belonging to this SVE
    dispersion = rt_buf.dispersion[sve.swin]

    # This is where the result goes, important to have the
    # view to the Jacobian array!
    result = @views rt_buf.jacobians[sve].I[rt_buf.indices[sve.swin]]

    if sve.coefficient_order == 0
        # For order 0, the Jacobian is simply 1.0, plus accounting
        # for unit differences between SVE and RT buffer
        unit_fac = (1.0 * rt_buf.radiance_unit / get_unit(sve))

        for i in eachindex(result)
             result[i] = unit_fac
        end
    else
        # Take the dispersion ww
        @views result[:] = dispersion.ww[:]

        # Subtract the reference ww, but account for a possible
        # difference in units from SVE to dispersion
        @views result[:] .-= ustrip(dispersion.ww_unit,
            sve.swin.ww_reference * sve.swin.ww_unit)

        # Exponentiate
        for i in eachindex(result)
            result[i] = result[i] ^ sve.coefficient_order
        end

        # Account for unit change between SVE and radiances in
        # RT buffer
        @views result[:] .*=
            (1.0 * rt_buf.radiance_unit / get_unit(sve) /
            dispersion.ww_unit ^ sve.coefficient_order)

    end

end
