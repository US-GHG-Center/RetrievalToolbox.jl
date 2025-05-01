"""
Returns the name of this ZeroLevelOffset state vector element
as a string.

$(SIGNATURES)

"""
function get_name(sve::ZeroLevelOffsetPolynomialSVE)
    return "ZeroLevelOffsetPolynomialSVE ($(sve.coefficient_order)) " *
        "[$(sve.swin.window_name)]"
end

"""
Pretty printing for ZeroLevelOffsetSVE types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sve::ZeroLevelOffsetPolynomialSVE)

    println(io, "Zero level offset polynomial SVE")
    println(io, "Coefficient order: $(sve.coefficient_order)")
    println(io, "First guess:       $(sve.first_guess)")
    println(io, "Prior value:       $(sve.prior_value)")
    println(io, "Prior covariance:  $(sve.prior_covariance)")
    println(io, "Current value:     $(sve.iterations[end])")

end

"""
Brief pretty printing for SurfacePressureSVE

$(SIGNATURES)
"""
function show(io::IO, sve::ZeroLevelOffsetPolynomialSVE)

    print(io, "$(get_name(sve)): $(sve.iterations[end])")

end

function calculate_jacobian_before_isrf(sve::ZeroLevelOffsetPolynomialSVE)

    # Partial derivative of forward model w.r.t. ZLO
    # does not need calculation before ISRF application.
    return false

end

"""
    Calculates the per-sample zero-level offset as given by the state vector element
    `sve` and indexed by the dispersion `dispersion`.
"""
function calculate_zlo!(
    output::AbstractVector,
    sve::ZeroLevelOffsetPolynomialSVE,
    dispersion::AbstractDispersion
)

    # This is the spectral window which the SVE references
    swin = sve.swin
    # These are the dispersion/detector indices for this spectral window
    indices = dispersion.index

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

        # Back convert to required units, and add to radiance
        output[idx] += this_value * this_unit_conversion

    end
end


"""
Calculates the Jacobian for a `ZeroLevelOffsetPolynomialSVE`

$(TYPEDSIGNATURES)

# Details

The ZLO radiance offset at wavelength of spectral sample `s` (``λ_s``) is generally
calculated as (and analogously for wavenumbers, ``ν_s``)
```math
\\sum_i c_i \\cdot (λ_s - λ_\\text{ref})^i
```
where i runs from 0 up to O - 1, where O is the order of the polynomial. See also
`apply_radiance_correction!` for this SVE type.

The partial derivative of the radiance ``I`` with respect to the polynomial coefficient
``c_i`` is then
```math
∂I/∂c_i = (λ_s - λ_\\text{ref})^i.
```

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
