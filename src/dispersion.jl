"""
Creates a spectral grid from a `SimplePolynomialDispersion` object
and returns the indices and the corresponding spectral grid.

$(TYPEDSIGNATURES)
"""
function calculate_grid_from_dispersion(
    disp::SimplePolynomialDispersion,
    spectral_window::AbstractSpectralWindow
    )

    full_grid = Polynomial(disp.coefficients).(disp.detector_samples)

    idx_start = searchsortedfirst(full_grid, spectral_window.ww_min)
    idx_stop = searchsortedfirst(full_grid, spectral_window.ww_max) + 1

    return idx_start:idx_stop, full_grid[idx_start:idx_stop]

end

"""
Updates the dispersion object `disp` given a state vector `sv`,
also considering an optional instrument Doppler factor. The
instrument Doppler factor reflects the relative movement between
(for example), the measurement location on Earth and the spacecraft.

$(TYPEDSIGNATURES)

"""
function update_dispersion!(
    disp::SimplePolynomialDispersion,
    sv::AbstractStateVector
    )

    # Determines whether `disp` actually needs updating via new values from
    # the state vector.
    needs_updating = false

    # Repace dispersion coefficients from state vector, and
    # be aware of units!
    for (idx, sve) in StateVectorIterator(sv, DispersionPolynomialSVE)
            # Array index that this coefficient order (0-based) belongs to..
            o = sve.coefficient_order + 1

        # Only update the dispersion coefficients
        if sve.dispersion === disp

            val = ustrip(disp.ww_unit, get_current_value_with_unit(sve))

            disp.coefficients[o] = val

            needs_updating = true

        end


    end

    # Update dispersion wavelengths/wavenumbers and indices!
    if needs_updating
        @debug "Dispersion attached to $(disp.spectral_window) updated!"
        update_dispersion!(disp)
    end

end





"""
Update the dispersion according to its coefficients.

$(TYPEDSIGNATURES)

# Details

During retrievals, the dispersion as given by polynomial coefficients, might change from
the first-guess values. This function provides a quick and low-allocating way of updating
the spectral samples according to the coefficients stored within. The user first must
manually update the coefficients

```
  disp.coefficients = ...
```

and then call this function in order to update `disp.ww` and `disp.index`.
"""
function update_dispersion!(
    disp::SimplePolynomialDispersion
    )

    # Rebind for convenience
    window = disp.spectral_window

    # Empty out index
    empty!(disp.index)
    # Empty out wavelengths
    empty!(disp.ww)

    # Build polynomial from coefficients
    p = Polynomial(disp.coefficients)
    full_grid = p.(disp.detector_samples)

    idx_start = -1
    idx_stop = -1

    # Decreasing or increasing?
    # (maybe the wavelength grid comes in decreasing order,
    #  due to the way the instrument optics work..)
    if full_grid[2] > full_grid[1]
        # Increasing wavelength
        idx_start = searchsortedfirst(full_grid, window.ww_min)
        idx_stop = searchsortedfirst(full_grid, window.ww_max) + 1
    else
        # Decreasing wavelength
        idx_start = searchsortedlast(full_grid, window.ww_max, rev=true)
        idx_stop = searchsortedlast(full_grid, window.ww_min, rev=true)
    end

    idx_start = max(idx_start, 1)
    idx_stop = min(idx_stop, length(disp.detector_samples))

    # Fill in new values
    for i in idx_start:idx_stop
        append!(disp.index, disp.detector_samples[i])
        append!(disp.ww, full_grid[i])
    end

end



"""
$(SIGNATURES)

Calculates the analytic dispersion polynomial Jacobian

# Details

Analytic computation of ∂I/dc_i, where `c_i` is the coefficient of the dispersion
polynomial order `i`, i = 0 meaning a flat dispersion, i = 1 is a linearly varying one
etc.
"""
function calculate_dispersion_polynomial_jacobian!(
    output::Vector,
    inst_buf::InstrumentBuffer,
    sve::DispersionPolynomialSVE,
    ISRF::TableISRF,
    disp::AbstractDispersion,
    data::AbstractVector,
    swin::AbstractSpectralWindow;
    doppler_factor=0.0
    )

    @assert length(data) == length(swin.ww_grid) (
        "Hi-res vector and hires vector grid must be same size!"
    )

    # Zero-out the output
    output[:] .= 0

    # Depending on λ/ν, we use a different effective Doppler formula
    if disp.ww_unit isa Unitful.LengthUnits
        doppler_term = 1.0 / (1.0 + doppler_factor)
    elseif disp.ww_unit isa Unitful.WavenumberUnits
        doppler_term = 1.0 + doppler_factor
    end

    # Cast a view into a temp array for use as the Δλ array
    this_ww_delta = @view inst_buf.tmp2[1:size(ISRF.ww_delta, 1)]

    # Calculate the unit conversion factor between ISRF wavelength and
    # spectral window wavelength
    unit_fac = 1.0 * ISRF.ww_delta_unit / swin.ww_unit |> NoUnits

    for i in eachindex(disp.index)

        if i == 1
            continue
        end

        this_wl = disp.ww[i] * doppler_term
        this_l1b_idx = disp.index[i - 1]

        tmp1 = this_l1b_idx ^ sve.coefficient_order

        idx_first = searchsortedfirst(
            swin.ww_grid,
            this_wl + unit_fac * ISRF.ww_delta[1, this_l1b_idx]
        )
        idx_last = searchsortedfirst(
            swin.ww_grid,
            this_wl + unit_fac * ISRF.ww_delta[end, this_l1b_idx]
        )

        if idx_first <= 1
            continue
        end

        if idx_last > length(swin.ww_grid)
            continue
        end


        this_ISRF = @view inst_buf.tmp1[1:idx_last - idx_first + 1]
        @views this_ISRF[:] .= 0.0

        @views inst_buf.tmp2[:] .= 0.0
        @views @. this_ww_delta[:] = (
            unit_fac * ISRF.ww_delta[:, this_l1b_idx] + this_wl
        )

        this_relative_response = @view ISRF.relative_response[:, this_l1b_idx]
        this_hires_wl = @view swin.ww_grid[idx_first:idx_last]

        # Interpolate ISRF to grid!
        pwl_value_1d!(
            this_ww_delta,
            this_relative_response,
            this_hires_wl,
            this_ISRF
        )

        # Re-zero tmp2
        @views inst_buf.tmp2[:] .= 0.0

        # Calculate dISRF / dwavelength
        for i in 1:length(this_ISRF) - 1
            inst_buf.tmp2[i] = this_ISRF[i + 1] - this_ISRF[i]
            inst_buf.tmp2[i] /= swin.ww_grid[idx_first + i + 1] - swin.ww_grid[idx_first + i]
        end

        dISRF_dwavelength = @view inst_buf.tmp2[1:idx_last - idx_first]
        this_data = @view data[idx_first:idx_last - 1]

        # Not sure why the minus here is needed (it is!), maybe check the math
        # again at some point.
        output[this_l1b_idx] = -avx_dot(
            dISRF_dwavelength,
            this_data
        ) * tmp1 / sum(this_ISRF)

    end

    return true

end
