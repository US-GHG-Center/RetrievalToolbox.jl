"""
$(TYPEDSIGNATURES)

Creates a spectral grid from a `SimplePolynomialDispersion` object and returns the indices
and the corresponding spectral grid.
"""
function calculate_grid_from_dispersion(
    disp::SimplePolynomialDispersion,
    window::AbstractSpectralWindow
    )

    full_grid = Polynomial(disp.coefficients).(disp.detector_samples)

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

    return idx_start:idx_stop, full_grid[idx_start:idx_stop]

end

"""
$(TYPEDSIGNATURES)

Updates the dispersion object `disp` given a state vector `sv`, also considering an
optional instrument Doppler factor. The instrument Doppler factor reflects the relative
movement between (for example), the measurement location on Earth and the spacecraft.
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
        @debug "[DISP] Dispersion attached to $(disp.spectral_window) updated!"
        update_dispersion!(disp)
    end

end





"""
$(TYPEDSIGNATURES)

Updates the dispersion according to its coefficients.

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
$(TYPEDSIGNATURES)

Calculates the analytic dispersion polynomial Jacobian for a `TableISRF` instrument
response function.

# Details

Analytic computation of ∂I/∂c_i, where `c_i` is the coefficient of the dispersion
polynomial order `i`, i = 0 meaning a flat dispersion, i = 1 is a linearly varying one
etc.
"""
function calculate_dispersion_polynomial_jacobian!(
    inst_buf::InstrumentBuffer,
    sve::DispersionPolynomialSVE,
    ISRF::TableISRF,
    data::AbstractVector;
    doppler_factor=0.0,
    Ndelta::Union{Nothing, Number}=nothing,
    )

    disp = sve.dispersion
    swin = disp.spectral_window

    # Estimate the best stepsize to use for partial derivative
    # ∂ISRF / ∂λ (or ∂ISRF / ∂ν):

    if isnothing(Ndelta)
        NΔ = (
            minimum(diff(disp.ww)) / minimum(diff(swin.ww_grid))
        ) |> round |> Int
    else
        NΔ = Ndelta
    end

    @debug "[DISPERSION] NΔ = $(NΔ)"

    # Call low-level high performance function
    return _calculate_dispersion_polynomial_jacobian_table!(
        inst_buf,
        sve.coefficient_order,
        ISRF,
        NΔ,
        swin.ww_unit,
        swin.ww_grid,
        disp.ww_unit,
        disp.ww,
        disp.index,
        data, # this should be called radiance or jacobian really
        doppler_factor
    )

end

function _calculate_dispersion_polynomial_jacobian_table!(
    inst_buf::InstrumentBuffer,
    order::Integer,
    ISRF::TableISRF,
    NΔ::Integer,
    hires_unit,
    hires_ww,
    lores_unit,
    lores_ww,
    lores_index,
    data,
    doppler_factor
)
    @assert NΔ > 0 "NΔ must be > 0! (NΔ = $(NΔ))"
    @assert length(data) == length(hires_ww) (
        "Hi-res vector and hires vector grid must be same size!"
    )

    output = inst_buf.low_res_output
    # Zero-out the output
    output[:] .= 0

    # Depending on λ/ν, we use a different effective Doppler formula
    if lores_unit isa Unitful.LengthUnits
        doppler_term = 1.0 / (1.0 + doppler_factor)
    elseif lores_unit isa Unitful.WavenumberUnits
        doppler_term = 1.0 + doppler_factor
    end

    # Cast a view into a temp array for use as the Δλ array
    this_ww_delta = @view inst_buf.tmp2[1:size(ISRF.ww_delta, 1)]

    # Calculate the unit conversion factor between ISRF wavelength and
    # spectral window wavelength
    unit_fac = 1.0 * ISRF.ww_delta_unit / hires_unit |> NoUnits

    for i in eachindex(lores_index)

        if i == 1
            continue
        end

        this_ww = lores_ww[i] * doppler_term
        this_l1b_idx = lores_index[i - 1]

        tmp1 = this_l1b_idx ^ order

        idx_first = searchsortedfirst(
            hires_ww,
            this_ww + unit_fac * ISRF.ww_delta[1, this_l1b_idx]
        )
        idx_last = searchsortedfirst(
            hires_ww,
            this_ww + unit_fac * ISRF.ww_delta[end, this_l1b_idx]
        )

        if idx_first <= 1
            continue
        end

        if idx_last > length(hires_ww)
            continue
        end


        this_ISRF = @view inst_buf.tmp1[1:idx_last - idx_first + 1]
        @views this_ISRF[:] .= 0

        @views inst_buf.tmp2[:] .= 0
        @views @. this_ww_delta[:] = (
            unit_fac * ISRF.ww_delta[:, this_l1b_idx] + this_ww
        )

        this_relative_response = @view ISRF.relative_response[:, this_l1b_idx]
        this_hires_wl = @view hires_ww[idx_first:idx_last]

        # Interpolate ISRF to grid!
        pwl_value_1d!(
            this_ww_delta,
            this_relative_response,
            this_hires_wl,
            this_ISRF
        )

        # Re-zero tmp2
        @views inst_buf.tmp2[:] .= 0

        # Calculate dISRF / dwavelength
        # (make this a central difference..)
        for i in 1:length(this_ISRF) - NΔ
            inst_buf.tmp2[i] = this_ISRF[i + NΔ] - this_ISRF[i]
            inst_buf.tmp2[i] /= this_hires_wl[i + NΔ] - this_hires_wl[i]
        end

        dISRF_dww = @view inst_buf.tmp2[1:length(this_ISRF) - 1]
        this_data = @view data[idx_first:idx_last - 1]

        # Not sure why the minus here is needed (it is!), maybe check the math
        # again at some point.
        output[this_l1b_idx] = -avx_dot(
            dISRF_dww,
            this_data
        ) * tmp1 / sum(this_ISRF)

    end

    return true

end

function calculate_dispersion_polynomial_jacobian!(
    inst_buf::InstrumentBuffer,
    sve::DispersionPolynomialSVE,
    ISRF::GaussISRF,
    data;
    doppler_factor=0.0,
    extend=4.0
    )

    # Zero out the output
    inst_buf.low_res_output[:] .= 0

    # Extract the correct dispersion object and spectral window from SVE
    disp = sve.dispersion
    swin = disp.spectral_window

    # Get the sigma of this ISRF in units of the spectral window
    σ = FWHM_to_sigma(ISRF.FWHM) * ISRF.FWHM_unit |> swin.ww_unit |> ustrip

    # Dispatch to lowlevel function to deal with type instability better
    return _calculate_dispersion_polynomial_jacobian_gauss!(
        inst_buf.low_res_output,
        sve.coefficient_order,
        σ,
        swin.ww_unit,
        swin.ww_grid,
        disp.ww_unit,
        disp.ww,
        disp.index,
        data, # this should be called radiance or jacobian really
        doppler_factor,
        extend
    )


end



function _calculate_dispersion_polynomial_jacobian_gauss!(
    output::Vector,
    order,
    σ,
    hires_unit,
    hires_ww,
    lores_unit,
    lores_ww,
    lores_index,
    data,
    doppler_factor,
    extend
    )

    @assert length(data) == length(hires_ww) (
        "Hi-res vector and hires vector grid must be same size!"
    )

    # Zero-out the output
    output[:] .= 0

    # buffer for ILS integration
    buffer = σ * extend
    # Unit conversion factor
    unit_fac = 1.0 * hires_unit / lores_unit
    unit_fac_stripped = ustrip(unit_fac)

    # Depending on λ/ν, we use a different effective Doppler formula
    doppler_term = 1.0
    if lores_unit isa Unitful.LengthUnits
        doppler_term = 1.0 / (1.0 + doppler_factor)
    elseif lores_unit isa Unitful.WavenumberUnits
        doppler_term = 1.0 + doppler_factor
    end

    for i in eachindex(lores_ww)

        if i == 1
            continue
        end

        this_l1b_idx = lores_index[i - 1]
        this_ww = unit_fac_stripped * lores_ww[i] / doppler_term

        # Need this factor to calculate the Jacobian for this SVE
        tmp_val = this_l1b_idx ^ order

        idx_left = Int(searchsortedfirst(hires_ww, (this_ww - buffer)))
        idx_right = Int(searchsortedlast(hires_ww, (this_ww + buffer)))

        if idx_left < 1
            continue
        end
        if idx_right > length(hires_ww)
            continue
        end

        running_sum = 0.0
        running_sum_isrf = 0.0

        for ii in idx_left:idx_right - 1

            ww_left = hires_ww[ii]
            ww_right = hires_ww[ii + 1]

            # Calculate dISRF / dwavelength - we can do this analytically
            ISRF_val_left = 1 / (σ * sqrt(2*pi)) * exp(-0.5 * (ww_left - this_ww)^2 / σ^2)
            ∂ISRF_val_left = ISRF_val_left * (ww_left - this_ww) / σ^2

            ISRF_val_right = 1 / (σ * sqrt(2*pi)) * exp(-0.5 * (ww_right - this_ww)^2 / σ^2)
            ∂ISRF_val_right = ISRF_val_right * (ww_right - this_ww) / σ^2

            running_sum += 0.5 * (
                ∂ISRF_val_left * data[ii] +
                ∂ISRF_val_right * data[ii + 1]
            ) * (ww_right - ww_left)

            running_sum_isrf += 0.5 * (
                ISRF_val_left +
                ISRF_val_right
            ) * (ww_right - ww_left)

        end

        output[this_l1b_idx + 1] = tmp_val * running_sum / running_sum_isrf

    end

    return true

end