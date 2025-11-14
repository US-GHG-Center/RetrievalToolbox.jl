"""

$(TYPEDSIGNATURES)
Applies the `TableISRF`-type instrument response function on some vector `data`, usually
intended to be the high-resolution model radiance or Jacobian. The output is stored in
`inst_buf.low_res_output`. Returns `true` if the calculation did not encounter any errors.

## Details
Application of the ISRF is in accordance with standard literature. The radiance ``I[i]``
at instrument-level on a particular spectral sample `i` at wavelength ``\\lambda_i`` is
given by

``I[i] = \\int I_\\mathrm{model}(\\lambda_i - \\Delta\\lambda) \\cdot
\\mathrm{ISRF}_i(\\Delta\\lambda) \\; d\\lambda``

where ``\\mathrm{ISRF}_i`` is the pre-calculated ISRF evaluated at the relative wavelength
``\\Delta\\lambda``. This implementation replaces the integral by a fast numerical
integration via the trapezoid rule.

The output is stored in `inst_buf.low_res_output`.

## Instrument Doppler shift

This function can take the `doppler_factor` argument, which is needed to compute the
appropriate Doppler shift caused by the relative motion between, for example, the
measurement footprint on the ground (for downlooking instruments) and the instrument along
the line of sight. The Doppler factor is ``v / c`` where ``c`` is the speed of light. This
function does not derive the appropriate Doppler shift, it is up to users to perform the
appropriate calculation.
"""
function apply_isrf_to_spectrum!(
    inst_buf::InstrumentBuffer,
    ISRF::TableISRF,
    disp::AbstractDispersion,
    data::AbstractVector,
    spectral_window::AbstractSpectralWindow;
    doppler_factor=0.0
    )

    # In Julia, we must encapsulate certain calculations in functions for explicit/clear
    # type inference. So all the heavy lifting of the "convolution" function is done
    # inside this following call.

    # This is also the place where we must account for units
    success = _apply_tableisrf_to_spectrum_lowlevel!(
        inst_buf.low_res_output,
        ISRF.ww_delta_unit,
        ISRF.ww_delta,
        ISRF.relative_response,
        spectral_window.ww_unit,
        spectral_window.ww_grid,
        data,
        disp.ww_unit,
        disp.ww,
        disp.index,
        doppler_factor,
        inst_buf.tmp1,
        inst_buf.tmp2,
    )

    return success

end

"""
$(TYPEDSIGNATURES)

Low-level function to performantly apply an instrument ISRF to a high-resolution model
spectrum, in order to obtain a spectrum comparable to that measured by an instrument.

"""
function _apply_tableisrf_to_spectrum_lowlevel!(
    lores_data, # Vector - Output!
    ISRF_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits}, # Length unit of the ISRF Δλ variable
    ISRF_ww_delta, # Array [delta ww, L1B index]
    ISRF_relative_response, # Array [delta w, L1B index],
    hires_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
    hires_ww, # Vector, hi-res wavelength/wavenumber
    hires_data, # Vector, hi-res spectrum
    lores_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
    lores_ww, # Vector, Dispersion wavelength/wavenumber
    lores_idx, # Vector, Dispersion index
    doppler_factor,
    tmp1,
    tmp2,
    )

    # Depending on λ/ν, we use a different effective Doppler formula
    if hires_unit isa Unitful.LengthUnits
        doppler_term = 1.0 / (1.0 + doppler_factor)
    elseif hires_unit isa Unitful.WavenumberUnits
        doppler_term = 1.0 + doppler_factor
    end


    # Fill with NaN's, so it is immediately clear which
    # indices have not been used.
    @views lores_data[:] .= NaN

    # We do not fill the buffer completely (usually),
    # so let's keep track of how many wavelengths we must process.
    N_low = length(lores_ww)

    # In the next steps, we mix wavelengths from the dispersion object
    # and the ISRF object, so let's make sure we get the units agreed

    unit_fac = 1.0 * ISRF_unit / lores_unit

    @views tmp1[:] .= 0.0
    @views tmp1[1:N_low] = lores_ww[:] .+ unit_fac * ISRF_ww_delta[1, lores_idx]
    idx_first_all = searchsortedfirst.(Ref(hires_ww), tmp1)

    @views tmp1[:] .= 0.0
    @views tmp1[1:N_low] = lores_ww[:] .+ unit_fac * ISRF_ww_delta[end, lores_idx]
    idx_last_all = searchsortedfirst.(Ref(hires_ww), tmp1)

    # Preallocate placeholder array for less in-loop allocations
    all_sizes = @. idx_last_all - idx_first_all + 1

    @views tmp2[:] .= 0.0
    this_ww_delta = @view tmp2[1:size(ISRF_ww_delta, 1)]

    idx_max = length(hires_ww)
    idx_min = 1

    for i in eachindex(lores_idx)

        this_l1b_idx = lores_idx[i]
        this_ww = lores_ww[i] / doppler_term

        idx_first = idx_first_all[i]
        idx_last = idx_last_all[i]

        if (idx_first >= idx_max) || (idx_last >= idx_max)
            @error "[INST] Sorry - idx_first larger than idx_max or last larger then max"
            @error "[INST] Sample $(i) at $(this_ww), idx_first: $(idx_first), idx_max: $(idx_max)"
            @error "[INST] Sample $(i) at $(this_ww), idx_last: $(idx_last), idx_max: $(idx_max)"
            @error "[INST] Try increasing the buffer of the spectral window!"
            return false
        end

        @views this_ww_delta[:] = unit_fac * ISRF_ww_delta[:, this_l1b_idx] .+ this_ww

        # Create a view to the placeholder array in which we stick the
        # interpolated ISRF.
        this_ISRF = @view tmp1[1:all_sizes[i]]
        @views this_ISRF[:] .= 0.0

        this_relative_response = @view ISRF_relative_response[:, this_l1b_idx]
        this_hires_ww = @view hires_ww[idx_first:idx_last]

        # Interpolate and store result in `this_ISRF`
        pwl_value_1d!(
            this_ww_delta,
            this_relative_response,
            this_hires_ww,
            this_ISRF
        )

        # Calculate the final pixel value by multiplying in-place
        # and then do the integration via trapezoidal rule.

        this_ISRF_sum = _trapz(this_hires_ww, this_ISRF)

        # Multiply in the high-resolution spectrum
        @views this_ISRF[:] .*= hires_data[idx_first:idx_last]
        lores_data[lores_idx[i]] = _trapz(this_hires_ww, this_ISRF) / this_ISRF_sum

    end

    # Successful run!
    return true

end



"""
$(TYPEDSIGNATURES)

Applies the `GaussISRF`-type instrument response function on some vector `data`, usually
intended to be the high-resolution model radiance or Jacobian. Integration limits are
given by the `extend` keyword, it defines how many ``\\sigma`` (derived fom `ISRF.FWHM`)
the integration limits are on both sides. The output is stored in
`inst_buf.low_res_output`. Returns `true` if the calculation did not encounter any errors.

## Details

Application of the ISRF is in accordance with standard literature. The radiance ``I[i]``
at instrument-level on a particular spectral sample `i` at wavelength ``\\lambda_i`` is
given by

``I[i] = \\int I_\\mathrm{model}(\\lambda_i - \\Delta\\lambda) \\cdot
\\mathrm{ISRF}(\\Delta\\lambda) \\; d\\lambda``

where ``\\mathrm{ISRF}`` is the Gaussian ISRF evaluated at the relative wavelength
``\\Delta\\lambda``. This implementation replaces the integral by a fast numerical
integration.

## Instrument Doppler shift

This function can take the `doppler_factor` argument, which is needed to compute the
appropriate Doppler shift caused by the relative motion between, for example, the
measurement footprint on the ground (for downlooking instruments) and the instrument along
the line of sight. The Doppler factor is ``v / c`` where ``c`` is the speed of light. This
function does not derive the appropriate Doppler shift, it is up to users to perform the
appropriate calculation.
"""
function apply_isrf_to_spectrum!(
    inst_buf::InstrumentBuffer,
    ISRF::GaussISRF,
    disp::AbstractDispersion,
    data::AbstractVector,
    swin::AbstractSpectralWindow;
    doppler_factor=0.0,
    extend=4.0 # How many σ's to move away from the center for integration?
    )

    # Convert to σ here, so we can pass a number to the lowlevel function
    # (if you include this code in the function itself, it produces a lot of
    #  additional allocations.)
    σ = FWHM_to_sigma(ISRF.FWHM) * ISRF.FWHM_unit |> swin.ww_unit |> ustrip

    success = _apply_gaussisrf_to_spectrum_lowlevel!(
        inst_buf.low_res_output,
        σ,
        extend,
        swin.ww_unit,
        swin.ww_grid,
        data,
        disp.ww_unit,
        disp.ww,
        disp.index,
        doppler_factor,
    )

    return success

end

function _apply_gaussisrf_to_spectrum_lowlevel!(
    lores_data, # Vector - Output!
    σ::Number,
    extend::Number,
    hires_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
    hires_ww, # Vector, hi-res wavelength/wavenumber
    hires_data, # Vector, hi-res spectrum
    lores_unit::Union{Unitful.LengthUnits, Unitful.WavenumberUnits},
    lores_ww, # Vector, Dispersion wavelength/wavenumber
    lores_idx, # Vector, Dispersion index
    doppler_factor,
)

    lores_data[:] .= 0

    # Depending on λ/ν, we use a different effective Doppler formula
    if hires_unit isa Unitful.LengthUnits
        doppler_term = 1.0 / (1.0 + doppler_factor)
    elseif hires_unit isa Unitful.WavenumberUnits
        doppler_term = 1.0 + doppler_factor
    end

    unit_fac = 1.0 * hires_unit / lores_unit

    for i in eachindex(lores_idx)

        this_ww = unit_fac * lores_ww[i] / doppler_term # Create spectral value, include Doppler

        idx_left = searchsortedfirst(hires_ww, this_ww - extend * σ)
        idx_right = searchsortedlast(hires_ww, this_ww + extend * σ) + 1

        running_sum = 0
        running_sum_isrf = 0

        for ii in idx_left:idx_right-1

            ww_left = hires_ww[ii]
            ww_right = hires_ww[ii + 1]

            # Implementation of the trapezoidal rule for this integration
            isrf_left = 1 / (σ * sqrt(2*pi)) *
                exp(-0.5 * (ww_left - lores_ww[i])^2 / σ^2)
            isrf_right = 1 / (σ * sqrt(2*pi)) *
                exp(-0.5 * (ww_right - lores_ww[i])^2 / σ^2)

            running_sum += 0.5 * (
                hires_data[ii] * isrf_left + hires_data[ii+1] * isrf_right
            ) * (ww_right - ww_left)

            running_sum_isrf +=  0.5 * (
                isrf_left + isrf_right
            ) * (ww_right - ww_left)

        end

        lores_data[lores_idx[i]] = running_sum / running_sum_isrf

    end

    # Success
    return true

end