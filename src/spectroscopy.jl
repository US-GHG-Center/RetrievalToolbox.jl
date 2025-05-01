"""
Loads a complete ABSCO table into memory and creates an
ABSCOSpectroscopy3D/4D object, depending on whether the
ABSCO file contains H2O broadening information.

$(SIGNATURES)

# Details

Wavenumber dimension and the corresponding axis for the
coefficient array are flipped to give wavelengths in
increasing order. User can supply a scale factor which
is multiplied into the entire table.
"""
function load_ABSCO_spectroscopy(
    fname::String;
    spectral_unit=:Wavelength,
    scale_factor=1.0f0
    )

    if !isfile(fname)
        @error "ABSCO file location: $(fname) not a valid file!"
        return nothing
    end

    h5 = h5open(fname, "r")

    # We hard-code the units here. ABSCO tables have so far
    # only come in units of wavenumber (1/cm). Let's just
    # turm them into microns. Same goes for the temparature
    # and pressure axes.
    # Maybe in the future they'll be in netCDF format..

    if spectral_unit == :Wavelength
        ww_unit = u"μm"
        ww = 1.0f4 ./ h5["Wavenumber"][:]
    elseif spectral_unit == :Wavenumber
        ww_unit = u"cm^-1"
        ww = h5["Wavenumber"][:]
    else
        @error "Spectral unit must be either :Wavenumber or :Wavelength"
    end

    temperatures_unit = u"K"
    temperatures = h5["Temperature"][:,:]

    pressures_unit = u"Pa"
    pressures = h5["Pressure"][:]


    # ABSCO tables come in either 3 (T, p, wavenumber) or
    # 4 (T, p, wavenumber, H2O VMR) dimensions. Find out
    # which one we have here.

    if haskey(h5, "Broadener_Index")
        has_H2O = true
        broadener_index = h5["Broadener_Index"][1]
        broadener_key = "Broadener_$(broadener_index)_VMR"
        broadener_vmrs = h5[broadener_key][:]

        @debug "This ABSCO table has a broadener: " * broadener_index
    else
        has_H2O = false
    end


    # We have one gas per ABSCO table file. Find
    # out which gas we are using here..

    absorption_key = ""
    for key in keys(h5)
        if occursin("Absorption", key)
            absorption_key = key
        end
    end

    if haskey(h5[absorption_key], "gas_name")
        gas_name = read(HDF5.attributes(h5[absorption_key])["gas_name"])
    else
        gas_name = ""
    end

    # ABSCO table dimension ordering:
    #  - Wavelength
    #  - Broadener VMR (optional)
    #  - Temperature
    #  - Pressure


    # We should use the same number type here, so we must convert:
    scale_factor = convert(eltype(h5[absorption_key]), scale_factor)
    cross_section_unit = u"cm^2" # actually, cm^2 per molecule, but we omit the molecule
    cross_section = read(h5[absorption_key]) .* scale_factor

    close(h5)

    # If we want to store in wavelength order, we need to reverse the cross_section
    if spectral_unit == :Wavelength
        reverse!(cross_section, dims=1)
        reverse!(ww)
    end

    if has_H2O

        return ABSCOSpectroscopy4D(
            fname,
            gas_name,
            scale_factor,
            ww,
            ww_unit,
            temperatures,
            temperatures_unit,
            pressures,
            pressures_unit,
            broadener_vmrs,
            cross_section,
            cross_section_unit
        )

    else

        return ABSCOSpectroscopy3D(
            fname,
            gas_name,
            scale_factor,
            ww,
            ww_unit,
            temperatures,
            temperatures_unit,
            pressures,
            pressures_unit,
            cross_section,
            cross_section_unit
        )

    end

end


"""
Given two wavelength arrays wl1 and wl2, this function computes
the indices that sort wl2 into wl1. Similar to a broadcast
searchsorted (think: searchsorted.(Ref(wl1), wl2)), however
faster if both wl1 and wl2 are sorted. This is used to e.g.
find the indices of the spectral window wavelength relative
to spectroscopy wavelengths.

$(TYPEDSIGNATURES)
"""
function find_wavelength_indices(wl1, wl2)

    output_idx = -ones(Int, length(wl2))
    first_idx = 1
    last_wl_idx = 1

    @inbounds for i in eachindex(wl2)
        if wl2[i] > wl1[1]
            first_idx = i
            break
        end
    end


    @inbounds for i in first_idx:length(wl2)

        if wl2[i] < wl1[1]
            output_idx[i] = -1
            continue
        end

        if wl2[i] > wl1[end]
            output_idx[i] = -1
            continue
        end

        for j in last_wl_idx:(length(wl1) - 1)
            if (wl2[i] >= wl1[j]) & (wl2[i] <= wl1[j + 1])

                output_idx[i] = j
                last_wl_idx = j
                break

            end
        end

    end

    return output_idx

end


"""
A wrapper function which retrieves the absorption cross_section
for a vector of wavelengths, and corresponding (scalar) values
of temperature, pressure, H2O broadening VMR. It reads the
needed data from a 4D ABSCO object.

$(SIGNATURES)

# Details

This calls the underlying, more explicit function where the
spectroscopy arguments are split up. In Julia, making those
arguments explicit provides a very large performance gain,
since the compiler can infer the types.
"""
function get_absorption_coefficient_value_at!(
    output::AbstractVector,
    spec::ABSCOSpectroscopy4D,
    wl::AbstractVector,
    absco_wl_idx_left::Vector{Int},
    is_matched_wl::Bool,
    p::Number,
    T::Number,
    H2O::Number
    )

    _get_absorption_coefficient_value_at!(
        output,
        spec.ww,
        spec.pressures,
        spec.temperatures,
        spec.broadener_vmrs,
        spec.cross_section,
        wl,
        absco_wl_idx_left,
        is_matched_wl,
        p,
        T,
        H2O
    )

end


function _get_absorption_coefficient_value_at!(
    output::AbstractVector,
    spec_ww::AbstractVector,
    spec_pressures::AbstractVector,
    spec_temperatures::AbstractArray,
    spec_H2O::AbstractVector,
    spec_cross_section::AbstractArray,
    wl::AbstractVector,
    absco_wl_idx_left::Vector{Int},
    is_matched_wl::Bool,
    p,
    T,
    H2O
)

    idx_r_p = 1
    idx_l_p = 1

    idx_lr_T = 1
    idx_ll_T = 1

    idx_rl_T = 1
    idx_rr_T = 1

    idx_l_H2O = 1
    idx_r_H2O = 1

    p_d = zero(eltype(spec_pressures))
    T_d_l = zero(eltype(spec_temperatures))
    T_d_r = zero(eltype(spec_temperatures))
    H2O_d = zero(eltype(spec_H2O))

    if p < spec_pressures[1]
        idx_l_p = 1
    elseif p > spec_pressures[end]
        idx_l_p = size(spec_pressures, 1) - 1
    else
        for j in 1:length(spec_pressures)-1
            @inbounds if (p >= spec_pressures[j]) & (p <= spec_pressures[j+1])
                idx_l_p = j
                break
            end
        end
    end
    idx_r_p = idx_l_p + 1

    if T < spec_temperatures[1, idx_l_p]
        idx_ll_T = 1
    elseif T > spec_temperatures[end, idx_l_p]
        idx_ll_T = size(spec_temperatures, 1) - 1
    else
        for j in 1:size(spec_temperatures, 1) - 1
            @inbounds if (T >= spec_temperatures[j, idx_l_p]) &
                (T <= spec_temperatures[j+1, idx_l_p])
                idx_ll_T = j
                break
            end
        end
    end
    idx_lr_T = idx_ll_T + 1

    if T < spec_temperatures[1, idx_r_p]
        idx_rl_T = 1
    elseif T > spec_temperatures[end, idx_r_p]
        idx_rl_T = size(spec_temperatures, 1) - 1
    else
        for j in 1:size(spec_temperatures, 1)-1
            @inbounds if (T >= spec_temperatures[j, idx_r_p]) &
                (T <= spec_temperatures[j+1, idx_r_p])
                idx_rl_T = j
                break
            end
        end
    end
    idx_rr_T = idx_rl_T + 1


    if H2O < spec_H2O[1]
        idx_l_H2O = 1
    elseif H2O > spec_H2O[end]
        idx_l_H2O = size(spec_H2O, 1) - 1
    else
        @inbounds for j in 1:length(spec_H2O)
            if (H2O >= spec_H2O[j]) & (H2O <= spec_H2O[j+1])
                idx_l_H2O = j
                break
            end
        end
    end
    idx_r_H2O = idx_l_H2O + 1


    H2O_d = (H2O - spec_H2O[idx_l_H2O]) /
        (spec_H2O[idx_r_H2O] - spec_H2O[idx_l_H2O])
    p_d = (p - spec_pressures[idx_l_p]) /
        (spec_pressures[idx_r_p] - spec_pressures[idx_l_p])
    T_d_l = (T - spec_temperatures[idx_ll_T, idx_l_p]) /
        (spec_temperatures[idx_lr_T, idx_l_p] - spec_temperatures[idx_ll_T, idx_l_p])
    T_d_r = (T - spec_temperatures[idx_rl_T, idx_r_p]) /
        (spec_temperatures[idx_rr_T, idx_r_p] - spec_temperatures[idx_rl_T, idx_r_p])


    one_minus_H2O_d = one(typeof(H2O_d)) - H2O_d
    one_minus_p_d = one(typeof(p_d)) - p_d
    one_minus_T_d_l = one(typeof(T_d_l)) - T_d_l
    one_minus_T_d_r = one(typeof(T_d_r)) - T_d_r


    if !is_matched_wl

        @turbo for wl_idx in eachindex(wl)

            idx_l_wl = absco_wl_idx_left[wl_idx]
            idx_r_wl = idx_l_wl + 1
            wl_d = (wl[wl_idx] - spec_ww[idx_l_wl]) /
                (spec_ww[idx_r_wl] -  spec_ww[idx_l_wl])
            one_minus_wl_d = (1.0 - wl_d)

            C3_111 = spec_cross_section[idx_l_wl, idx_l_H2O, idx_ll_T, idx_l_p]
            C3_211 = spec_cross_section[idx_l_wl, idx_r_H2O, idx_ll_T, idx_l_p]
            C3_121 = spec_cross_section[idx_l_wl, idx_l_H2O, idx_lr_T, idx_l_p]
            C3_112 = spec_cross_section[idx_l_wl, idx_l_H2O, idx_rl_T, idx_r_p]
            C3_221 = spec_cross_section[idx_l_wl, idx_r_H2O, idx_lr_T, idx_l_p]
            C3_122 = spec_cross_section[idx_l_wl, idx_l_H2O, idx_rr_T, idx_r_p]
            C3_212 = spec_cross_section[idx_l_wl, idx_r_H2O, idx_rl_T, idx_r_p]
            C3_222 = spec_cross_section[idx_l_wl, idx_r_H2O, idx_rr_T, idx_r_p]

            C3_111 *= one_minus_wl_d
            C3_211 *= one_minus_wl_d
            C3_121 *= one_minus_wl_d
            C3_112 *= one_minus_wl_d
            C3_221 *= one_minus_wl_d
            C3_122 *= one_minus_wl_d
            C3_212 *= one_minus_wl_d
            C3_222 *= one_minus_wl_d

            # C3 = (h2o, temp, pressure)
            C3_111 += spec_cross_section[idx_r_wl, idx_l_H2O, idx_ll_T, idx_l_p] * wl_d
            C3_211 += spec_cross_section[idx_r_wl, idx_r_H2O, idx_ll_T, idx_l_p] * wl_d
            C3_121 += spec_cross_section[idx_r_wl, idx_l_H2O, idx_lr_T, idx_l_p] * wl_d
            C3_112 += spec_cross_section[idx_r_wl, idx_l_H2O, idx_rl_T, idx_r_p] * wl_d
            C3_221 += spec_cross_section[idx_r_wl, idx_r_H2O, idx_lr_T, idx_l_p] * wl_d
            C3_122 += spec_cross_section[idx_r_wl, idx_l_H2O, idx_rr_T, idx_r_p] * wl_d
            C3_212 += spec_cross_section[idx_r_wl, idx_r_H2O, idx_rl_T, idx_r_p] * wl_d
            C3_222 += spec_cross_section[idx_r_wl, idx_r_H2O, idx_rr_T, idx_r_p] * wl_d

            C2_11 = C3_111 * one_minus_H2O_d + C3_211 * H2O_d
            C2_12 = C3_112 * one_minus_H2O_d + C3_212 * H2O_d
            C2_21 = C3_121 * one_minus_H2O_d + C3_221 * H2O_d
            C2_22 = C3_122 * one_minus_H2O_d + C3_222 * H2O_d

            C1_1 = C2_11 * one_minus_T_d_l + C2_21 * T_d_l
            C1_2 = C2_12 * one_minus_T_d_r + C2_22 * T_d_r

            output[wl_idx] = one_minus_p_d * C1_1 + p_d * C1_2



        end

        # Set all out-of-ABSCO-range points to zero, and over-write
        # whatever result was produced above.
        @inbounds for wl_idx in eachindex(wl)
            if absco_wl_idx_left[wl_idx] == -1
                output[wl_idx] = 0
            end
        end

    else

        @turbo for wl_idx in eachindex(wl)

            idx_l_wl = absco_wl_idx_left[wl_idx]

            C3_111 = spec_cross_section[idx_l_wl, idx_l_H2O, idx_ll_T, idx_l_p]
            C3_211 = spec_cross_section[idx_l_wl, idx_r_H2O, idx_ll_T, idx_l_p]
            C3_121 = spec_cross_section[idx_l_wl, idx_l_H2O, idx_lr_T, idx_l_p]
            C3_112 = spec_cross_section[idx_l_wl, idx_l_H2O, idx_rl_T, idx_r_p]
            C3_221 = spec_cross_section[idx_l_wl, idx_r_H2O, idx_lr_T, idx_l_p]
            C3_122 = spec_cross_section[idx_l_wl, idx_l_H2O, idx_rr_T, idx_r_p]
            C3_212 = spec_cross_section[idx_l_wl, idx_r_H2O, idx_rl_T, idx_r_p]
            C3_222 = spec_cross_section[idx_l_wl, idx_r_H2O, idx_rr_T, idx_r_p]

            C2_11 = C3_111 * one_minus_H2O_d + C3_211 * H2O_d
            C2_12 = C3_112 * one_minus_H2O_d + C3_212 * H2O_d
            C2_21 = C3_121 * one_minus_H2O_d + C3_221 * H2O_d
            C2_22 = C3_122 * one_minus_H2O_d + C3_222 * H2O_d

            C1_1 = C2_11 * one_minus_T_d_l + C2_21 * T_d_l
            C1_2 = C2_12 * one_minus_T_d_r + C2_22 * T_d_r

            output[wl_idx] = one_minus_p_d * C1_1 + p_d * C1_2


        end

    end

end



"""
A wrapper function which retrieves the absorption cross_section
for a vector of wavelengths, and corresponding (scalar) values
of temperature, pressure, H2O broadening VMR. It reads the
needed data from a 4D ABSCO object.

$(SIGNATURES)

# Details

This calls the underlying, more explicit function where the
spectroscopy arguments are split up. In Julia, making those
arguments explicit provides a very large performance gain,
since the compiler can infer the types.
"""
function get_absorption_coefficient_value_at!(
    output::AbstractVector,
    spec::ABSCOSpectroscopy3D,
    wl::AbstractVector,
    absco_wl_idx_left::Vector{Int},
    is_matched_wl::Bool,
    p::Number,
    T::Number,
    H2O::Number,
    )

    _get_absorption_coefficient_value_at!(
        output,
        spec.ww,
        spec.pressures,
        spec.temperatures,
        spec.cross_section,
        wl,
        absco_wl_idx_left,
        is_matched_wl,
        p,
        T
    )

end

"""
Low-level function that calculates absorption cross_section

$(TYPEDSIGNATURES)
"""
function _get_absorption_coefficient_value_at!(
    output::AbstractVector,
    spec_ww::AbstractVector,
    spec_pressures::AbstractVector,
    spec_temperatures::AbstractArray,
    spec_cross_section::AbstractArray,
    wl::AbstractVector,
    absco_wl_idx_left::Vector{Int},
    is_matched_wl::Bool,
    p,
    T
    )

    idx_r_p = 1
    idx_l_p = 1

    idx_lr_T = 1
    idx_ll_T = 1

    idx_rl_T = 1
    idx_rr_T = 1

    p_d = zero(eltype(spec_pressures))
    T_d_l = zero(eltype(spec_temperatures))
    T_d_r = zero(eltype(spec_temperatures))

    if p < spec_pressures[1]
        idx_l_p = 1
    elseif p > spec_pressures[end]
        idx_l_p = size(spec_pressures, 1) - 1
    else
        for j in 1:length(spec_pressures)-1
            @inbounds if (p >= spec_pressures[j]) & (p <= spec_pressures[j+1])
                idx_l_p = j
                break
            end
        end
    end
    idx_r_p = idx_l_p + 1

    if T < spec_temperatures[1, idx_l_p]
        idx_ll_T = 1
    elseif T > spec_temperatures[end, idx_l_p]
        idx_ll_T = size(spec_temperatures, 1) - 1
    else
        for j in 1:size(spec_temperatures, 1) - 1
            @inbounds if (T >= spec_temperatures[j, idx_l_p]) &
                (T <= spec_temperatures[j+1, idx_l_p])
                idx_ll_T = j
                break
            end
        end
    end
    idx_lr_T = idx_ll_T + 1

    if T < spec_temperatures[1, idx_r_p]
        idx_rl_T = 1
    elseif T > spec_temperatures[end, idx_r_p]
        idx_rl_T = size(spec_temperatures, 1) - 1
    else
        for j in 1:size(spec_temperatures, 1)-1
            @inbounds if (T >= spec_temperatures[j, idx_r_p]) &
                (T <= spec_temperatures[j+1, idx_r_p])
                idx_rl_T = j
                break
            end
        end
    end
    idx_rr_T = idx_rl_T + 1


    p_d = (p - spec_pressures[idx_l_p]) / (spec_pressures[idx_r_p] - spec_pressures[idx_l_p])
    T_d_l = (T - spec_temperatures[idx_ll_T, idx_l_p]) /
        (spec_temperatures[idx_lr_T, idx_l_p] - spec_temperatures[idx_ll_T, idx_l_p])
    T_d_r = (T - spec_temperatures[idx_rl_T, idx_r_p]) /
        (spec_temperatures[idx_rr_T, idx_r_p] - spec_temperatures[idx_rl_T, idx_r_p])

    one_minus_p_d = one(typeof(p_d)) - p_d
    one_minus_T_d_l = one(typeof(T_d_l)) - T_d_l
    one_minus_T_d_r = one(typeof(T_d_r)) - T_d_r

    if !is_matched_wl

        @turbo for wl_idx in eachindex(wl)

            idx_l_wl = absco_wl_idx_left[wl_idx]
            idx_r_wl = idx_l_wl + 1

            wl_d = (wl[wl_idx] - spec_ww[idx_l_wl]) /
                (spec_ww[idx_r_wl] -  spec_ww[idx_l_wl])
            one_minus_wl_d = (1.0 - wl_d)


            C2_11 = spec_cross_section[idx_l_wl, idx_ll_T, idx_l_p]
            C2_21 = spec_cross_section[idx_l_wl, idx_lr_T, idx_l_p]
            C2_12 = spec_cross_section[idx_l_wl, idx_rl_T, idx_r_p]
            C2_22 = spec_cross_section[idx_l_wl, idx_rr_T, idx_r_p]

            C2_11 *= one_minus_wl_d
            C2_21 *= one_minus_wl_d
            C2_12 *= one_minus_wl_d
            C2_22 *= one_minus_wl_d

            # C3 = (h2o, temp, pressure)
            C2_11 += spec_cross_section[idx_r_wl, idx_ll_T, idx_l_p] * wl_d
            C2_21 += spec_cross_section[idx_r_wl, idx_lr_T, idx_l_p] * wl_d
            C2_12 += spec_cross_section[idx_r_wl, idx_rl_T, idx_r_p] * wl_d
            C2_22 += spec_cross_section[idx_r_wl, idx_rr_T, idx_r_p] * wl_d

            C1_1 = C2_11 * one_minus_T_d_l + C2_21 * T_d_l
            C1_2 = C2_12 * one_minus_T_d_r + C2_22 * T_d_r

            output[wl_idx] = one_minus_p_d * C1_1 + p_d * C1_2

        end

        # Set all out-of-ABSCO-range points to zero, and over-write
        # whatever result was produced above.
        @inbounds for wl_idx in eachindex(wl)
            if absco_wl_idx_left[wl_idx] == -1
                output[wl_idx] = 0
            end
        end

    else

        @turbo for wl_idx in eachindex(wl)

            idx_l_wl = absco_wl_idx_left[wl_idx]

            C2_11 = spec_cross_section[idx_l_wl, idx_ll_T, idx_l_p]
            C2_21 = spec_cross_section[idx_l_wl, idx_lr_T, idx_l_p]
            C2_12 = spec_cross_section[idx_l_wl, idx_rl_T, idx_r_p]
            C2_22 = spec_cross_section[idx_l_wl, idx_rr_T, idx_r_p]

            C1_1 = C2_11 * one_minus_T_d_l + C2_21 * T_d_l
            C1_2 = C2_12 * one_minus_T_d_r + C2_22 * T_d_r

            output[wl_idx] = one_minus_p_d * C1_1 + p_d * C1_2


        end

    end

end
