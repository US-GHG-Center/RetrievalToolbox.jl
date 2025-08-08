#=
Notes regarding the handling of units:

=#



"""
$(TYPEDSIGNATURES)

Loads a complete ABSCO table into memory and creates an ABSCOSpectroscopy3D/4D object,
depending on whether the ABSCO file contains H2O broadening information.

# Details

Wavenumber dimension and the corresponding axis for the coefficient array are flipped to
give wavelengths in increasing order. User can supply a scale factor which is multiplied
into the entire table.
"""
function load_ABSCO_spectroscopy(
    fname::String;
    spectral_unit=:Wavelength,
    scale_factor=1.0f0,
    distributed=false
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

        @debug "[SPEC] This ABSCO table has a broadener: " * broadener_index
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

    # We handle things a little differently if we use SharedArrays for distributed
    # computing.

    if distributed
        # Create a SharedArray
        cross_section = SharedArray{eltype(h5[absorption_key])}(
            size(h5[absorption_key])...)
        # .. fill with values (we know this is 4D)
        cross_section[:,:,:,:] = h5[absorption_key][:,:,:,:]
    else
        # "normal" array
        cross_section = read(h5[absorption_key]) .* scale_factor
    end

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
$(TYPEDSIGNATURES)

Loads a complete ABSCO table into memory and creates an ABSCOSpectroscopy3D/4D object,
depending on whether the ABSCO file contains H2O broadening information.

# Details

Wavenumber dimension and the corresponding axis for the coefficient array are flipped to
give wavelengths in increasing order. User can supply a scale factor which is multiplied
into the entire table.
"""
function load_ABSCOAER_spectroscopy(
    fname::String;
    spectral_unit=:Wavelength,
    scale_factor=1.0,
    distributed=false
    )

    if !isfile(fname)
        @error "ABSCO file location: $(fname) not a valid file!"
        return nothing
    end

    nc = NCDataset(fname, "r")

    # Parse the gas name from the nc attribute
    gas_name = split(lstrip(split(nc.attrib["description"], "for")[2]))[1] |> string

    # Read the spectral unit inside the file
    absco_spec_unit = nc["Spectral_Grid"].attrib["units"]

    if (spectral_unit == :Wavelength) && (absco_spec_unit == "cm-1")
        ww_unit = u"μm"
        ww = 1e4 ./ nc["Spectral_Grid"].var[:]
    elseif (spectral_unit == :Wavenumber) && (absco_spec_unit == "cm-1")
        ww_unit = u"cm^-1"
        ww = nc["Spectral_Grid"].var[:]
    else
        @error "Spectral unit must be either :Wavenumber or :Wavelength, and " *
        "the source units must be `cm-1`."
    end

    # IMPORTANT!
    # In the ABSCO files, the atmosphere is flipped compared to our convention, which
    # has the top of the atmosphere as the first element.

    temperatures_unit = uparse(nc["T_layer"].attrib["units"])
    temperatures_array = nc["T_layer"].var[:,end:-1:1]

    pressures_unit = uparse(nc["P_level"].attrib["units"])
    pressures_array = nc["P_level"].var[end:-1:1]

    if nc["Cross_Section"].attrib["units"] == "cm2/molecule"
        cross_section_unit = u"cm^2"
    else
        @error "Problem with cross section unit - I do not understand the supplied unit!"
    end

    # For now, we only understand broadening due to water vapor - but only if available!
    # Remember - we have to flip the cross section table also, so that it follows our
    # convention with the first element being the lowest-pressure one.

    if "H2O_VMR" in NCDatasets.listVar(nc.ncid)
        broadener_vmrs = nc["H2O_VMR"].var[:]
        cross_section = nc["Cross_Section"].var[:,end:-1:1,:,:]
        close(nc)

        return ABSCOAERSpectroscopy4D(
            fname,
            gas_name,
            scale_factor,
            ww,
            ww_unit,
            temperatures_array,
            temperatures_unit,
            pressures_array,
            pressures_unit,
            broadener_vmrs,
            cross_section,
            cross_section_unit
        )

    else

        cross_section = nc["Cross_Section"].var[end:-1:1,:,:]
        close(nc)

        return ABSCOAERSpectroscopy3D(
            fname,
            gas_name,
            scale_factor,
            ww,
            ww_unit,
            temperatures_array,
            temperatures_unit,
            pressures_array,
            pressures_unit,
            cross_section,
            cross_section_unit
        )

    end



end

"""
$(TYPEDSIGNATURES)

Given two wavelength or wavenumber arrays ww1 and ww2, this function computes the indices
that sort ww2 into ww1. This provides the same functionality to a broadcast searchsorted
(think: searchsorted.(Ref(ww1), ww2)), however faster if both ww1 and ww2 are sorted. This
is used to e.g. find the indices of the spectral window wavelength relative to
spectroscopy wavelengths. For example, for inputs ww1 = [100., 110., 120., 140., 145.],
ww2 = [105., 111., 142.], the function will return [1,2,4], signifying that ww2[1] can be
inserted between ww1[1] and ww1[2].

NOTE! This function does *not* check if ww1 and ww2 are sorted and will return meaningless
results if either ww1 or ww2 are not sorted!
"""
function _find_ww_indices(ww1, ww2)

    output_idx = -ones(Int, length(ww2))
    first_idx = 1
    last_ww_idx = 1

    @inbounds for i in eachindex(ww2)
        if ww2[i] > ww1[1]
            first_idx = i
            break
        end
    end


    @inbounds for i in first_idx:length(ww2)

        if ww2[i] < ww1[1]
            output_idx[i] = -1
            continue
        end

        if ww2[i] > ww1[end]
            output_idx[i] = -1
            continue
        end

        for j in last_ww_idx:(length(ww1) - 1)
            if (ww2[i] >= ww1[j]) & (ww2[i] <= ww1[j + 1])

                output_idx[i] = j
                last_ww_idx = j
                break

            end
        end

    end

    return output_idx

end


"""
$(TYPEDSIGNATURES)

A wrapper function which retrieves the absorption cross_section for a vector of
wavelengths, and corresponding (scalar) values of temperature, pressure, H2O broadening
VMR. It reads the needed data from a 4D ABSCO object.

# Details

This calls the underlying, more explicit function where the spectroscopy arguments are
split up. In Julia, making those arguments explicit provides a very large performance
gain, since the compiler can infer the types.
"""
function get_cross_section_value_at!(
    output::AbstractVector,
    spec::ABSCOSpectroscopy4D,
    wl::AbstractVector,
    absco_wl_idx_left::Vector{Int},
    is_matched_wl::Bool,
    p::Unitful.Pressure,
    T::Unitful.Temperature,
    H2O::Number
    )

    _get_absco_cross_section_value_at!(
        output,
        spec.ww,
        spec.pressures,
        spec.temperatures,
        spec.broadener_vmrs,
        spec.cross_section,
        wl,
        absco_wl_idx_left,
        is_matched_wl,
        p |> spec.pressures_unit |> ustrip,
        T |> spec.temperatures_unit |> ustrip,
        H2O
    )

end

"""
$(TYPEDSIGNATURES)

This function is a high-performance, explicit interpolation routine that returns the
absorption cross section values of an 4D ABSCO table for a given set of coordinates in
spectral, pressure, temperature and water vapor space. The spectral coordinates must be
a vector, since in most applications we want to obtain the cross sections for all relevant
wavelengths (or wavenumbers) at once.

This function differs from a simple 4D-interpolation due to the irregular grid of the
ABSCO tables, whose temperature axis depends on the the pressure coordinate.
"""
function _get_absco_cross_section_value_at!(
    output::AbstractVector,
    spec_ww::AbstractVector,
    spec_pressures::AbstractVector,
    spec_temperatures::AbstractArray,
    spec_H2O::AbstractVector,
    spec_cross_section::AbstractArray,
    ww::AbstractVector,
    absco_ww_idx_left::Vector{Int},
    is_matched_ww::Bool,
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


    if !is_matched_ww

        @turbo for ww_idx in eachindex(ww)

            idx_l_ww = absco_ww_idx_left[ww_idx]
            idx_r_ww = idx_l_ww + 1
            ww_d = (ww[ww_idx] - spec_ww[idx_l_ww]) /
                (spec_ww[idx_r_ww] -  spec_ww[idx_l_ww])
            one_minus_ww_d = (1.0 - ww_d)

            C3_111 = spec_cross_section[idx_l_ww, idx_l_H2O, idx_ll_T, idx_l_p]
            C3_211 = spec_cross_section[idx_l_ww, idx_r_H2O, idx_ll_T, idx_l_p]
            C3_121 = spec_cross_section[idx_l_ww, idx_l_H2O, idx_lr_T, idx_l_p]
            C3_112 = spec_cross_section[idx_l_ww, idx_l_H2O, idx_rl_T, idx_r_p]
            C3_221 = spec_cross_section[idx_l_ww, idx_r_H2O, idx_lr_T, idx_l_p]
            C3_122 = spec_cross_section[idx_l_ww, idx_l_H2O, idx_rr_T, idx_r_p]
            C3_212 = spec_cross_section[idx_l_ww, idx_r_H2O, idx_rl_T, idx_r_p]
            C3_222 = spec_cross_section[idx_l_ww, idx_r_H2O, idx_rr_T, idx_r_p]

            C3_111 *= one_minus_ww_d
            C3_211 *= one_minus_ww_d
            C3_121 *= one_minus_ww_d
            C3_112 *= one_minus_ww_d
            C3_221 *= one_minus_ww_d
            C3_122 *= one_minus_ww_d
            C3_212 *= one_minus_ww_d
            C3_222 *= one_minus_ww_d

            # C3 = (h2o, temp, pressure)
            C3_111 += spec_cross_section[idx_r_ww, idx_l_H2O, idx_ll_T, idx_l_p] * ww_d
            C3_211 += spec_cross_section[idx_r_ww, idx_r_H2O, idx_ll_T, idx_l_p] * ww_d
            C3_121 += spec_cross_section[idx_r_ww, idx_l_H2O, idx_lr_T, idx_l_p] * ww_d
            C3_112 += spec_cross_section[idx_r_ww, idx_l_H2O, idx_rl_T, idx_r_p] * ww_d
            C3_221 += spec_cross_section[idx_r_ww, idx_r_H2O, idx_lr_T, idx_l_p] * ww_d
            C3_122 += spec_cross_section[idx_r_ww, idx_l_H2O, idx_rr_T, idx_r_p] * ww_d
            C3_212 += spec_cross_section[idx_r_ww, idx_r_H2O, idx_rl_T, idx_r_p] * ww_d
            C3_222 += spec_cross_section[idx_r_ww, idx_r_H2O, idx_rr_T, idx_r_p] * ww_d

            C2_11 = C3_111 * one_minus_H2O_d + C3_211 * H2O_d
            C2_12 = C3_112 * one_minus_H2O_d + C3_212 * H2O_d
            C2_21 = C3_121 * one_minus_H2O_d + C3_221 * H2O_d
            C2_22 = C3_122 * one_minus_H2O_d + C3_222 * H2O_d

            C1_1 = C2_11 * one_minus_T_d_l + C2_21 * T_d_l
            C1_2 = C2_12 * one_minus_T_d_r + C2_22 * T_d_r

            output[ww_idx] = one_minus_p_d * C1_1 + p_d * C1_2



        end

        # Set all out-of-ABSCO-range points to zero, and over-write
        # whatever result was produced above.
        @inbounds for ww_idx in eachindex(ww)
            if absco_ww_idx_left[ww_idx] == -1
                output[ww_idx] = 0
            end
        end

    else

        @turbo for ww_idx in eachindex(ww)

            idx_l_ww = absco_ww_idx_left[ww_idx]

            C3_111 = spec_cross_section[idx_l_ww, idx_l_H2O, idx_ll_T, idx_l_p]
            C3_211 = spec_cross_section[idx_l_ww, idx_r_H2O, idx_ll_T, idx_l_p]
            C3_121 = spec_cross_section[idx_l_ww, idx_l_H2O, idx_lr_T, idx_l_p]
            C3_112 = spec_cross_section[idx_l_ww, idx_l_H2O, idx_rl_T, idx_r_p]
            C3_221 = spec_cross_section[idx_l_ww, idx_r_H2O, idx_lr_T, idx_l_p]
            C3_122 = spec_cross_section[idx_l_ww, idx_l_H2O, idx_rr_T, idx_r_p]
            C3_212 = spec_cross_section[idx_l_ww, idx_r_H2O, idx_rl_T, idx_r_p]
            C3_222 = spec_cross_section[idx_l_ww, idx_r_H2O, idx_rr_T, idx_r_p]

            C2_11 = C3_111 * one_minus_H2O_d + C3_211 * H2O_d
            C2_12 = C3_112 * one_minus_H2O_d + C3_212 * H2O_d
            C2_21 = C3_121 * one_minus_H2O_d + C3_221 * H2O_d
            C2_22 = C3_122 * one_minus_H2O_d + C3_222 * H2O_d

            C1_1 = C2_11 * one_minus_T_d_l + C2_21 * T_d_l
            C1_2 = C2_12 * one_minus_T_d_r + C2_22 * T_d_r

            output[ww_idx] = one_minus_p_d * C1_1 + p_d * C1_2


        end

    end

end



"""
$(TYPEDSIGNATURES)

A wrapper function which retrieves the absorption cross_section for a vector of
wavelengths, and corresponding (scalar) values of temperature, pressure, H2O broadening
VMR. It reads the needed data from a 3D ABSCO object. Note that H2O must be supplied, even
though it is not used - this is to make the function interface compatible with the
ABSCO 4D table.

# Details

This calls the underlying, more explicit function where the spectroscopy arguments are
split up. In Julia, making those arguments explicit provides a very large performance
gain, since the compiler can infer the types.
"""
function get_cross_section_value_at!(
    output::AbstractVector,
    spec::ABSCOSpectroscopy3D,
    ww::AbstractVector,
    absco_ww_idx_left::Vector{Int},
    is_matched_ww::Bool,
    p::Unitful.Pressure,
    T::Unitful.Temperature,
    H2O::Number,
    )

    _get_absco_cross_section_value_at!(
        output,
        spec.ww,
        spec.pressures,
        spec.temperatures,
        spec.cross_section,
        ww,
        absco_ww_idx_left,
        is_matched_ww,
        p |> spec.pressures_unit |> ustrip,
        T |> spec.temperatures_unit |> ustrip,
    )

end

"""
$(TYPEDSIGNATURES)

This function is a high-performance, explicit interpolation routine that returns the
absorption cross section values of an 3D ABSCO table for a given set of coordinates in
spectral, pressure, temperature and water vapor space. The spectral coordinates must be
a vector, since in most applications we want to obtain the cross sections for all relevant
wavelengths (or wavenumbers) at once.

This function differs from a simple 3D-interpolation due to the irregular grid of the
ABSCO tables, whose temperature axis depends on the the pressure coordinate.
"""
function _get_absco_cross_section_value_at!(
    output::AbstractVector,
    spec_ww::AbstractVector,
    spec_pressures::AbstractVector,
    spec_temperatures::AbstractArray,
    spec_cross_section::AbstractArray,
    ww::AbstractVector,
    absco_ww_idx_left::Vector{Int},
    is_matched_ww::Bool,
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

    if !is_matched_ww

        @turbo for ww_idx in eachindex(ww)

            idx_l_ww = absco_ww_idx_left[ww_idx]
            idx_r_ww = idx_l_ww + 1

            ww_d = (ww[ww_idx] - spec_ww[idx_l_ww]) /
                (spec_ww[idx_r_ww] -  spec_ww[idx_l_ww])
            one_minus_ww_d = (1.0 - ww_d)


            C2_11 = spec_cross_section[idx_l_ww, idx_ll_T, idx_l_p]
            C2_21 = spec_cross_section[idx_l_ww, idx_lr_T, idx_l_p]
            C2_12 = spec_cross_section[idx_l_ww, idx_rl_T, idx_r_p]
            C2_22 = spec_cross_section[idx_l_ww, idx_rr_T, idx_r_p]

            C2_11 *= one_minus_ww_d
            C2_21 *= one_minus_ww_d
            C2_12 *= one_minus_ww_d
            C2_22 *= one_minus_ww_d

            # C3 = (h2o, temp, pressure)
            C2_11 += spec_cross_section[idx_r_ww, idx_ll_T, idx_l_p] * ww_d
            C2_21 += spec_cross_section[idx_r_ww, idx_lr_T, idx_l_p] * ww_d
            C2_12 += spec_cross_section[idx_r_ww, idx_rl_T, idx_r_p] * ww_d
            C2_22 += spec_cross_section[idx_r_ww, idx_rr_T, idx_r_p] * ww_d

            C1_1 = C2_11 * one_minus_T_d_l + C2_21 * T_d_l
            C1_2 = C2_12 * one_minus_T_d_r + C2_22 * T_d_r

            output[ww_idx] = one_minus_p_d * C1_1 + p_d * C1_2

        end

        # Set all out-of-ABSCO-range points to zero, and over-write
        # whatever result was produced above.
        @inbounds for ww_idx in eachindex(ww)
            if absco_ww_idx_left[ww_idx] == -1
                output[ww_idx] = 0
            end
        end

    else

        @turbo for ww_idx in eachindex(ww)

            idx_l_ww = absco_ww_idx_left[ww_idx]

            C2_11 = spec_cross_section[idx_l_ww, idx_ll_T, idx_l_p]
            C2_21 = spec_cross_section[idx_l_ww, idx_lr_T, idx_l_p]
            C2_12 = spec_cross_section[idx_l_ww, idx_rl_T, idx_r_p]
            C2_22 = spec_cross_section[idx_l_ww, idx_rr_T, idx_r_p]

            C1_1 = C2_11 * one_minus_T_d_l + C2_21 * T_d_l
            C1_2 = C2_12 * one_minus_T_d_r + C2_22 * T_d_r

            output[ww_idx] = one_minus_p_d * C1_1 + p_d * C1_2


        end

    end

end




"""
$(TYPEDSIGNATURES)

A wrapper function which retrieves the absorption cross_section for a vector of
wavelengths, and corresponding (scalar) values of temperature, pressure, H2O broadening
VMR. It reads the needed data from a 4D ABSCOAER object.

# Details

This calls the underlying, more explicit function where the spectroscopy arguments are
split up. In Julia, making those arguments explicit provides a very large performance
gain, since the compiler can infer the types.
"""
function get_cross_section_value_at!(
    output::AbstractVector,
    spec::ABSCOAERSpectroscopy4D,
    wl::AbstractVector,
    absco_wl_idx_left::Vector{Int},
    is_matched_wl::Bool,
    p::Unitful.Pressure,
    T::Unitful.Temperature,
    H2O::Number
    )

    _get_abscoaer_cross_section_value_at!(
        output,
        spec.ww,
        spec.pressures,
        spec.temperatures,
        spec.broadener_vmrs,
        spec.cross_section,
        wl,
        absco_wl_idx_left,
        is_matched_wl,
        p |> spec.pressures_unit |> ustrip,
        T |> spec.temperatures_unit |> ustrip,
        H2O
    )

end


"""
$(TYPEDSIGNATURES)

This function is a high-performance, explicit interpolation routine that returns the
absorption cross section values of an 4D ABSCOAER table for a given set of coordinates in
spectral, pressure, temperature and water vapor space. The spectral coordinates must be
a vector, since in most applications we want to obtain the cross sections for all relevant
wavelengths (or wavenumbers) at once.

"""
function _get_abscoaer_cross_section_value_at!(
    output::AbstractVector,
    spec_ww::AbstractVector,
    spec_pressures::AbstractVector,
    spec_temperatures::AbstractArray,
    spec_H2O::AbstractVector,
    spec_cross_section::AbstractArray,
    ww::AbstractVector,
    absco_ww_idx_left::Vector{Int},
    is_matched_ww::Bool,
    p,
    T,
    H2O
)

    idx_p = 1 # no interpolation in pressure space
    idx_l_H2O = 1
    idx_l_T = 1

    # There is no interpolation in pressure space, we simply pick the layer in which
    # `p` falls.
    if p < spec_pressures[1]
        idx_p = 1
    elseif p > spec_pressures[end]
        idx_p = size(spec_pressures, 1) - 1
    else
        for j in 1:length(spec_pressures) - 1
            @inbounds if (p >= spec_pressures[j]) & (p <= spec_pressures[j+1])
                idx_p = j
                break
            end
        end
    end

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

    # The T grid has NaNs! We must figure out where they are
    finite_T = isfinite.(spec_temperatures[:, idx_p])
    jmin = findfirst(finite_T)
    jmax = findlast(finite_T)

    if T < spec_temperatures[jmin, idx_p]
        idx_l_T = jmin
    elseif T > spec_temperatures[jmax, idx_p]
        idx_l_T = jmax - 1
    else
        for j in jmin:jmax - 1
            @inbounds if (T >= spec_temperatures[j, idx_p]) &
                (T <= spec_temperatures[j+1, idx_p])
                idx_l_T = j
                break
            end
        end
    end
    idx_r_T = idx_l_T + 1


    T_d = (T - spec_temperatures[idx_l_T, idx_p]) /
        (spec_temperatures[idx_r_T, idx_p] - spec_temperatures[idx_l_T, idx_p])
    H2O_d = (H2O - spec_H2O[idx_l_H2O]) /
        (spec_H2O[idx_r_H2O] - spec_H2O[idx_l_H2O])

    one_minus_T_d = one(typeof(T_d)) - T_d
    one_minus_H2O_d = one(typeof(H2O_d)) - H2O_d


    if !is_matched_ww


        @turbo for ww_idx in eachindex(ww)

            idx_l_ww = absco_ww_idx_left[ww_idx]
            idx_r_ww = idx_l_ww + 1

            ww_d = (ww[ww_idx] - spec_ww[idx_l_ww]) /
                (spec_ww[idx_r_ww] -  spec_ww[idx_l_ww])
            one_minus_ww_d = (1.0 - ww_d)


            # produce the ww-interpolated values first
            C2_11 = spec_cross_section[idx_l_H2O, idx_p, idx_l_T, idx_l_ww]
            C2_21 = spec_cross_section[idx_r_H2O, idx_p, idx_l_T, idx_l_ww]
            C2_12 = spec_cross_section[idx_l_H2O, idx_p, idx_r_T, idx_l_ww]
            C2_22 = spec_cross_section[idx_r_H2O, idx_p, idx_r_T, idx_l_ww]

            C2_11 *= one_minus_ww_d
            C2_21 *= one_minus_ww_d
            C2_12 *= one_minus_ww_d
            C2_22 *= one_minus_ww_d

            C2_11 = spec_cross_section[idx_l_H2O, idx_p, idx_l_T, idx_r_ww] * ww_d
            C2_21 = spec_cross_section[idx_r_H2O, idx_p, idx_l_T, idx_r_ww] * ww_d
            C2_12 = spec_cross_section[idx_l_H2O, idx_p, idx_r_T, idx_r_ww] * ww_d
            C2_22 = spec_cross_section[idx_r_H2O, idx_p, idx_r_T, idx_r_ww] * ww_d

            # then move on to get the interpolation in H2O and T

            # first collapse H2O
            C1_1 = C2_11 * one_minus_H2O_d + C2_21 * H2O_d
            C1_2 = C2_12 * one_minus_H2O_d + C2_22 * H2O_d

            # then T
            output[ww_idx] = one_minus_T_d * C1_1 + T_d * C1_2

        end

        # Set all out-of-ABSCO-range points to zero, and over-write
        # whatever result was produced above.
        @inbounds for ww_idx in eachindex(ww)
            if absco_ww_idx_left[ww_idx] == -1
                output[ww_idx] = 0
            end
        end

    else

        @turbo for ww_idx in eachindex(ww)

            idx_l_ww = absco_ww_idx_left[ww_idx]

            # indices: H2O, T
            C2_11 = spec_cross_section[idx_l_H2O, idx_p, idx_l_T, idx_l_ww]
            C2_21 = spec_cross_section[idx_r_H2O, idx_p, idx_l_T, idx_l_ww]
            C2_12 = spec_cross_section[idx_l_H2O, idx_p, idx_r_T, idx_l_ww]
            C2_22 = spec_cross_section[idx_r_H2O, idx_p, idx_r_T, idx_l_ww]

            # first collapse H2O
            C1_1 = C2_11 * one_minus_H2O_d + C2_21 * H2O_d
            C1_2 = C2_12 * one_minus_H2O_d + C2_22 * H2O_d

            # then T
            output[ww_idx] = one_minus_T_d * C1_1 + T_d * C1_2

        end

    end

end


"""
$(TYPEDSIGNATURES)

A wrapper function which retrieves the absorption cross_section for a vector of
wavelengths, and corresponding (scalar) values of temperature, pressure, H2O broadening
VMR. It reads the needed data from a 3D ABSCOAER object.

# Details

This calls the underlying, more explicit function where the spectroscopy arguments are
split up. In Julia, making those arguments explicit provides a very large performance
gain, since the compiler can infer the types.
"""
function get_cross_section_value_at!(
    output::AbstractVector,
    spec::ABSCOAERSpectroscopy3D,
    wl::AbstractVector,
    absco_wl_idx_left::Vector{Int},
    is_matched_wl::Bool,
    p::Unitful.Pressure,
    T::Unitful.Temperature,
    H2O::Number
    )

    _get_abscoaer_cross_section_value_at!(
        output,
        spec.ww,
        spec.pressures,
        spec.temperatures,
        spec.cross_section,
        wl,
        absco_wl_idx_left,
        is_matched_wl,
        p |> spec.pressures_unit |> ustrip,
        T |> spec.temperatures_unit |> ustrip,
    )

end


"""
$(TYPEDSIGNATURES)

This function is a high-performance, explicit interpolation routine that returns the
absorption cross section values of an 4D ABSCOAER table for a given set of coordinates in
spectral, pressure, temperature and water vapor space. The spectral coordinates must be
a vector, since in most applications we want to obtain the cross sections for all relevant
wavelengths (or wavenumbers) at once.

"""
function _get_abscoaer_cross_section_value_at!(
    output::AbstractVector,
    spec_ww::AbstractVector,
    spec_pressures::AbstractVector,
    spec_temperatures::AbstractArray,
    spec_cross_section::AbstractArray,
    ww::AbstractVector,
    absco_ww_idx_left::Vector{Int},
    is_matched_ww::Bool,
    p,
    T
)

    idx_p = 1 # no interpolation in pressure space
    idx_l_T = 1

    # There is no interpolation in pressure space, we simply pick the layer in which
    # `p` falls.
    if p < spec_pressures[1]
        idx_p = 1
    elseif p > spec_pressures[end]
        idx_p = size(spec_pressures, 1) - 1
    else
        for j in 1:length(spec_pressures) - 1
            @inbounds if (p >= spec_pressures[j]) & (p <= spec_pressures[j+1])
                idx_p = j
                break
            end
        end
    end

    # The T grid has NaNs! We must figure out where they are
    finite_T = isfinite.(spec_temperatures[:, idx_p])
    jmin = findfirst(finite_T)
    jmax = findlast(finite_T)

    if T < spec_temperatures[jmin, idx_p]
        idx_l_T = jmin
    elseif T > spec_temperatures[jmax, idx_p]
        idx_l_T = jmax - 1
    else
        for j in jmin:jmax - 1
            @inbounds if (T >= spec_temperatures[j, idx_p]) &
                (T <= spec_temperatures[j+1, idx_p])
                idx_l_T = j
                break
            end
        end
    end
    idx_r_T = idx_l_T + 1


    T_d = (T - spec_temperatures[idx_l_T, idx_p]) /
        (spec_temperatures[idx_r_T, idx_p] - spec_temperatures[idx_l_T, idx_p])

    one_minus_T_d = one(typeof(T_d)) - T_d


    if !is_matched_ww

        @turbo for ww_idx in eachindex(ww)

            idx_l_ww = absco_ww_idx_left[ww_idx]
            idx_r_ww = idx_l_ww + 1

            ww_d = (ww[ww_idx] - spec_ww[idx_l_ww]) /
                (spec_ww[idx_r_ww] -  spec_ww[idx_l_ww])
            one_minus_ww_d = (1.0 - ww_d)

            # produce the ww-interpolated values first
            C1_1 = spec_cross_section[idx_p, idx_l_T, idx_l_ww]
            C1_2 = spec_cross_section[idx_p, idx_r_T, idx_l_ww]

            C1_1 *= one_minus_ww_d
            C1_2 *= one_minus_ww_d

            C1_1 = spec_cross_section[idx_p, idx_l_T, idx_r_ww] * ww_d
            C1_2 = spec_cross_section[idx_p, idx_r_T, idx_r_ww] * ww_d

            # then move on to get the interpolation in T
            output[ww_idx] = one_minus_T_d * C1_1 + T_d * C1_2

        end

        # Set all out-of-ABSCO-range points to zero, and over-write
        # whatever result was produced above.
        @inbounds for ww_idx in eachindex(ww)
            if absco_ww_idx_left[ww_idx] == -1
                output[ww_idx] = 0
            end
        end

    else

        for ww_idx in eachindex(ww)

            idx_l_ww = absco_ww_idx_left[ww_idx]

            # indices: (p), T, (λ, ν)
            C1_1 = spec_cross_section[idx_p, idx_l_T, idx_l_ww]
            C1_2 = spec_cross_section[idx_p, idx_r_T, idx_l_ww]

            # then T
            output[ww_idx] = one_minus_T_d * C1_1 + T_d * C1_2
            if isnan(output[ww_idx])
                output[ww_idx] = 0
            end


        end

    end



end