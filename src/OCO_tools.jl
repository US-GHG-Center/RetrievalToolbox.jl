"""
Extracts the frame number from an integer-valued
OCO sounding ID using the simple divide-by-ten method
and then grabbing the remainder; the last digit of the
sounding ID represents the footprint number.


$(TYPEDSIGNATURES)

"""
function OCO_get_footprint_from_sounding_id(snid::Integer)

    return snid % 10

end

"""
Converts the band index (1,2,3) to its respective
'human-readable' name, (e.g. 'o2'), that is used in
various variable names.

$(TYPEDSIGNATURES)

"""
function OCO_band_index_to_string(idx::Integer)

    bdict = Dict(
        1 => "o2",
        2 => "weak_co2",
        3 => "strong_co2"
    )

    return bdict[idx]

end


"""
Returns a slicer object to read data from OCO L2 files in
the usual (Footprint, Frame) shape. If `sounding_id` is an
integer-valued sounding ID, the returned slicer object will
be the position of the sounding within the (Footprint, Frame)
array shape, as returned by the `OCO_get_idx_from_sounding_id`
function.
If `sounding_id` is `nothing`, the function returns a `(:,:)`
slicer, such that **all** scenes will can be extracted. Requires
a valid `HDF5.File` object to look into and expects the
`SoundingGeometry/sounding_id` field within that file.

$(TYPEDSIGNATURES)

# Example usage


"""
function OCO_one_or_all(h5::HDF5.File, sounding_id::Union{Integer, Nothing})

    if sounding_id isa Integer
        @debug "Pulling atmosphere for a single ID: $(sounding_id)"
        # Get the specific array index for the requested sounding ID
        idx = OCO_get_idx_from_sounding_id(h5, sounding_id).I
    elseif isnothing(sounding_id)
        @debug "Pull all atmoshperes."
        # This is a 2D array slicer to grab everything
        idx = (:,:)
    else
        @error "`sounding_id` must be an Integer or a Nothing."
    end

    return idx

end


function OCO_get_idx_from_sounding_id(h5::HDF5.File, sounding_id::Integer)

    sounding_idx = findall(h5["/SoundingGeometry/sounding_id"][:,:] .== sounding_id)
    @assert length(sounding_idx) == 1 "Sounding ID $(sounding_id) not found!"

    return sounding_idx[1]

end

function OCO_get_footprint_and_frame_from_snid(
    h5::HDF5.File;
    sounding_id::Union{Integer, Nothing}=nothing
    )

    idx = OCO_one_or_all(h5, sounding_id)

    if sounding_id isa Integer
        return Dict(sounding_id => (idx[1], idx[2]))
    else
        sounding_ids = h5["/SoundingGeometry/sounding_id"][:,:]
        fp_frame_dict = Dict{Integer, Tuple{Integer, Integer}}()

        for idx2 in CartesianIndices(sounding_ids)
            fp_frame_dict[sounding_ids[idx2]] = idx2.I
        end
    end

    return fp_frame_dict
end

"""


"""
function OCO_pull_atmosphere_from_met(
    met_h5::HDF5.File;
    sounding_id::Union{Integer, Nothing}=nothing
    )

    # Grab indexer object
    idx = OCO_one_or_all(met_h5, sounding_id)

    if "Meteorology" in keys(met_h5)

        _h5 = met_h5["/Meteorology/"]
        _h5_aer = met_h5["/Aerosol/"]
        # Take data out of MET file

        psurf = _h5["surface_pressure_met"][idx...]u"Pa"
        ptropo = _h5["blended_tropopause_pressure_met"][idx...]u"Pa"
        qflag = _h5["meteorology_flag"][idx...]

        pressure_levels_met = _h5["vector_pressure_levels_met"][:,idx...]u"Pa"
        temperature_levels_met = _h5["temperature_profile_met"][:,idx...]u"K"
        specific_humidity_levels_met = _h5["specific_humidity_profile_met"][:,idx...]

        windspeed_u_met = _h5["windspeed_u_met"][idx...]u"m/s"
        windspeed_v_met = _h5["windspeed_v_met"][idx...]u"m/s"

        # ?, height, width, aod
        aerosol_gauss_params = _h5_aer["composite_aod_gaussian_met"][:,:, idx...]
        aerosol_sort = _h5_aer["composite_aod_sort_index_met"][:, idx...]

    elseif "ECMWF" in keys(met_h5)

        _h5 = met_h5["/ECMWF/"]

        psurf = _h5["surface_pressure_ecmwf"][idx...]u"Pa"
        #ptropo = _h5["blended_tropopause_pressure_ecwmf"][idx]u"Pa"
        ptropo = similar(psurf)
        ptropo[:,:] .= 0.0u"Pa"

        pressure_levels_met = _h5["vector_pressure_levels_ecmwf"][:,idx...]u"Pa"
        temperature_levels_met = _h5["temperature_profile_ecmwf"][:,idx...]u"K"
        specific_humidity_levels_met = _h5["specific_humidity_profile_ecmwf"][:,idx...]

        windspeed_u_met = _h5["windspeed_u_ecmwf"][idx...]u"m/s"
        windspeed_v_met = _h5["windspeed_v_ecmwf"][idx...]u"m/s"

    end


    # Now push the retrieved data into a dictionary
    # and return.

    if sounding_id isa Integer

        # If only one atmosphere is requested,
        # we return one object
        return Dict(sounding_id =>
                Dict(
                "qflag_met" => qflag,
                "pressure_levels_met" => pressure_levels_met,
                "temperature_profile_met" => temperature_levels_met,
                "specific_humidity_profile_met" => specific_humidity_levels_met,
                "surface_pressure_met" => psurf,
                "tropopause_pressure_met" => ptropo,
                "windspeed_u_met" => windspeed_u_met,
                "windspeed_v_met" => windspeed_v_met,
                "aerosol_gauss_params_met" => aerosol_gauss_params,
                "aerosol_sort_met" => aerosol_sort,
            )
        )

    elseif isnothing(sounding_id)

        sounding_ids = met_h5["/SoundingGeometry/sounding_id"][:,:]
        met_dict = Dict{Integer, Dict{String, Any}}()

        for idx2 in CartesianIndices(sounding_ids)

            snid = sounding_ids[idx2]
            met_dict[snid] = Dict{String, Any}()

            met_dict[snid]["qflag_met"] = qflag[idx2]

            met_dict[snid]["surface_pressure_met"] = psurf[idx2]
            met_dict[snid]["tropopause_pressure_met"] = ptropo[idx2]

            met_dict[snid]["pressure_levels_met"] =
                @view pressure_levels_met[:, idx2]
            met_dict[snid]["temperature_profile_met"] =
                @view temperature_levels_met[:, idx2]
            met_dict[snid]["specific_humidity_profile_met"] =
                @view specific_humidity_levels_met[:, idx2]

            met_dict[snid]["windspeed_u_met"] = windspeed_u_met[idx2]
            met_dict[snid]["windspeed_v_met"] = windspeed_v_met[idx2]
            met_dict[snid]["aerosol_gauss_params_met"] = aerosol_gauss_params[:,:,idx2]
            met_dict[snid]["aerosol_sort_met"] = aerosol_sort[:,idx2]

        end

        return met_dict

    end

end

function OCO_pull_sounding_time(
    h5::HDF5.File;
    sounding_id::Union{Integer, Nothing}=nothing
    )

    idx = OCO_one_or_all(h5, sounding_id)

    this_time = @. DateTime(1993, 1, 1) + Dates.Second(
        round(h5["SoundingGeometry/sounding_time_tai93"][idx...])
    )

    if sounding_id isa Integer

        return Dict(sounding_id => this_time)

    else

        sounding_ids = h5["/SoundingGeometry/sounding_id"][:,:]
        time_dict = Dict{Integer, DateTime}()

        for idx2 in CartesianIndices(sounding_ids)

            time_dict[sounding_ids[idx2]] = this_time[idx2]

        end

        return time_dict


    end


end


function OCO_pull_sounding_location(
    h5::HDF5.File;
    sounding_id::Union{Integer, Nothing}=nothing
    )

    idx = OCO_one_or_all(h5, sounding_id)

    #this_time = @. DateTime(1993, 1, 1) + Dates.Second(
    #    round(h5["SoundingGeometry/sounding_time_tai93"][idx...]))

    lon = h5["SoundingGeometry/sounding_longitude"][idx...]u"°" |> ustrip
    lat = h5["SoundingGeometry/sounding_latitude"][idx...]u"°" |> ustrip
    alt = h5["SoundingGeometry/sounding_altitude"][idx...]

    if sounding_id isa Integer

        return Dict(sounding_id =>
            EarthLocation(
                lon, lat, alt, u"m"
            )
        )

    elseif isnothing(sounding_id)

        sounding_ids = h5["/SoundingGeometry/sounding_id"][:,:]
        loc_dict = Dict{Integer, EarthLocation}()

        for idx2 in CartesianIndices(sounding_ids)

            loc_dict[sounding_ids[idx2]] = EarthLocation(
                lon[idx2], lat[idx2], alt[idx2], u"m"
            )

        end

        return loc_dict

    end

end

function OCO_pull_CO2_prior(
    h5::HDF5.File;
    sounding_id::Union{Integer, Nothing}=nothing
    )

    idx = OCO_one_or_all(h5, sounding_id)

    # Make these profiles in "ppm"
    if haskey(h5, "CO2Prior")
        co2pr = h5["/CO2Prior/co2_prior_profile_cpr"][:,idx...] .* 1e6 * u"ppm"
    else
        co2pr = nothing
    end
    if haskey(h5, "CH4Prior")
        ch4pr = h5["/CH4Prior/ch4_prior_profile_cpr"][:,idx...] .* 1e9 * u"ppb"
    else
        ch4pr = nothing
    end

    if sounding_id isa Integer

        return Dict(sounding_id =>
            Dict(
                "co2_prior" => co2pr,
                "ch4_prior" => ch4pr
            )
        )

    elseif isnothing(sounding_id)


        sounding_ids = h5["/SoundingGeometry/sounding_id"][:,:]
        cpr_dict = Dict{Integer, Dict{String, Any}}()

        for idx2 in CartesianIndices(sounding_ids)

            snid = sounding_ids[idx2]

            cpr_dict[snid] = Dict()
            if !isnothing(co2pr)
                cpr_dict[snid]["co2_prior"] = co2pr[:,idx2]
            end
            if !isnothing(ch4pr)
                cpr_dict[snid]["ch4_prior"] = ch4pr[:,idx2]
            end
        end

        return cpr_dict

    end

end

function OCO_pull_observer(
    h5::HDF5.File;
    sounding_id::Union{Integer, Nothing}=nothing
    )

    idx = OCO_one_or_all(h5, sounding_id)

    # Viewing angles in OCO files are already in degrees.
    this_vza = h5["SoundingGeometry/sounding_zenith"][idx...]
    this_vaa = h5["SoundingGeometry/sounding_azimuth"][idx...]

    this_velocity = h5["FrameGeometry/spacecraft_velocity"][:,idx[2]]
    this_position = h5["FrameGeometry/spacecraft_position"][:,idx[2]]

    if sounding_id isa Integer

        return Dict(sounding_id =>
            SatelliteObserver(
                this_vza,
                this_vaa,
                this_velocity,
                this_position
            )
        )

    else

        sounding_ids = h5["/SoundingGeometry/sounding_id"][:,:]

        obs_dict = Dict{Integer, SatelliteObserver}()

        for idx2 in CartesianIndices(sounding_ids)

            obs_dict[sounding_ids[idx2]] = SatelliteObserver(
                this_vza[idx2],
                this_vaa[idx2],
                this_velocity[:, idx2.I[2]],
                this_position[:, idx2.I[2]]
            )

        end

        return obs_dict

    end

end

function OCO_calculate_noise(
    rad::AbstractVector,
    snr_coef::AbstractArray,
    mms::Number
        )

    noise = similar(rad)
    @turbo for i in eachindex(noise)
        noise[i] = mms / 100.0 * sqrt(
            abs(100.0 * rad[i]) / mms * snr_coef[1,i]^2
            + snr_coef[2,i]^2
        )
    end

    return noise

end