"""
Pretty printing for OCO solar model types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sm::OCOHDFSolarModel)

    println(io, "OCO HDF solar model from file: $(sm.file_name)")

end

"""
Brief pretty printing for OCO solar model types

$(SIGNATURES)
"""
function show(io::IO, sm::OCOHDFSolarModel)

    print(io, "OCOHDFSolarModel: $(sm.file_name)")

end



"""
Calculates the down-sampled solar spectrum at the
high-resolution grid specified within the spectral
window `swin` and saves it in `rt.hires_solar`.

$(TYPEDSIGNATURES)

# Details

In this function, the solar Doppler shift is considered via
sampling the original solar spectrum at acccording wavelengths.
The Doppler factor is defined as the relative velocity between
the observation point on Earth (for Earth-Backscatter spectra),
and a negative sign indicates the observation point moving closer
to the sun, as a fraction of the speed of light in vacuum.
"""
function calculate_solar_irradiance!(
    rt::AbstractRTMethod,
    swin::AbstractSpectralWindow,
    solar_model::AbstractSolarModel;
    doppler_factor=nothing
)

    # If no Doppler factor is provided, we could
    # calculate it right here using rt.scene.location
    # and rt.scene.time.
    if isnothing(doppler_factor)

        doppler_factor = calculate_solar_doppler_shift(
            rt.scene.location,
            rt.scene.time
        )

    end

    # Doppler effect depends on the spectral unit.

    if solar_model.ww_unit isa Unitful.LengthUnits
        # Very cheeky way of moving from λ -> λ * (1 + Doppler)
        @turbo for i in eachindex(solar_model.ww)
            solar_model.ww[i] *= (1 + doppler_factor)
        end
    elseif solar_model.ww_unit isa Unitful.WavenumberUnits
        # Very cheeky way of moving from ν -> ν / (1 + Doppler)
        @turbo for i in eachindex(solar_model.ww)
            solar_model.ww[i] /= (1 + doppler_factor)
        end
    end

    #=
    We must sample the solar spectrum at Doppler-influenced
    wavelengths and produce following result:
        y = transmittance(λ * (1 + Doppler)) * continuum(λ * (1 + Doppler))
    and in a final step, multiply by the (per-wavelength) scaler
        y = y * solar_scaler
    =#

    if solar_model isa OCOHDFSolarModel
        # Sample the solar spectrum at our Doppler-influenced
        # retrieval wavelength grid
        # (and calculate continuum * transmittance in one step)
        pwl_value_1d_axb!(
            solar_model.ww,
            solar_model.transmittance,
            solar_model.continuum,
            swin.ww_grid,
            rt.hires_solar.I,
        )

    elseif solar_model isa TSISSolarModel
        # Sample the TSIS irradiance (includes both transmittance
        # and continuum) and store in hires_solar.I

        pwl_value_1d!(
            solar_model.ww,
            solar_model.irradiance,
            swin.ww_grid,
            rt.hires_solar.I
        )

    elseif solar_model isa ListSolarModel

        pwl_value_1d!(
            solar_model.ww,
            solar_model.irradiance,
            swin.ww_grid,
            rt.hires_solar.I
        )

    else

        throw(error("Solar model type $(typeof(solar_model)) not implemented!"))

    end

    # Apply the solar scaler
    @turbo for i in eachindex(rt.hires_solar.I)
        rt.hires_solar.S[i,1] *= rt.solar_scaler[i]
    end

    # Restore original solar model grid
    if solar_model.ww_unit isa Unitful.LengthUnits
        # Very cheeky way of moving back from λ * (1 + Doppler) -> λ
        @turbo for i in eachindex(solar_model.ww)
            solar_model.ww[i] /= (1 + doppler_factor)
        end
    elseif solar_model.ww_unit isa Unitful.WavenumberUnits
        # Very cheeky way of moving back from ν * (1 + Doppler) -> ν
        @turbo for i in eachindex(solar_model.ww)
            solar_model.ww[i] *= (1 + doppler_factor)
        end
    end

end

"""
    Converts solar model data from "ph/s/m^2/µm" to "W/m^2/µm"
"""
function convert_solar_model_to_W!(s::OCOHDFSolarModel)

    if s.irradiance_unit == u"ph/s/m^2/µm"
        @debug "[SOLAR] Solar model in units of ph/s/m2/µm - converting to W/m2/µm!"
        @views s.continuum[:] .*= ustrip.(Ref(u"W"),
            1.0u"s^-1" .* SPEED_OF_LIGHT ./ (s.ww[:] .* u"µm") .* PLANCK
        )

        s.irradiance_unit = u"W/m^2/µm"

    else

        @warn "Units already in W/m^2/µm - skipping!"

    end

end

"""
    Converts solar model data from "W/m^2/nm" to "ph/s/m^2/µm"
"""
function convert_solar_model_to_photons!(s::TSISSolarModel)

    if s.irradiance_unit == u"W/m^2/nm"
        @debug "[SOLAR] Solar model in units of W/m^2/nm - converting to ph/s/m^2/µm!"
        @views s.irradiance[:] .*= ustrip.(u"m^-2 * µm^-1", # We want this in 1/m2 1/µm
        s.irradiance_unit * 1.0u"s" ./ SPEED_OF_LIGHT .* (s.ww[:] .* u"nm") ./ PLANCK
        )

        s.irradiance_unit = u"ph/s/m^2/µm"

    else

        @warn "Units already in ph/s/m^2/µm - skipping!"

    end

    return true

end


"""
Reads a JPL/OCO-type solar model HDF5 file and returns
a `OCOHDFSolarModel` object.

$(TYPEDSIGNATURES)
"""
function OCOHDFSolarModel(
    filename::String,
    band_number::Integer;
    spectral_unit=:Wavelength
    )

    @assert isfile(filename) "File $(filename) is not a regular file!"

    @debug "[SOLAR] Opening up Solar HDF file $(filename)"
    h5 = h5open(filename, "r")

    h5g_c = h5["Solar/Continuum/Continuum_$(band_number)"]
    h5g_a = h5["Solar/Absorption/Absorption_$(band_number)"]

    # Grab Absorption and Continuum
    # (values are hard-coded in the files, but not the units, so
    #  these could theoretically change in the future)

    if spectral_unit == :Wavelength

        # This is in microns
        ww_unit = u"µm"
        rad_unit = u"ph/s/m^2/µm"

        cont_ww = (1e4 ./ h5g_c["wavenumber"][:][end:-1:1])
        trans_ww = (1e4 ./ h5g_a["wavenumber"][:][end:-1:1])
        # Continuum-level values are provided in units of
        # ph/s/m2/µm, hence no further conversion between spectral radiance
        # per wavelength and spectral radiance per wavenumber is necessary.
        cont_val = h5g_c["spectrum"][:][end:-1:1]
        trans_val = h5g_a["spectrum"][:][end:-1:1]

    elseif spectral_unit == :Wavenumber

        ww_unit = u"cm^-1"
        rad_unit = u"W/m^2/cm^-1"

        # Need microns for calculations/conversions
        cont_microns = (1e4 ./ h5g_c["wavenumber"][:])

        # This is in cm^-1
        cont_ww = h5g_c["wavenumber"][:]
        trans_ww = h5g_a["wavenumber"][:]

        # Here, we need to convert from ph/s/m2/µm into
        # W/m2/cm^-1

        cont_val = h5g_c["spectrum"][:]
        # 1) convert ph/s into W
        @views cont_val[:] .*= ustrip.(Ref(u"W"),
            1.0u"s^-1" .* SPEED_OF_LIGHT ./ (cont_microns .* u"µm") .* PLANCK
        )

        #2) convert  W/m2/µm into W/m2/cm^-1
        @views cont_val[:] ./= (1e4 ./ cont_microns) .^ 2

        trans_val = h5g_a["spectrum"][:]

    end


    # Build a continuum polynomial through fitting the data
    cont_poly = Polynomials.fit(
        cont_ww,
        cont_val,
        4, #length(cont_val) - 1
    )
    # Evaluate the polynomial at all spectral window wavelengths
    # (this is reasonable since the continuum is so smooth)
    cont_resampled = cont_poly.(trans_ww)

    # Close up HDF file, all done
    close(h5)

    # Return solar model object
    return OCOHDFSolarModel(
        filename,
        band_number,
        trans_ww,
        trans_val,
        cont_resampled,
        ww_unit, # Wavelength unit
        rad_unit# Radiance unit
    )

end

"""
    Reads the full TSIS file into memory
"""
function TSISSolarModel(
    fname::String;
    spectral_unit=:Wavelength
    )

    irradiance = ncread(fname, "SSI")
    # The TSIS file has units such as
    # W m-2 nm-1, so we need to add the caret
    # such that Unitful understands them

    # Read the irradiance unit from the NetCDF file
    irradiance_u  = replace(ncgetatt(fname, "SSI", "units"), "-" => "^-", " " => "*")
    irradiance_unit = uparse(irradiance_u)

    # Read the wavelength unit from the NetCDF file
    wavelength = ncread(fname, "Vacuum Wavelength")
    wavelength_unit = uparse(ncgetatt(fname, "Vacuum Wavelength", "units"))

    # Turn wavelength into microns
    ww = ustrip.(Ref(u"µm"), wavelength * wavelength_unit)
    ww_unit = u"µm"

    # Convert wavelength to µm
    irradiance = ustrip.(Ref(u"W/m^2/µm"), irradiance * irradiance_unit)


    if spectral_unit == :Wavelength
        # Nothing to do here

    elseif spectral_unit == :Wavenumber

        # Turn irradiance into W/m^2/cm^-1 and reverse
        # array order to make them in increasing wavenumbers
        irradiance ./= (1e4 ./ ww) .^ 2
        irradiance = irradiance[end:-1:1]
        # Set the new irradiance units
        irradiance_unit = u"W/m^2/cm^-1"

        # Turn µm into cm^-1 and reverse array
        ww = 1e4 ./ ww[end:-1:1]
        ww_unit = u"cm^-1"

    else
        @error "Unknown parameter for spectral_unit: $(spectral_unit)"
    end

    # Create solar model object
    return TSISSolarModel(
        fname,
        ww,
        ww_unit,
        irradiance,
        irradiance_unit
    )

end


"""
Calculates the solar Doppler shift factor for a given
`EarthLocation` and `DateTime`.

$(TYPEDSIGNATURES)

"""
function calculate_solar_doppler_shift(
    loc::EarthLocation,
    time::DateTime,
    )

    @warn "This is not yet implemented."
    return 0.0

end
