"""
$(TYPEDSIGNATURES)

Pretty printing for EarthAtmosphereBuffer
"""
function show(io::IO, ::MIME"text/plain", buf::EarthAtmosphereBuffer)
    # Just print this for now
    println(io, "EarthAtmosphereBuffer")
    println(io, "=====================")
    println(io, "Windows: ")
    for swin in buf.spectral_window
        println("   $(swin)")
    end
end

"""
$(TYPEDSIGNATURES)

Pretty printing for Beer-Lambert RT
"""
function show(io::IO, buf::EarthAtmosphereBuffer)
    # Just print this for now
    println(io, "EarthAtmosphereBuffer")

end


"""
$(TYPEDSIGNATURES)

Helper function to populate an `EarthAtmosphereBuffer`, which also includes an
`EarthAtmosphere` and the corresponding `OpticalProperties` with correctly sized arrays.

# Details

TODO: a lot of documentation needs to go here.

"""
function EarthAtmosphereBuffer(
    sv::AbstractStateVector,
    spectral_windows,
    surface_types::Vector{Tuple},
    atmospheric_elements,
    solar_models::Dict{<:AbstractSpectralWindow,<:AbstractSolarModel},
    RT_models::Vector{Symbol},
    RadType::Type{<:Radiance},
    rt_buf::AbstractRTBuffer,
    inst_buf,
    N_level::Integer,
    N_met_level::Integer,
    T::Type{<:AbstractFloat}
    )

    N_layer = N_level - 1
    N_met_layer = N_met_level - 1

    # We can supply both a single spectral window
    # or a list of them. Maybe this should be a dict?
    if spectral_windows isa Vector
        swins = spectral_windows
    else
        swins = [spectral_windows]
    end

    if atmospheric_elements isa Vector
        atm_elements = atmospheric_elements
    else
        atm_elements = [atmospheric_elements]
    end

    # Create an empty EarthAtmosphere
    # This is shared between ALL spectral windows
    atmosphere = create_empty_EarthAtmosphere(
        N_level,
        N_met_level,
        T,
        pressure_unit=u"Pa",
        met_pressure_unit=u"Pa",
        temperature_unit=u"K",
        specific_humidity_unit=u"kg/kg",
        altitude_unit=u"m",
        gravity_unit=u"m/s^2"
    )

    # Atmosphere elements into `atmosphere` object
    for atm in atm_elements
        push!(atmosphere.atm_elements, atm)
    end

    surfaces = Dict{SpectralWindow, AbstractSurface}()

    # Loop through all surfaces provided and create appropriate objects
    for (i_surf, swin) in enumerate(swins)

        # Lambertian surface: second argument is the polynomial order we wish to use
        if surface_types[i_surf][1] in [:Lambertian, :Lambert]

            # We use a BRDF kernel if the RT model of choice is XRTM
            if RT_models[i_surf] == :XRTM
                kern = LambertianPolynomialKernel(swin,
                    repeat([0.0], surface_types[i_surf][2])
                )

                surfaces[swin] = BRDFSurface([kern])
            # We use the non-BRDF Lambertian when doing Beer-Lambert RT
            elseif RT_models[i_surf] == :BeerLambert
                surfaces[swin] = LambertianPolynomialSurface(
                    swin, zeros(surface_types[i_surf][2])
                )
            end
        end

        # RPV BRDF surface: second argument is asymmetry, third is anisotropy
        if surface_types[i_surf][1] in [:RPV, :Rahman]
            # Create a single-kernel BRDF surface:
            kern = RPVPolynomialKernel(swin,
                repeat([0.0], surface_types[i_surf][2]),
                surface_types[i_surf][3], # hotspot
                surface_types[i_surf][4], # asymmetry
                surface_types[i_surf][5], # anisotropy
            )

            surfaces[swin] = BRDFSurface([kern])
        end

        if surface_types[i_surf][1] in [:NoSurface, :nothing, :none]
            surfaces[swin] = NoSurface()
        end
    end

    # Create a satellite observer
    # Shared between all spectral windows
    observer = SatelliteObserver(
        zero(T), zero(T),
        zeros(T, 3), zeros(T, 3)
    )

    # Create an Earth location
    # Shared between all spectral windows
    location = EarthLocation(
        zero(T), zero(T), zero(T), u"m"
    )

    # Create an Earth scene
    # Shared between all spectral windows
    earth_scene = EarthScene(
        atmosphere,
        surfaces,
        observer,
        location,
        zero(T),
        zero(T),
        DateTime(0)
    )

    # Allocate arrays for optical properties, which
    # are derived from the atmospheric elemements

    # These are created for each spectral window!

    optical_properties = Dict{
        AbstractSpectralWindow,
        EarthAtmosphereOpticalProperties
    }()


    for swin in swins

        # Create an empty EarthAtmosphereOpticalProperties object and push into dict
        # (Note that we must supply the explicit number type T here, otherwise the
        #  default will be Float64)
        optical_properties[swin] = create_empty_EarthAtmosphereOpticalProperties(
            swin, earth_scene, sv, RadType, T
        )

    end

    # Allocate RT objects for each spectral window
    rt = Dict{AbstractSpectralWindow, AbstractRTMethod}()

    for (i_swin,swin) in enumerate(swins)

        # Create a high-resolution solar irradiance object
        hires_solar = RadType(T, swin.N_hires)

        # Create a high-resolution instrument radiance object
        hires_radiance = RadType(T, swin.N_hires)

        # .. and the associated Jacobians for each SVE.
        if sv isa RetrievalStateVector
            hires_jacobians = Dict(
                sve => RadType(T, swin.N_hires)
                for sve in sv.state_vector_elements
            )
        elseif sv isa ForwardModelStateVector
            hires_jacobians = nothing
        end

        # An array to hold the per-wavelength solar continuum scaling
        # factor values.
        solar_scaler = ones(T, swin.N_hires)

        if RT_models[i_swin] == :BeerLambert

            this_rt = BeerLambertRTMethod(
                earth_scene,
                optical_properties[swin],
                solar_models[swin],
                sv,
                hires_solar,
                hires_radiance,
                hires_jacobians,
                solar_models[swin].irradiance_unit / u"sr",
                solar_scaler
            )

        elseif RT_models[i_swin] == :XRTM


            # If we are doing Jacobians, we also need arrays to hold
            # weighting functions.

            #=
                How many weighting functions will we need?

                At this point, when the buffer is created, we do not know which
                atmospheric elements are part of the atmosphere, what type of observer
                will be used, and so on. Hence, we must make a guess as to how many
                weighting functions will be used.

                A decent guess for the maximum number that we would ever need can be based
                off the number of layers. If we need per-layer derivatives for optical
                depth and single-scatter albedo, we would need 2xN_lay weighting
                functions.

                We only allocate this array once during the creation of the buffer, so
                there is no notable performance drawback by making this placeholder
                slightly larger than needed.
            =#

            # N_wfunctions = 2 * N_layer + 2 * N_aerosol + 2 * N_surface kernels?
            N_aero_sv = length(filter(is_aerosol_SVE, sv.state_vector_elements))
            N_wfunctions = 2 * N_layer + N_aero_sv * N_layer + 5

            if sv isa RetrievalStateVector
                hires_wfunctions = [RadType(T, swin.N_hires)
                    for i in 1:N_wfunctions]
            elseif sv isa ForwardModelStateVector
                hires_wfunctions = nothing
            end

            this_rt = MonochromaticRTMethod(
                :XRTM,
                Dict[], # Supply an empty dictionary vector for options at first
                earth_scene,
                optical_properties[swin],
                solar_models[swin],
                sv,
                hires_solar,
                hires_radiance,
                hires_jacobians,
                hires_wfunctions,
                Dict{Any, Vector{Int}}(),
                solar_models[swin].irradiance_unit / u"sr",
                solar_scaler
            )



        else
            error("RT model $(RT_models[i_swin]) not implemented!")
            return nothing
        end

        rt[swin] = this_rt

    end

    return EarthAtmosphereBuffer(
        swins, # List of spectral windows
        earth_scene, # Earth scene object
        optical_properties, # Optical property dictionary
        rt, # Radiative transfer dictionary
        rt_buf, # Radiative transfer buffer
        inst_buf # Instrument model buffer
    )

end

"""
Calculates the indices of an `EarthAtmosphereBuffer`/`AbstractRTBuffer` type to assign
forward model output (radiances) from specific spectral windows to a retrieval-wide array
of radiances. The order in which those are concatenated into the RT buffer is determined
by the order in the `buf.spectral_window` vector.

$(TYPEDSIGNATURES)
"""
function calculate_indices!(buf::EarthAtmosphereBuffer)

    # We loop through all spectral windows present in the buffer
    idx_start = 1

    for swin in buf.spectral_window

        # Clear out old data
        empty!(buf.rt_buf.indices[swin])

        # How many elements do we need for this spectral window?
        L = length(buf.rt_buf.dispersion[swin].index)
        @debug "[BUFFER] $(swin) has $(L) elements"

        idx_stop = idx_start+L-1
        # Store the values into the empty vector
        append!(buf.rt_buf.indices[swin], collect(idx_start:idx_stop))

        @debug "[BUFFER] Added indices for $(swin) => $(idx_start):$(idx_stop)"
        # Set the new start index
        idx_start += L


    end

    #=
    @info "Indices: "
    for swin in buf.spectral_window
        @info "$(swin.window_name): $(buf.rt_buf.indices[swin][1:10])..."
    end
    =#
end



"""
$(TYPEDSIGNATURES)

Takes the optical properties and some other scene-relevant data and dumps them into
an HDF5 file. This is intended for e.g. de-bugging so that an external code could take
these (total) optical properties and perform an independent RT calculation that could
be compared to the output generated by RetrievalToolbox.
"""
function dump_to_hdf5(
    buf::EarthAtmosphereBuffer,
    fname::String
    )

    h5open(fname, "w") do h5out

        # Window-independent quantities
        h5out["SZA"] = [buf.scene.solar_zenith]
        h5out["SAA"] = [buf.scene.solar_azimuth]
        h5out["VZA"] = [buf.scene.observer.viewing_zenith]
        h5out["VAA"] = [buf.scene.observer.viewing_azimuth]

        swin_list = [swin.window_name for swin in keys(buf.rt)]
        h5out["windows"] = swin_list
        h5out["n_layer"] = [buf.scene.atmosphere.N_layer]

        # Figure out how many coefs there are at most..
        max_coef = [size(rt.optical_properties.total_coef)[1]
            for (swin, rt) in buf.rt] |> maximum

        h5out["max_coef"] = [max_coef]


        for (swin, rt) in buf.rt

            @info "Exporting $(swin.window_name)"
            h5g = HDF5.create_group(h5out, swin.window_name)

            # Layer optical inputs: τ, ω and β
            h5g["tau"] = rt.optical_properties.total_tau
            h5g["omega"] = rt.optical_properties.total_omega
            h5g["coef"] = rt.optical_properties.total_coef

            # Surface optical inputs:

            # Get the surface first
            surface = get_surface(rt.scene, swin)
            # How many BRDF kernels do we have?
            N_surf = length(surface.kernels)

            surface_array = zeros(swin.N_hires, N_surf)

            kernel_list = String[]

            # Loop over all BRDF kernels in this surface
            for (k, kernel) in enumerate(surface.kernels)

                # Grab the "XRTM"-name for this kernel
                push!(kernel_list, get_XRTM_name(kernel))

                for i in 1:swin.N_hires
                    surface_array[i,k] = evaluate_surface_at_idx(kernel, i)
                end
            end

            # Export the surface
            h5g["surf"] = surface_array
            h5g["kernels"] = kernel_list


        end # End RT/spectral window loop

    end # h5out will be closed here..

end