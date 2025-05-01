"""
$(TYPEDSIGNATURES)

Helper function to populate an **EarthAtmosphereBuffer**, which also includes an
**EarthAtmosphere** and the corresponding **OpticalProperties** with correctly sized
arrays.

# Details

TODO: a lot of documentation needs to go here.

"""
function EarthAtmosphereBuffer(
    sv,
    spectral_windows,
    surface_types,
    atmospheric_elements,
    solar_models,
    RT_models::Vector{Symbol},
    RadType::Type{<:Radiance},
    rt_buf,
    inst_buf,
    N_level,
    N_met_level,
    ::Type{T}
    ) where {T <: AbstractFloat}

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
        u"Pa",
        u"Pa",
        u"K",
        u"kg/kg",
        u"m",
        u"m/s^2"
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

        # TODO delegate this into its own function: create_empty_EarthOpticalProperties

        # Dictionaries to keep track of gas optical depths and derivatives
        gas_tau = Dict{GasAbsorber{T}, Array{T, 2}}()
        gas_derivatives = Dict{GasAbsorber{T}, Dict{String, AbstractArray}}()
        # Dictionaries to keep track of aerosol optical depths and single-scatter albedo
        aerosol_tau = Dict{AbstractAerosolType, Array{T, 2}}()
        aerosol_omega = Dict{AbstractAerosolType, Array{T, 2}}()

        # Rayleigh extinction optical depth
        rayleigh_tau = zeros(T, swin.N_hires, N_layer)
        # Total optical depth and single-scatter albedo
        total_tau = zeros(T, swin.N_hires, N_layer)
        total_omega = zeros(T, swin.N_hires, N_layer)

        for atm in atm_elements

            # Gases
            if atm isa GasAbsorber
                # per-gas layer-resolved optical depth
                gas_tau[atm] = zeros(T, swin.N_hires, N_layer)

                # related partial derivatives w.r.t.
                # surface pressure, temperature, VMRs
                # -> zero-arrays are created later when
                #    looping through the state vector elements

                # Note: this only creates an empty Dict, and is
                # filled just below if needed.
                gas_derivatives[atm] = Dict{String, AbstractArray}()

                if sv isa RetrievalStateVector
                    for _sve in sv.state_vector_elements

                        if _sve isa SurfacePressureSVE
                            gas_derivatives[atm]["dTau_dpsurf"] =
                                zeros(T, swin.N_hires)
                        end

                        if _sve isa GasLevelScalingFactorSVE
                            gas_derivatives[atm]["dTau_dVMR"] =
                                zeros(T, swin.N_hires, N_layer, 2)
                        end

                        if _sve isa GasVMRProfileSVE
                            gas_derivatives[atm]["dTau_dVMR"] =
                                zeros(T, swin.N_hires, N_layer, 2)
                        end

                        if _sve isa TemperatureOffsetSVE
                            gas_derivatives[atm]["dTau_dT"] =
                                zeros(T, swin.N_hires, N_layer)
                        end

                    end
                end

            # End of atm-loop
            end


            # TODO: add other atmospheric elements
            # (aerosols etc.)
            if atm isa AbstractAerosolType

                aerosol_tau[atm] = zeros(T, swin.N_hires, N_layer)
                aerosol_omega[atm] = zeros(T, swin.N_hires, N_layer)

            end

        end

        # Allocate dry and wet air molecule numbers
        nair_dry = zeros(T, N_layer)
        nair_wet = zeros(T, N_layer)

        # Allocate vector(s) to hold temporary values
        tmp_Nhi1 = zeros(T, swin.N_hires)
        tmp_Nhi2 = zeros(T, swin.N_hires)
        tmp_Nlay1 = zeros(T, N_layer)
        tmp_Nlay2 = zeros(T, N_layer)

        # Determine the maximal number of phasefunction expansion
        # coefficients found in all aerosols of this scene.
        max_coef = 1
        for atm in atm_elements
            if atm isa AbstractAerosolType
                if atm.optical_property.max_coefs > max_coef
                    max_coef = atm.optical_property.max_coefs
                end
            end
        end

        # If there are no aerosols, but there is RayleighScattering, adjust here
        if findanytype(atm_elements, RayleighScattering)
            # Bump up to 3 coefficients
            max_coef = max(max_coef, 3)
        end

        # Allocate matrix to hold phasefunction expansion coefficient
        # for this spectral window.
        # NOTE/TODO
        # Currently, we only use one per band (no spectral variation)
        # NOTE/TODO
        # Expand this for polarized RT
        if (length(keys(aerosol_tau)) > 0) | findanytype(atm_elements, RayleighScattering)

            # Check if we are doing polarization or not
            if RadType === ScalarRadiance
                n_elem = 1
            elseif RadType === VectorRadiance
                n_elem = 6
            end

            total_coef = zeros(T, max_coef, n_elem, N_layer)
            tmp_coef = zeros(T, max_coef, n_elem, N_layer)
        else
            # No aerosols, no Rayleigh scattering => no need for coefficients
            total_coef = nothing
            tmp_coef = nothing
        end

        # Allocate optical property object
        this_opt_prop = EarthAtmosphereOpticalProperties(
            swin,
            gas_tau,
            gas_derivatives,
            aerosol_tau,
            aerosol_omega,
            rayleigh_tau,
            total_tau,
            total_omega,
            total_coef,
            nair_dry,
            nair_wet,
            tmp_Nhi1,
            tmp_Nhi2,
            tmp_Nlay1,
            tmp_Nlay2,
            tmp_coef
        )

        # push into dict
        optical_properties[swin] = this_opt_prop

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
            N_wfunctions = 2 * N_layer + 20

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
        @debug "$(swin) has $(L) elements"

        idx_stop = idx_start+L-1
        # Store the values into the empty vector
        append!(buf.rt_buf.indices[swin], collect(idx_start:idx_stop))

        @debug "Added indices for $(swin) => $(idx_start):$(idx_stop)"
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