@testset "Forward run MWE 1" begin


    wavelength_start = 761.0u"nm"
    wavelength_end = 762.0u"nm"
    gas_name = "O2"
    N_level = 10

    my_source_atmosphere = create_example_atmosphere("US-midwest-summer", N_level)

    if isfile("demo_spectroscopy.jld2")

        my_spectroscopy = jldopen("demo_spectroscopy.jld2")["my_spectroscopy"]

    else

        my_spectroscopy = create_ABSCO_from_HITRAN(
            gas_name,
            wavelength_start - 0.5u"nm", # put some buffer around our requested wavelengths
            wavelength_end + 0.5u"nm";
            wavelength_output=true,
            h2o_broadening=false
        )

        jldsave("demo_spectroscopy.jld2"; my_spectroscopy)
    end

    my_gas = GasAbsorber(
        gas_name,
        my_spectroscopy,
        fill(20.95, N_level), # Oxygen VMR is ~20.95%
        Unitful.percent # Unit of VMR profile
    )

    my_spectral_window = spectralwindow_from_ABSCO(
        "$(gas_name)_window",
        wavelength_start |> ustrip,
        wavelength_end |> ustrip,
        (wavelength_start + wavelength_end) / 2 |> ustrip,
        0.2, # Buffer length in same units as everything else
        my_spectroscopy,
        unit(wavelength_start)
    )

    my_dispersion = SimplePolynomialDispersion(
        [wavelength_start, 0.01u"nm"],
        1:100,
        my_spectral_window
    )

    my_state_vector = ForwardModelStateVector();

    N1 = 7_000
    N2 = 150

    my_instrument_buffer = InstrumentBuffer(
        zeros(Float64, N1),
        zeros(Float64, N1),
        zeros(Float64, N2),
    )

    my_rt_buffer = ScalarRTBuffer(
        Dict(my_spectral_window => my_dispersion), # Already a SpectralWindows -> Dispersion dictionary
        ScalarRadiance(Float64, N2), # Hold the radiance - we use ScalarRadiance because we don't need polarization
        nothing,
        Dict(my_spectral_window => zeros(Int, 0)), # Hold the detector indices
        Unitful.NoUnits # Radiance units for the forward model
    )

    my_buffer = EarthAtmosphereBuffer(
        my_state_vector,
        my_spectral_window,
        [(:Lambert, 1)],
        [my_gas],
        Dict(my_spectral_window => UnitSolarModel()),
        [:BeerLambert],
        ScalarRadiance, # Use ScalarRadiance for high-res radiance calculations
        my_rt_buffer,
        my_instrument_buffer,
        N_level, # The number of retrieval or RT pressure levels
        my_source_atmosphere.N_met_level, # The number of meteorological pressure levels, as given by the atmospheric inputs
        Float64 # The chosen Float data type (e.g. Float16, Float32, Float64)
    )

    calculate_indices!(my_buffer)

    my_buffer.scene.atmosphere.met_pressure_levels[:] =
        my_source_atmosphere.met_pressure_levels[:]
    my_buffer.scene.atmosphere.specific_humidity_levels[:] =
        my_source_atmosphere.specific_humidity_levels[:]
    my_buffer.scene.atmosphere.temperature_levels[:] =
        my_source_atmosphere.temperature_levels[:]

    my_buffer.scene.atmosphere.pressure_levels[:] = logrange(10, 1000_00, N_level) |> collect

    calculate_altitude_and_gravity!(my_buffer.scene)

    my_buffer.scene.surfaces[my_spectral_window].coefficients[1] = 0.25

    calculate_earth_optical_properties!(
        my_buffer.rt[my_spectral_window].optical_properties,
        my_buffer.scene,
        my_state_vector
    )

    # Some optical properties must be > 0
    total_tau = sum(my_buffer.rt[my_spectral_window].optical_properties.total_tau, dims=2)
    @test any(total_tau .> 0)

    calculate_solar_irradiance!(
        my_buffer.rt[my_spectral_window],
        my_spectral_window,
        my_buffer.rt[my_spectral_window].solar_model
    )

    calculate_radiances_and_jacobians!(my_buffer.rt[my_spectral_window])

    # Can't have all zeros! RT must produce some values > 0
    @test any(my_buffer.rt[my_spectral_window].hires_radiance .> 0)

    my_isrf = GaussISRF(0.05, u"nm")

    isrf_success = apply_isrf_to_spectrum!(
        my_instrument_buffer,
        my_isrf,
        my_dispersion,
        my_buffer.rt[my_spectral_window].hires_radiance.I
    )

    @test isrf_success == true

end