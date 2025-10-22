@testset "Example atmospheres" begin

    atmos_names = list_example_atmospheres()

    # Check if we have at least one example atmosphere
    @test length(atmos_names) > 0

    # We will do the following tests for ALL atmospheres!
    for atmos in atmos_names

        # Create empty atmosphere with default unit types
        atm = create_example_atmosphere(atmos, 10)
        # Make sure the type is right
        @test atm isa EarthAtmosphere
        # Check gravity calculation (must monotonically decrease with leve)
        @test all(diff(atm.gravity_levels) .> 0)
        # Check altitude levels (must monotonically increase with level)
        @test all(diff(atm.altitude_levels) .< 0)

    end

end


@testset "Gravity" begin

    # See if gravity functions work

    gn = JPL_gravity(0.0, 0.0u"m")
    # Check unit
    @test gn isa Unitful.Acceleration
    # Check that the result is what we expect it to be (roughly)
    @test abs(gn - 9.7803267715u"m/s^2" |> u"m/s^2" |> ustrip) < 1e-3
    # Gravity needs to be smaller poleward
    gn_pole = JPL_gravity(89.0, 0.0u"m")
    @test gn_pole < gn
    # Gravity needs to be smaller when higher up
    gn_higher = JPL_gravity(0.0, 4000.0u"m")
    @test gn > gn_higher

end


@testset "Pressure schemes" begin

    plevels = create_UoL_pressure_grid(1000.0u"hPa", 123.0u"hPa")
    # Make sure they are in correct order
    @test all(diff(plevels) .> 0u"Pa")

    plevels = create_ACOS_pressure_grid(1001.0u"hPa")
    # Make sure they are in correct order
    @test all(diff(plevels) .> 0u"Pa")

end

@testset "Atmosphere object" begin

    # Try creating an empty atmosphere
    atm = create_empty_EarthAtmosphere(6, 11, Float64)

    # Fill met p with values
    met_p = [1., 10., 50., 100., 300., 400., 500., 600., 800., 900., 1000.] * u"hPa"
    ingest!(atm, :met_pressure_levels, met_p)

    # Calculate z from p using simple barometric formula
    # (discouraged otherwise!)
    met_z = @. 44330.0u"m" * (1.0 - met_p / met_p[end])
    ingest!(atm, :altitude_levels, met_z)

    # Calculate g from z!
    calculate_gravity_from_z!(atm)

    # Add some very simple temperature profile
    met_t = LinRange(210.0u"K", 250u"K", atm.N_met_level)
    ingest!(atm, :temperature_levels, met_t)

end


@testset "SIF" begin

    # Try creating a SIF object
    sif = SIFRadiance(
        1.0e-7,
        u"W/m^2/sr/µm",
        750.0,
        u"nm"
    )

    # Get some SIF radiance value at arbitrary wavelength
    get_SIF_radiance(
        sif,
        0.755u"µm"
    )

    # Get some SIF radiance value at arbitrary wavenumber
    get_SIF_radiance(
        sif,
        13000.0u"cm^-1"
    )

    # Make sure SIF is zero far away from expected emission range
    @test get_SIF_radiance(sif, 0.01u"µm") ≈ 0.
    @test get_SIF_radiance(sif, 10.0u"µm") ≈ 0.
end