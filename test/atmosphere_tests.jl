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
    @test abs(gn - 9.7803267715u"m/s^2" |> u"m/s^2" |> ustrip) < 1e-3

    gn_pole = JPL_gravity(89.0, 0.0u"m")
    @test gn_pole < gn

    gn_higher = JPL_gravity(0.0, 4000.0u"m")
    @test gn > gn_higher

end