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