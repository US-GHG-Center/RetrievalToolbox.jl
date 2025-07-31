using RetrievalToolbox
using Test

@testset "RetrievalToolbox tests" begin

    @testset "Utils tests" begin
        include("utils_tests.jl")
    end
    @testset "Spectroscopy tests" begin
        include("spectroscopy_tests.jl")
    end
    @testset "Atmosphere tests" begin
        include("atmosphere_tests.jl")
    end

end