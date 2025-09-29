using RetrievalToolbox
using Documenter
using Test
using Unitful

@testset "RetrievalToolbox tests" begin

    @testset "Doctests" begin

        # Make sure module is loaded for all doctests
        DocMeta.setdocmeta!(
            RetrievalToolbox,
            :DocTestSetup,
            :(using RetrievalToolbox);
            recursive=true
        )

        doctest(RetrievalToolbox; manual=false)
    end

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