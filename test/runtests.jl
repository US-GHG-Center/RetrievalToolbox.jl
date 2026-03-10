using RetrievalToolbox
using Documenter
using JLD2
using Test
using Unitful

@testset "RetrievalToolbox tests" begin

    @testset "Doctests" begin

        # Make sure module is loaded for all doctests
        DocMeta.setdocmeta!(
            RetrievalToolbox,
            :DocTestSetup,
            :(using RetrievalToolbox; using Unitful);
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
    @testset "Compound tests" begin
        include("compound_tests.jl")
    end

end