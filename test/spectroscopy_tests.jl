@testset "Spectroscopy utils" begin

    ##
    ww1 = rand(100); sort!(ww1);
    ww2 = rand(10); sort!(ww2);
    result = searchsortedfirst.(Ref(ww1), ww2) .- 1
    @test all(RetrievalToolbox._find_ww_indices(ww1, ww2) .== result)

end