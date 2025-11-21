@testset "Spectroscopy utils" begin

    ww1 = rand(100); sort!(ww1)
    ww2 = rand(10); sort!(ww2)
    ## Remove all elements of ww2 that are lower than the smallest element of ww1
    ww2 = ww2[ww2 .> ww1[1]]
    ## Remove all elements of ww2 that are larger than the largest element of ww1
    ww2 = ww2[ww2 .< ww1[end]]

    result = searchsortedfirst.(Ref(ww1), ww2) .- 1
    @test all(RetrievalToolbox._find_ww_indices(ww1, ww2) .== result)

end