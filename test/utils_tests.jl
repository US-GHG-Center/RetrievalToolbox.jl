@testset "Misc utils" begin

    ##
    @test maxmin([1,2,4,5,8,9]) == 8

    ##
    check_set = Any[Int32(1), 2.0f0, 3im, -99]
    @test findanytype(check_set, Float64) == false
    @test findanytype(check_set, Float32) == true
    @test findany(check_set, -99) == true
    # This evaluates to true!! (-99.0 == -99 .. wow!)
    @test findany(check_set, -99.0) == true
    @test findany(check_set, -10) == false

    ##
    x1 = rand(100)
    x2 = rand(100)
    x1[23] = NaN
    @test check_for_not_finite(x1) == true
    @test check_for_not_finite(x2) == false

    ## AVX functions!
    nrow = 10
    ncol = 8

    # Run these tests for different data types
    for T in [Int32, Int64, Float16, Float32, Float64]

        A = zeros(T, nrow, ncol)
        B = zeros(T, nrow)

        for i in axes(A,1)
            for j in axes(A,2)
                A[i,j] = i * j
            end
        end

        @test avx_sum(A) == sum(1:nrow) * sum(1:ncol)

        avx_sum_along_columns!(B, A)
        @test B == [sum(1:ncol) * one(T) * i for i in 1:nrow]

        B[:] .= 0
        avx_sum_along_columns_between!(B, A, 3, 5)
        @test B == [sum(3:5) * one(T) * i for i in 1:nrow]

    end

    ## Vector contained in vector!
    for T in [Int16, Int32, Int, Float16, Float32, Float64]
        vecA = convert.(Ref(T), [1,2,3,4,5,6,7,8,9])
        vecB = convert.(Ref(T), [3,4,5,6])
        vecC = convert.(Ref(T), [9,10,11,12])

        @test vector_contained_in_vector(vecA, vecB) == 3
        @test vector_contained_in_vector(vecA, vecC) == -1
    end


    ## Levels to layers
    lev_end = 10
    levels = collect(0:1.0:lev_end)
    layers = collect(0.5:1.0:lev_end-0.5)
    layers_new = zeros(lev_end)

    @test levels_to_layers(levels) == layers
    levels_to_layers!(layers_new, levels)
    @test layers_new == layers


end