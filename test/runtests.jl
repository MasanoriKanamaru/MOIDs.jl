using MOIDs
using Test

@testset "MOIDs.jl" begin

    ## Comparison with the example at https://github.com/mkretlow/MOID.jl
    ceres = [2.7691652, 0.0760091,  73.59764,  80.30553,10.59407]
    urania = [2.3655722, 0.127581 ,  87.42605, 307.46872, 2.09575]
    @test moid(ceres...,urania...) ≈ 0.24521440655831864
end
