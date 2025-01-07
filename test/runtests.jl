using MOIDs
using Test

@testset "MOIDs.jl" begin
        
    elements = [
        [" 1 Ceres"     , 2.7691652, 0.0760091,  73.59764,  80.30553, 10.59407],
        ["29 Amphitrite", 2.5541136, 0.0726956,  63.36319, 356.34176,  6.08252],
        ["30 Urania"    , 2.3655722, 0.127581 ,  87.42605, 307.46872,  2.09575],
        ["50 Virginia"  , 2.6487939, 0.2859856, 200.08054, 173.52874,  2.83822],
        ["51 Nemausa"   , 2.3658354, 0.0675594,   2.58053, 175.9785 ,  9.97718],
    ]
    
    ##= Expected results =##
    # Asteroid-1 : Asteroid-2    : MOID [AU]
    # 1 Ceres    : 1 Ceres         0.0
    # 1 Ceres    : 29 Amphitrite   0.15677463452737
    # 1 Ceres    : 30 Urania       0.24521440655832
    # 1 Ceres    : 50 Virginia     0.08934734026105
    # 1 Ceres    : 51 Nemausa      0.35972678460706
    
    result = [0.00000000000001, 0.15677463452737, 0.24521440655832, 0.08934734026105, 0.35972678460706]
    
    println()
    println(rpad("Asteroid-1", 15), rpad("Asteroid-2", 15), rpad("MOID (expected) [au]", 24), rpad("MOID (calucalted) [au]", 24), rpad("Relative error", 24))
    
    for i in eachindex(elements)
        tmp = moid(elements[1][2:6]..., elements[i][2:6]...)
        relative_error = abs(tmp - result[i]) / result[i]
        println(rpad(elements[1][1], 14), rpad(elements[i][1], 14), rpad(result[i], 24), rpad(tmp, 24), rpad(relative_error, 24))
        @test isapprox(tmp, result[i], atol=1e-5)
    end
    
    println()
end
