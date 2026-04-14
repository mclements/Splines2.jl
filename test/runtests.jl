include("set_up_tests.jl")

@testset "Aqua" begin
    Aqua.test_all(Splines2; ambiguities=false)
end

@testset "bspline" include("bspline.jl")
@testset "nspline" include("nspline.jl")
@testset "ispline" include("ispline.jl")
@testset "mspline" include("mspline.jl")
@testset "StatsModels" include("statsmodels.jl")
