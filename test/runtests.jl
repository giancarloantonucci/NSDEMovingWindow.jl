using NSDEMovingWindow
using NSDERungeKutta
using LinearAlgebra
using Test
using Aqua
using RecipesBase: RecipesBase
using NSDEBase: NSDEBase

include("support.jl")

@testset "NSDEMovingWindow" begin
    @testset "Aqua" begin
        Aqua.test_all(
            NSDEMovingWindow; piracies=(; treat_as_own=[NSDEMovingWindow.RecipesBase.recipetype])
        )
    end

    include("construction.jl")
    include("integration.jl")
    include("strategies.jl")
    include("solutions_and_statistics.jl")
    include("reuse_and_failures.jl")
    include("reliability.jl")
end
