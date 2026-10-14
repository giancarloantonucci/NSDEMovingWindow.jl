@testset "construction validation" begin
    parareal = good_parareal()
    @test_throws ArgumentError MoWi(parareal; τ=1.0, Δτ=2.0)   # Δτ > τ
    @test_throws ArgumentError MoWi(parareal; τ=1.0, Δτ=0.0)   # no forward shift
    @test_throws ArgumentError MoWi(parareal; τ=0.0)           # no window
    # Zoom over the unweighted ψ₁ would be a silent no-op: refuse loudly.
    @test_throws ArgumentError MoWi(parareal; τ=1.0, adaptive=ZoomParameters())
    weighted = bad_parareal(criterion=WeightedDefect(ϵ=1e-5, weight=FixedWeight(2.0)))
    @test MoWi(weighted; τ=1.0, adaptive=ZoomParameters()) isa MoWi
    # Zoom requires a fixed weight: within-window updates would overwrite it
    # on the next sweep — refuse the pairing loudly at construction.
    remeasured = bad_parareal(criterion=WeightedDefect(ϵ=1e-5, weight=NeighbourWeight()))
    @test_throws ArgumentError MoWi(remeasured; τ=1.0, adaptive=ZoomParameters())
    @test MoWi(remeasured; τ=1.0) isa MoWi                          # fine without Zoom
    @test MoWi(remeasured; τ=1.0, adaptive=LeapParameters()) isa MoWi # Leap never touches w
    # The solver is deep-copied: adaptive runs must not mutate the user's.
    mowi = MoWi(parareal; τ=1.0)
    @test mowi.parallelsolver !== parareal
end
