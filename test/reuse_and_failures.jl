@testset "MoWi reuse resets seeding and restart counters" begin
    p = Logistic(0.3, (0, 3))
    solver = Parareal(RK4(h=0.01), Euler(h=0.1); N=4, criterion=NoStop())
    mowi = MoWi(solver; τ=1.0, Δτ=0.5)
    cache = NSDEMovingWindow.MoWiCache(p, mowi)
    result = MoWiSolution(p, mowi)
    mowi(cache, result, p)
    q = remake(p; u0=[0.7])
    fill!(result.restarts, 99)
    cache.windowcache.makeGs .= false
    mowi(cache, result, q)
    fresh = solve(q, mowi)
    @test result.status == :converged
    @test result.restarts == fresh.restarts == zeros(Int, length(result))
    @test [boundaryvalues(w.lastiterate) for w in result.windows] == [boundaryvalues(w.lastiterate) for w in fresh.windows]
    @test result[1][1].u[1] == [0.7]
    @test p.u0 == [0.3]
    @test q.u0 == [0.7]
    @test result[2][1].u[1] != q.u0
    @test_throws DimensionMismatch mowi(cache, result, remake(q; u0=[0.7, 0.8]))
    @test_throws ArgumentError mowi(cache, result, q; nocollect=true)
end

@testset "Unrecoverable nonfinite window stops the run" begin
    p = Dahlquist([1.0], (0.0, 32.0); λ=-1000.0)
    solver = Parareal(RK4(h=0.001), Euler(h=0.02); N=8)
    for strategy in (nothing, ZoomParameters(R=2))
        mowi = MoWi(remake(solver; criterion=WeightedDefect()); τ=8.0, Δτ=4.0, adaptive=strategy)
        failed = @test_logs (:warn, r"stopped at a nonfinite") solve(p, mowi)
        @test failed.status == :nonfinite
        @test length(failed) == 1
        @test failed.restarts == [strategy === nothing ? 0 : 2]
        @test failed[1].diagnostics.status == :nonfinite
        @test_throws ArgumentError ensemblemean(failed)
        @test_throws ArgumentError ensemblevariance(failed)
        @test_throws ArgumentError ensemblesem(failed)
    end
end

@testset "Zoom remains responsive after many successes" begin
    p = Parareal(RK4(h=0.01), Euler(h=0.1); N=4, criterion=WeightedDefect(weight=FixedWeight(2.0)))
    for _ in 1:2000
        p = NSDEMovingWindow.zoom_solver(p, 0.25)
    end
    @test p.criterion.weight.w == 1.0
    @test NSDEMovingWindow.zoom_solver(p, 2.0).criterion.weight.w == 2.0
end

@testset "Adaptive strategy inputs" begin
    for bad in (0.0, -1.0, 1.0, 2.0, NaN, Inf)
        @test_throws ArgumentError StretchParameters(δτ⁻=bad)
        @test_throws ArgumentError LeapParameters(δΔτ⁻=bad)
    end
    for bad in (0.0, -1.0, 0.5, NaN, Inf)
        @test_throws ArgumentError StretchParameters(δτ⁺=bad)
        @test_throws ArgumentError LeapParameters(δΔτ⁺=bad)
    end
    for bad in (0.0, -1.0, 1.0, NaN, Inf)
        @test_throws ArgumentError ZoomParameters(δw⁻=bad)
    end
    for bad in (0.0, -1.0, 2.0, NaN, Inf)
        @test_throws ArgumentError ZoomParameters(δw⁺=bad)
    end
    @test_throws ArgumentError StretchParameters(R=-1)
    @test_throws ArgumentError LeapParameters(R=-1)
    @test_throws ArgumentError ZoomParameters(R=-1)
    @test_throws ArgumentError StretchParameters(2.0, 1.0, 1, false)
    @test_throws ArgumentError LeapParameters(2.0, 1.0, 1)
    @test_throws ArgumentError ZoomParameters(0.5, 0.25, 1)
end

@testset "MoWi accepts a changed RHS of the same type" begin
    p = IVP(-1.0, [1.0], (0.0, 3.0))
    q = IVP(-2.0, [0.5], (0.0, 3.0))
    solver = MoWi(Parareal(RK4(h=0.01), Euler(h=0.1); N=4, criterion=NoStop()); τ=1.0, Δτ=0.5)
    cache = NSDEMovingWindow.MoWiCache(p, solver)
    result = MoWiSolution(p, solver)
    solver(cache, result, p)
    solver(cache, result, q)
    fresh = solve(q, solver)
    @test [boundaryvalues(w.lastiterate) for w in result.windows] == [boundaryvalues(w.lastiterate) for w in fresh.windows]
    @test cache.windowcache.chunkproblems[1].rhs === q.rhs
    @test_throws ArgumentError solver(cache, result, IVP((u, t)->-u, [0.5], (0.0, 3.0)))
end
