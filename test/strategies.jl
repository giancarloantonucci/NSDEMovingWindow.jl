@testset "Stretch (AMoWi-1): restarts fire, τ shrinks, R caps" begin
    mowi = MoWi(bad_parareal(); τ=1.0, Δτ=0.5, adaptive=StretchParameters(R=3))
    solution = @test_logs (:warn, r"accepted WITHOUT convergence") match_mode=:any solve(
        problem, mowi
    )
    @test any(>(0), solution.restarts)
    @test all(≤(3), solution.restarts)                         # R caps every window
    @test length(unique(round.(windowlengths(solution); digits=10))) > 1 # τ demonstrably moved
    @test NSDEMovingWindow.windowspan(solution[end])[2] ≥ 3.0  # still covers the span
end

@testset "Leap (AMoWi-2): hooks, exemption, cap" begin
    # Unit level (a forced-failure end-to-end can shrink-spiral by design, so
    # the restart mechanics are pinned here deterministically). Floors are set
    # NON-BINDING (0.0) on purpose: this testset pins the hop/leap arithmetic;
    # the floor behaviour has its own testset below.
    state = NSDEMovingWindow.MoWiState(1.0, 0.5, 0.0, 0.0)
    NSDEMovingWindow.adapt_failure!(LeapParameters(), state, nothing)
    @test state.Δτ ≈ 0.125                                     # hop: ×δΔτ⁻
    NSDEMovingWindow.adapt_success!(LeapParameters(δΔτ⁺=100.0), state, nothing)
    @test state.Δτ ≈ 1.0                                       # leap, capped at τ
    @test first(NSDEMovingWindow.restart_start(LeapParameters(), 7.0, 2.0, 0.25, 0.0, 0.0)) ≈ 2.25 # retry the SHIFT
    @test first(NSDEMovingWindow.restart_start(StretchParameters(), 7.0, 2.0, 0.25, 0.0, 0.0)) ≈ 7.0 # same window start
    @test !NSDEMovingWindow.allowsrestart(LeapParameters(), 1)  # window 1 exempt
    @test NSDEMovingWindow.allowsrestart(LeapParameters(), 2)
    @test NSDEMovingWindow.allowsrestart(StretchParameters(), 1)
    @test !NSDEMovingWindow.allowsrestart(nothing, 2)           # fixed MoWi never restarts
    # Benign setup with aggressive growth: the leap saturates at τ.
    mowi = MoWi(good_parareal(); τ=1.0, Δτ=0.25, adaptive=LeapParameters(δΔτ⁺=4.0))
    solution = solve(problem, mowi)
    starts = windowstarts(solution)
    @test length(starts) ≥ 3
    @test starts[3] - starts[2] ≈ 1.0 atol = 1e-12             # capped shift ≡ τ
end

@testset "failure hooks floor their state (termination guards)" begin
    # Without floors, Δτ (Leap) and τ (Stretch) decay geometrically across
    # consecutively failing windows: the advance per window tends to zero and
    # mowi! never terminates. The floors are guards, not adaptive policy.
    parareal = good_parareal()
    N = parareal.parameters.N
    mkstate(parareal) = NSDEMovingWindow.MoWiState(
        1.0, 0.5, 1.0 / N, max(N * parareal.finesolver.stepsize.h, 1.0 / N)
    )
    st = mkstate(parareal)
    for _ in 1:100
        NSDEMovingWindow.adapt_failure!(LeapParameters(δΔτ⁻=0.1), st, parareal)
    end
    @test st.Δτ ≥ st.Δτmin             # floored at one ORIGINAL chunk length
    parareal2 = good_parareal()
    st = mkstate(parareal2)
    for _ in 1:100
        NSDEMovingWindow.adapt_failure!(StretchParameters(δτ⁻=0.1), st, parareal2)
    end
    @test st.τ ≥ st.τmin               # window floor holds
    @test st.Δτ ≥ st.Δτmin             # ABSOLUTE shift floor: τ₀/N, not the
    # shrunken τ/N — otherwise the advance decays to one fine step per window
end

@testset "Leap: terminates under persistent failure" begin
    parareal = bad_parareal()
    N = parareal.parameters.N
    tN = problem.tspan[2]
    mowi = MoWi(parareal; τ=1.0, Δτ=0.5, adaptive=LeapParameters(δΔτ⁻=0.1, R=2))
    solution = @test_logs (:warn, r"accepted WITHOUT convergence") match_mode=:any solve(
        problem, mowi
    )
    # The termination guard, in numbers: every ACCEPTED shift sits between the
    # ABSOLUTE floor τ₀/N and the overlap cap τ, so the window count is bounded.
    @test all(x -> 1.0 / N - 1e-12 ≤ x ≤ 1.0 + 1e-12, diff(windowstarts(solution)))
    @test length(solution) ≤ ceil(Int, (tN - 1.0) / (1.0 / N)) + 2
    @test solution.restarts[1] == 0    # window 1 exempt
    @test any(>(0), solution.restarts) # the machinery fired
    @test all(≤(2), solution.restarts) # R caps every window
    # NOT asserted: restarts == R everywhere. A reduced shift means more
    # overlap and better seeding, so a retry can CONVERGE before exhausting R
    # — that is recovery working, not a failed guard. bad_parareal's ϵ is
    # deliberately reachable (see its comment in test/support.jl).
end

@testset "Zoom (AMoWi-3): multiplicative, per restart, on the base weight" begin
    # Unit-level direction check on the hooks:
    parareal = MoWi(
        bad_parareal(criterion=WeightedDefect(ϵ=1e-5, weight=FixedWeight(2.0)));
        τ=1.0,
        adaptive=ZoomParameters(),
    ).parallelsolver
    state = NSDEMovingWindow.MoWiState(1.0, 1.0, 0.0, 0.0) # Zoom's hooks never touch the state
    parareal = NSDEMovingWindow.adapt_failure!(ZoomParameters(), state, parareal)
    @test parareal.criterion.weight.w ≈ 4.0                   # ×δw⁻ = ×2: looser
    parareal = NSDEMovingWindow.adapt_success!(ZoomParameters(), state, parareal)
    @test parareal.criterion.weight.w ≈ 1.0                   # ×δw⁺ = ×0.25: tighter
    # End-to-end: restarts fire; τ and Δτ never move; the weight has been driven.
    mowi = MoWi(
        bad_parareal(criterion=WeightedDefect(ϵ=1e-5, weight=FixedWeight(2.0)));
        τ=1.0,
        Δτ=0.5,
        adaptive=ZoomParameters(R=3),
    )
    solution = solve(problem, mowi)
    @test any(>(0), solution.restarts)
    @test all(l -> isapprox(l, 1.0; atol=1e-12), windowlengths(solution))
    @test diff(windowstarts(solution)) ≈ fill(0.5, length(solution) - 1) atol = 1e-12
    @test mowi.parallelsolver.criterion.weight.w == 2.0
    @test any(win.diagnostics.weights[end] != 2.0 for win in solution.windows)
end

@testset "nonfinite coarse starts restart without reading empty errors" begin
    p = Dahlquist([1.0], (0.0, 8.0); λ=-1000.0)
    parareal = Parareal(RK4(h=0.001), Euler(h=0.02); N=8)
    initial = solve(p, parareal)
    @test initial.diagnostics.status == :nonfinite
    @test isempty(initial.errors)
    solver = MoWi(parareal; τ=8.0, Δτ=4.0, adaptive=StretchParameters(δτ⁺=1.0, R=6))
    recovered = solve(p, solver)
    @test recovered.restarts[1] > 0
    @test all(win -> win.diagnostics.status in (:converged, :exact), recovered.windows)
    @test all(t -> all(isfinite, recovered(t)), 0.0:0.25:8.0)
    # Zoom cannot cure a nonfinite coarse map, but still exhausts its restart
    # budget cleanly and returns explicit NaNs rather than an unrelated BoundsError.
    zoom = MoWi(
        remake(parareal; criterion=WeightedDefect()); τ=8.0, Δτ=4.0, adaptive=ZoomParameters(R=2)
    )
    failed = @test_logs (:warn, r"stopped at a nonfinite") solve(p, zoom)
    @test failed.restarts == [2]
    @test failed[1].diagnostics.status == :nonfinite
    @test all(isnan, failed(4.0))
end

@testset "Zoom reuses unchanged budget parameters without repeated warnings" begin
    p = @test_logs (:warn, r"exceeds") Parareal(
        RK4(h=0.01), Euler(h=0.1); N=4, K=6, criterion=WeightedDefect()
    )
    state = NSDEMovingWindow.MoWiState(1.0, 1.0, 0.25, 0.25)
    @test_logs for _ in 1:5
        next = NSDEMovingWindow.adapt_failure!(ZoomParameters(), state, p)
        @test next.parameters === p.parameters
        p = NSDEMovingWindow.adapt_success!(ZoomParameters(), state, next)
    end
end
