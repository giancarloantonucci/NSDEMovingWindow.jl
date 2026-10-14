@testset "fixed MoWi: end-to-end against the fine solve" begin
    mowi = MoWi(good_parareal(); τ=1.0, Δτ=0.5)
    solution = solve(problem, mowi)
    @test length(solution) == 5                                # (3−1)/0.5 + 1
    @test all(==(0), solution.restarts)
    @test windowstarts(solution) ≈ 0.0:0.5:2.0                 # boundaries advance by Δτ
    @test NSDEMovingWindow.windowspan(solution[end])[2] ≥ 3.0  # full coverage
    for t in 0.0:0.3:3.0
        @test norm(solution(t) - finetruth(t)) < 1e-6
    end
    # Overlap consistency: neighbouring windows agree on their shared stretch.
    for t in 0.5:0.1:1.0
        @test norm(solution[2](t) - solution[1](t)) < 1e-6
    end
end

@testset "implicit coarse solver through MoWi: success and failure" begin
    # Success: backward Euler as the coarse propagator of every window.
    implicit = Parareal(
        finesolver,
        BackwardEuler(h=5e-2);
        parameters=PararealParameters(N=4, K=4),
        criterion=MeanRelativeUpdate(ϵ=1e-9),
    )
    solution = solve(problem, MoWi(implicit; τ=1.0, Δτ=0.5))
    @test length(solution) == 5
    @test all(==(0), solution.restarts)
    for t in 0.0:0.3:3.0
        @test norm(solution(t) - finetruth(t)) < 1e-6
    end
    # A failed coarse stage gives a nonfinite window. Fixed MoWi stops and
    # retains that window; adaptive strategies may retry it.
    blow = IVP((u, t) -> u .^ 2, [1.0], (0.0, 0.9))
    failing = Parareal(
        RK4(h=1e-3),
        BackwardEuler(h=0.3);
        parameters=PararealParameters(N=3, K=3),
        criterion=MeanRelativeUpdate(ϵ=1e-9),
    )
    failed = @test_logs (:warn, r"Chunk solve failed") (:warn, r"stopped at a nonfinite") match_mode=:any solve(
        blow, MoWi(failing; τ=0.9, Δτ=0.3)
    )
    @test failed.status == :nonfinite
end

@testset "complex states through MoWi" begin
    # KS pushes complex spectral states through the whole MoWi stack in the
    # application; CI pins the support here on a complex Dahlquist.
    λ = -0.3 + 1.5im
    zproblem = Dahlquist([1.0 + 0.0im], (0.0, 3.0); λ)
    zfine = solve(zproblem, finesolver)
    zparareal = Parareal(
        finesolver,
        RungeKutta4(h=5e-2);
        parameters=PararealParameters(N=4, K=4),
        criterion=MeanRelativeUpdate(ϵ=1e-10),
    )
    zsolution = solve(zproblem, MoWi(zparareal; τ=1.0, Δτ=0.5))
    u, _ = flatten(zsolution)
    @test eltype(first(u)) == ComplexF64                 # nothing silently realified
    for t in 0.0:0.4:3.0
        @test norm(zsolution(t) - zfine(t)) < 1e-7       # window hand-offs preserve ℂ accuracy
    end
end

@testset "strategies on a chaotic problem (smoke)" begin
    # Stretch/Leap/Zoom are pinned mechanically on Logistic with bad_parareal
    # above; this runs all three END-TO-END on Lorenz — the chaotic regime
    # actually deploys them in — at a CI-sized span.
    lproblem = Lorenz([2.0, 3.0, -14.0], (0.0, 6.0))
    ltruth = solve(lproblem, finesolver)
    lparareal(criterion) = Parareal(
        finesolver, RungeKutta4(h=5e-2); parameters=PararealParameters(N=4, K=4), criterion
    )
    runs = (
        (
            "Stretch",
            MoWi(
                lparareal(MeanRelativeUpdate(ϵ=1e-10)); τ=2.0, Δτ=1.0, adaptive=StretchParameters()
            ),
        ),
        (
            "Leap",
            MoWi(lparareal(MeanRelativeUpdate(ϵ=1e-10)); τ=2.0, Δτ=1.0, adaptive=LeapParameters()),
        ),
        # Zoom needs ψ₂ with a FROZEN base weight (the constructor refuses
        # adaptive weights); w = exp(0.9056), an estimated Lorenz Lyapunov rate.
        (
            "Zoom",
            MoWi(
                lparareal(WeightedDefect(ϵ=1e-10, weight=FixedWeight(exp(0.9056))));
                τ=2.0,
                Δτ=1.0,
                adaptive=ZoomParameters(),
            ),
        ),
    )
    for (name, mowi) in runs
        solution = solve(lproblem, mowi)
        @test NSDEMovingWindow.windowspan(solution[end])[2] ≥ 6.0  # full coverage
        @test all(≤(mowi.adaptive.R), solution.restarts)           # budget respected
        for t in 0.0:1.5:6.0
            # Pointwise tracking on a chaotic system, honestly bounded. The
            # per-window residual (∼ϵ scale) is amplified by LOCAL expansion
            # rates — Lorenz's ‖J‖ bursts far exceed the mean rate Λ ≈ 0.9,
            # so a mean-rate estimate undershoots by orders of magnitude
            # (measured: up to ∼4e-2 absolute by t = 6 at ϵ = 1e-10). The
            # scaled bound below asserts "still on the TRUE trajectory to
            # ∼5%" — which uncontrolled divergence, O(attractor diameter),
            # fails by more than an order of magnitude — with ×10-plus
            # headroom over the measured worst for cross-platform rounding.
            @test norm(solution(t) - ltruth(t)) < 5e-2 * (1 + norm(ltruth(t)))
        end
    end
end

@testset "kwargs reach the subroutine" begin
    mowi = MoWi(good_parareal(); τ=1.5, Δτ=1.5)
    solution = solve(problem, mowi; history=:iterates)
    @test solution[1].iterates !== nothing                     # the old dispatcher dropped this
end

@testset "Window settings are reset between solves" begin
    p = Logistic(0.3, (0.0, 2.0))
    for adaptive in (nothing, StretchParameters(), LeapParameters(), ZoomParameters())
        pr = Parareal(
            RK4(h=0.01),
            Euler(h=0.1);
            N=4,
            criterion=WeightedDefect(ϵ=1e-8, weight=FixedWeight(2.0)),
        )
        m = MoWi(pr; τ=0.5, Δτ=0.25, adaptive)
        a = solve(p, m)
        b = solve(p, m)
        @test length(a) == length(b)
        @test [x.diagnostics.status for x in a.windows] == [x.diagnostics.status for x in b.windows]
        @test [boundaryvalues(x.lastiterate) for x in a.windows] == [boundaryvalues(x.lastiterate) for x in b.windows]
        @test m.parallelsolver.criterion.weight.w == pr.criterion.weight.w == 2.0
        @test m.parallelsolver.finesolver.stepsize.h == 0.01
        @test m.parallelsolver.coarsesolver.stepsize.h == 0.1
    end
end
