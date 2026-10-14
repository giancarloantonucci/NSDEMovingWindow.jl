@testset "Complex ensemble observables and sample SEM" begin
    p = IVP(-1im, [1.0], (0.0, 2.0))
    solver = MoWi(Parareal(RK4(h=0.01), Euler(h=0.05); N=2, criterion=NoStop()); τ=0.6, Δτ=0.3)
    sol = solve(p, solver)
    X = NSDEMovingWindow.ensemblesamples(sol, first, 7)
    mean = vec(sum(X; dims=2)) / size(X, 2)
    spread = vec(sum(abs2, X .- mean; dims=2)) / size(X, 2)
    @test eltype(ensemblemean(sol; samples=7)) <: Complex
    @test ensemblemean(sol; samples=7) ≈ mean
    @test ensemblevariance(sol; samples=7) ≈ spread
    @test ensemblesem(sol; samples=7) ≈ sqrt.(spread / (size(X, 2) - 1))
    evaluations = Ref(0)
    ensemblevariance(sol, u -> (evaluations[] += 1; first(u)); samples=7)
    @test evaluations[] == 7count(!, sol.clipped)
    one_window = solve(IVP(-1im, [1.0], (0.0, 0.3)), solver)
    @test_throws ArgumentError ensemblesem(one_window)
    @test_throws ArgumentError ensemblemean(sol; samples=0)
    @test length(ensemblemean(sol; samples=1)) == 1
end

# Adapter deliberately has no `stepsize`. Fixed-window MoWi must depend on the
# solver protocol, not a Runge-Kutta field layout.
struct WindowAdapter{S} <: NSDEBase.AbstractInitialValueSolver
    solver::S
end
function NSDEBase.initialize_cache(p::NSDEBase.AbstractInitialValueProblem, s::WindowAdapter)
    return initialize_cache(p, s.solver)
end
function NSDEBase.initialize_solution(
    p::NSDEBase.AbstractInitialValueProblem, s::WindowAdapter; kwargs...
)
    return initialize_solution(p, s.solver; kwargs...)
end
function NSDEBase.solve!(c, sol, p::NSDEBase.AbstractInitialValueProblem, s::WindowAdapter)
    return solve!(c, sol, p, s.solver)
end
@testset "Generic fine adapters and window end caps" begin
    fine = WindowAdapter(RK4(h=0.01))
    par = Parareal(fine, Euler(h=0.05); N=2, criterion=NoStop())
    for adaptive in (nothing, StretchParameters(fineupdate=false, δτ⁺=10.0))
        sol = solve(Logistic(0.3, (0.0, 1.1)), MoWi(par; τ=0.5, Δτ=0.2, adaptive))
        @test NSDEMovingWindow.windowspan(sol[end])[2] == 1.1
        @test all(w -> NSDEMovingWindow.windowspan(w)[2] <= 1.1, sol.windows)
        @test sol.status == :converged
    end
    @test_throws ArgumentError MoWi(par; τ=0.5, adaptive=StretchParameters(fineupdate=true))
    @test_throws ArgumentError MoWi(
        remake(par; coarsesolver=fine); τ=0.5, adaptive=StretchParameters()
    )
end

@testset "Decimal windows snap to the final endpoint" begin
    par = Parareal(RK4(h=0.005), Euler(h=0.01); N=4, criterion=NoStop())
    for (t0, T, τ, Δτ, M) in (
        (0.0, 1.0, 0.1, 0.1, 10),
        (0.0, 3.0, 0.3, 0.3, 10),
        (0.0, 100.0, 0.1, 0.1, 1000),
        (0.0, 10.0, 1.0, 0.1, 91),
        (6000.0, 6001.0, 0.1, 0.1, 10),
    )
        sol = solve(Logistic(0.3, (t0, T)), MoWi(par; τ, Δτ))
        @test sol.status == :converged
        @test length(sol) == M
        @test NSDEMovingWindow.windowspan(sol[end])[2] === T
        @test !any(sol.clipped)
        @test all(w -> diff(collect(NSDEMovingWindow.windowspan(w)))[1] > τ / 2, sol.windows)
        @test all(isfinite, ensemblemean(sol; samples=3))
    end
end

@testset "Clipped window selection and vector observables" begin
    p = IVP(-1im, [1.0], (0.0, 2.0))
    mowi = MoWi(Parareal(RK4(h=0.01), Euler(h=0.05); N=2, criterion=NoStop()); τ=0.6, Δτ=0.3)
    sol = solve(p, mowi)
    @test sol.clipped == [false, false, false, false, false, true]
    X = NSDEMovingWindow.ensemblesamples(sol, first, 5)
    Xa = NSDEMovingWindow.ensemblesamples(sol, first, 5; include_clipped=true)
    @test size(X) == (5, 5)
    @test size(Xa) == (5, 6)
    @test X == Xa[:, 1:5]
    @test ensemblemean(sol; samples=5, include_clipped=true) ≈ vec(sum(Xa; dims=2)) / 6
    g(u) = [u[1], 2u[1]]
    for keep in (false, true)
        means = ensemblemean(sol, g; samples=5, include_clipped=keep)
        vars = ensemblevariance(sol, g; samples=5, include_clipped=keep)
        sems = ensemblesem(sol, g; samples=5, include_clipped=keep)
        mean1 = ensemblemean(sol; samples=5, include_clipped=keep)
        var1 = ensemblevariance(sol; samples=5, include_clipped=keep)
        sem1 = ensemblesem(sol; samples=5, include_clipped=keep)
        @test all(i -> means[i] ≈ [mean1[i], 2mean1[i]], 1:5)
        @test all(i -> vars[i] ≈ [var1[i], 4var1[i]], 1:5)
        @test all(i -> sems[i] ≈ [sem1[i], 2sem1[i]], 1:5)
    end
    short = solve(remake(p; tspan=(0.0, 0.2)), mowi)
    @test short.clipped == [true]
    @test_throws ArgumentError ensemblemean(short)
    @test length(ensemblemean(short; include_clipped=true)) == 100
    cache = NSDEMovingWindow.MoWiCache(p, mowi)
    result = MoWiSolution(p, mowi)
    mowi(cache, result, p)
    @test any(result.clipped)
    mowi(cache, result, remake(p; tspan=(0.0, 1.8)))
    @test !any(result.clipped)
end

@testset "Long decimal run with 256 chunks" begin
    N = 256
    p = IVP((u, t) -> ones(length(u)), [0.0], (0.0, 100.0))
    par = Parareal(
        RK4(h=0.1/N), Euler(h=0.1/N); N, criterion=WeightedDefect(ϵ=1e-10, weight=FixedWeight(1))
    )
    sol = solve(p, MoWi(par; τ=0.1, Δτ=0.1))
    @test sol.status == :converged
    @test length(sol) == 1000
    @test NSDEMovingWindow.windowspan(sol[end])[2] == 100.0
    @test !any(sol.clipped)
    @test all(
        m -> NSDEMovingWindow.windowspan(sol[m])[2] == NSDEMovingWindow.windowspan(sol[m + 1])[1],
        1:999,
    )
    @test sol[end](100.0)[1] ≈ 100.0 atol=1e-8
    @test all(isfinite, ensemblemean(sol; samples=3))
end

@testset "Compensated shifts and Leap rollback" begin
    start = 6000.0
    error = 0.0
    reference = BigFloat(start)
    for i in 1:20000
        shift = (0.1, 0.03, 0.2)[mod1(i, 3)]
        start, error = NSDEMovingWindow.window_shift(start, shift, error)
        reference += BigFloat(shift)
    end
    @test abs(BigFloat(start) - reference) <= eps(start)
    # A failed shift is replaced, not accumulated into the retry.
    previous = 0.6
    previous_error = 2e-17
    attempted, err = NSDEMovingWindow.window_shift(previous, 0.3, previous_error)
    retried, retry_error = NSDEMovingWindow.restart_start(
        LeapParameters(), attempted, previous, 0.1, err, previous_error
    )
    @test (retried, retry_error) == NSDEMovingWindow.window_shift(previous, 0.1, previous_error)
    @test NSDEMovingWindow.restart_start(
        StretchParameters(), attempted, previous, 0.1, err, previous_error
    ) == (attempted, err)
end
