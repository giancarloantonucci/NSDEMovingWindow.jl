# Executable companion to docs/src/guide.md.
# Run in an environment containing this package and the dependencies used below.

using NSDERungeKutta, NSDETimeParallel, NSDEMovingWindow
problem = Logistic(0.3, (0.0, 2.3))
par = Parareal(RK4(; h=0.01), RK4(; h=0.05); N=4, criterion=RelativeUpdate(; atol=1e-12, rtol=1e-9))
mowi = MoWi(par; τ=1.0, Δτ=0.5)
sol = solve(problem, mowi)
@assert sol.status == :converged
@assert length(sol) == 4
@assert sol.windows[end].lastiterate[end].t[end] == 2.3
@assert sol.clipped == [false, false, false, true]
nothing

stretched = solve(
    problem, MoWi(par; τ=1.0, Δτ=0.5, adaptive=StretchParameters(; δτ⁻=0.5, δτ⁺=1.1, R=3))
)
leaped = solve(problem, MoWi(par; τ=1.0, Δτ=0.5, adaptive=LeapParameters(; R=3)))
zpar = remake(par; criterion=WeightedDefect(; ϵ=1e-8, weight=FixedWeight(2.0)))
zoomed = solve(problem, MoWi(zpar; τ=1.0, Δτ=0.5, adaptive=ZoomParameters(; R=3)))
@assert all(s -> s.status == :converged, (stretched, leaped, zoomed))
nothing

μ = ensemblemean(sol, first; samples=20)
v = ensemblevariance(sol, first; samples=20)
se = ensemblesem(sol, first; samples=20)
vector_mean = ensemblemean(sol, u -> [u[1], u[1]^2]; samples=20)
μ_with_tail = ensemblemean(sol, first; samples=20, include_clipped=true)
@assert length(μ) == length(v) == length(se) == 20
@assert all(x -> x >= 0, v) && all(x -> x >= 0, se)
@assert all(x -> length(x) == 2, vector_mean)
nothing

cache = NSDEMovingWindow.MoWiCache(problem, mowi)
buffer = MoWiSolution(problem, mowi)
mowi(cache, buffer, remake(problem; u0=[0.4]))
@assert buffer.status == :converged

chaos = Lorenz([2.0, 3.0, -14.0], (0.0, 3.0))
cp = Parareal(RK4(; h=0.002), RK4(; h=0.02); N=4, criterion=RelativeUpdate(; rtol=1e-8))
csol = solve(chaos, MoWi(cp; τ=1.0, Δτ=0.5))
@assert csol.status == :converged
nothing
