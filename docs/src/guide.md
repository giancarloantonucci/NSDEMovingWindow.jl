# Guided tour: moving windows

MoWi solves successive, possibly overlapping windows with Parareal. The preceding
accepted window supplies the next initial state and interior overlap guesses.
Shorter windows can make each local solve easier, at the cost of more window
transitions. This is an algorithmic strategy, not a statistical-accuracy theorem.

Run this page in order or use `examples/guided_tour.jl`.

## Fixed windows first

```@example windows
using NSDERungeKutta, NSDETimeParallel, NSDEMovingWindow
problem = Logistic(0.3, (0.0, 2.3))
par = Parareal(RK4(h=0.01), RK4(h=0.05); N=4,
    criterion=RelativeUpdate(atol=1e-12, rtol=1e-9))
mowi = MoWi(par; τ=1.0, Δτ=0.5)
sol = solve(problem, mowi)
@assert sol.status == :converged
@assert length(sol) == 4
@assert sol.windows[end].lastiterate[end].t[end] == 2.3
@assert sol.clipped == [false, false, false, true]
nothing # hide
```

The requested final time is a hard boundary. The last window here is `[1.5,2.3]`,
shorter than the nominal length one. Genuine truncations are marked `clipped`;
roundoff-sized snaps are not. Compensated shifts prevent cumulative addition
drift. The exact scheduling tolerance is given on the [home page](index.md).
It is separate from a solver tolerance.

`sol.windows[m]` is a Parareal solution; inspect its status, work counts, residuals,
and seams. `sol(t)` selects the latest covering window. Adjacent selections can
have jumps, so this convenience view is not guaranteed to be a continuous fine
trajectory. A shift does not necessarily change orbit: with matching exact solves,
overlap values can agree. Finite errors and early acceptance can make them differ.

## Adapt one quantity at a time

```@example windows
stretched = solve(problem, MoWi(par; τ=1.0, Δτ=0.5,
    adaptive=StretchParameters(δτ⁻=0.5, δτ⁺=1.1, R=3)))
leaped = solve(problem, MoWi(par; τ=1.0, Δτ=0.5,
    adaptive=LeapParameters(R=3)))
zpar = remake(par; criterion=WeightedDefect(ϵ=1e-8, weight=FixedWeight(2.0)))
zoomed = solve(problem, MoWi(zpar; τ=1.0, Δτ=0.5,
    adaptive=ZoomParameters(R=3)))
@assert all(s -> s.status == :converged, (stretched, leaped, zoomed))
nothing # hide
```

| Strategy | Changes | What to check |
|---|---|---|
| Fixed | Neither length nor shift | Choose a window that works from the start |
| Stretch | Window length, shift, and coarse step | Fine step stays fixed unless `fineupdate=true` |
| Leap | Shift, capped by window length | First-window failure cannot be repaired by shifting |
| Zoom | Fixed discount base between attempts | Larger weights relax late-defect acceptance |

Zoom requires `WeightedDefect` with `FixedWeight`; it cannot compete with an
in-window secant update. See [Strategies](strategies.md) for factors, floors,
and an explicit probe-then-freeze example. These small runs demonstrate the API;
they are not evidence that one strategy is best on chaotic systems.

`R` limits retries per window. A finite window that exhausts its budget advances
with a warning and leaves the aggregate status `:budget_exhausted`. A nonfinite
window that cannot be recovered stops the run with `:nonfinite`; later windows
are not computed. MoWi's aggregate success is `:converged`, even when individual
windows finish as `:exact`. It does not have a separate successful `:exact` status.

## What the ensemble statistics mean

```@example windows
μ = ensemblemean(sol, first; samples=20)
v = ensemblevariance(sol, first; samples=20)
se = ensemblesem(sol, first; samples=20)
vector_mean = ensemblemean(sol, u -> [u[1], u[1]^2]; samples=20)
μ_with_tail = ensemblemean(sol, first; samples=20, include_clipped=true)
@assert length(μ) == length(v) == length(se) == 20
@assert all(x -> x >= 0, v) && all(x -> x >= 0, se)
@assert all(x -> length(x) == 2, vector_mean)
nothing # hide
```

At each relative position `s∈[0,1]`, these functions sample
`window(t_start+s*(t_end-t_start))` and average **across retained windows**.
The result is a vector over relative positions, not a single time average or
an estimate automatically associated with an invariant measure. Clipped windows
are excluded by default; `include_clipped=true` opts in. Stretch can produce
unequal full-window lengths too, so equal relative positions need not mean equal
elapsed times.

The observable may be scalar or a fixed-length numeric vector. Complex means are
complex; variance uses squared modulus and is real, componentwise for vectors.
Variance uses divisor `M`; SEM uses sample variance (`M-1`) divided by `M` before
taking the square root and requires at least two retained windows.

SEM assumes independent samples. Overlapping, sequentially seeded windows are
correlated; positive correlations can make it optimistic, but it is not a general
lower bound on uncertainty. Assess burn-in, dependence, observable choice, and
sampling separately. After `:nonfinite`, statistics raise an error. `flatten(sol)`
includes overlapping times repeatedly and does not provide uniform sampling of a
single global trajectory.

## Reuse and a chaotic problem

```@example windows
cache = NSDEMovingWindow.MoWiCache(problem, mowi)
buffer = MoWiSolution(problem, mowi)
mowi(cache, buffer, remake(problem; u0=[0.4]))
@assert buffer.status == :converged

chaos = Lorenz([2.0, 3.0, -14.0], (0.0, 3.0))
cp = Parareal(RK4(h=0.002), RK4(h=0.02); N=4,
    criterion=RelativeUpdate(rtol=1e-8))
csol = solve(chaos, MoWi(cp; τ=1.0, Δτ=0.5))
@assert csol.status == :converged
nothing # hide
```

`MoWiCache` is a qualified low-level constructor; this version does not implement
`initialize_cache(problem,mowi)` or `initialize_solution(problem,mowi)`.
For simpler solution-buffer reuse, `mowi(buffer, problem)` creates the cache for you.

Start with local accuracy, statuses, work counts, and unweighted seams. Compare
against a finer calculation on a tractable interval before extending the horizon.
The configuration is immutable; each run resets adaptation and restart counts.
For cache compatibility and pipeline warm-start limitations, see [Strategies](strategies.md).

To build intuition, vary window length with fine/coarse steps fixed, then vary the
shift. Record total fine/coarse calls and wall time, including restarts. Shorter
windows do not give free convergence: the horizon, discretisation, overlap, and
synchronization costs all matter.
