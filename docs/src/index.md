# NSDEMovingWindow.jl

Overlapping time windows with Parareal as the local solver.

## First example

```julia
using NSDERungeKutta, NSDETimeParallel, NSDEMovingWindow
problem = Logistic(0.3, (0.0, 2.3))
par = Parareal(RK4(h=0.01), RK4(h=0.05); N=4)
solution = solve(problem, MoWi(par; τ=1.0, Δτ=0.5))
@assert solution.status == :converged
ensemblemean(solution, first; samples=20)
```

Start with the [guided tour](guide.md), then consult:

- [Strategies](strategies.md)
- [API](api.md)

## Environment

These pages describe the current development checkout (Julia 1.10 or later).
Use compatible versions of the NSDE packages in one project. For local
checkouts, develop each package with `Pkg.develop(path="/path/to/package")`;
`Pkg.status()` reports what will actually load. A registered release may expose
an older API. Plotting examples additionally need Plots.

## Endpoint windows and ensemble statistics

The problem's final time is a hard boundary: MoWi does not integrate beyond it.
Accepted window shifts use compensated summation. Leap retries restore the
previous accepted start and its compensation before applying the shorter shift;
Stretch and Zoom retain the attempted start. Window ends use the same addition,
so non-overlapping windows meet without rounding gaps.

An endpoint gap within `max(64eps(max(abs(t0), abs(tN))), 1e-9*τ_current)` is
snapped to the final time. This is a scheduling tolerance, not an ODE error
tolerance. A larger overshoot is clipped and marked in `solution.clipped`. The
ensemble mean, variance and SEM exclude these genuinely shortened windows by
default. Pass `include_clipped=true` to include them explicitly. If no full
windows remain, the default raises an error. A tolerance-sized endpoint snap
retains the window.

All statistics sample each retained window at the same relative positions. For
Stretch, full windows can still differ in duration: this compares relative
positions, not a common elapsed time. The observable may return a scalar or a
fixed-length numeric vector; variance and SEM are then componentwise. SEM assumes
independence and is not corrected for correlations between overlapping windows.

## Reference

G. A. Antonucci (2025), *Solving chaotic ODEs with time-parallel algorithms*,
DPhil thesis, University of Oxford. DOI:
[10.5287/ora-xoevp0dez](https://doi.org/10.5287/ora-xoevp0dez). The thesis
provides the mathematical background; the API and implementation conventions
are specified by these package docs.
