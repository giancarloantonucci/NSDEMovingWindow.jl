# NSDEMovingWindow.jl

This is the documentation of [NSDEMovingWindow.jl](https://github.com/giancarloantonucci/NSDEMovingWindow.jl), a Julia package implementing the moving-window (MoWi) algorithm for long-time integration of chaotic systems, with time parallelization as a subroutine.

## Installation

From the Julia REPL,

```
]add https://github.com/giancarloantonucci/NSDEMovingWindow.jl
```

## Getting started

The time span is covered by overlapping windows; each window is solved by Parareal (via NSDETimeParallel.jl), seeded from the previous window's solution over the overlap. The output is an **ensemble** of short window solutions — each shift jumps orbit in a chaotic system, by design — and the package ships the matching statistics.

```julia
using NSDERungeKutta, NSDETimeParallel, NSDEMovingWindow

problem = Lorenz([2.0, 3.0, -14.0], (0.0, 100.0))
parareal = Parareal(RK4(h = 1e-4), RK4(h = 1e-2);
                    parameters = PararealParameters(N = 8),
                    tolerance = Tolerance(ϵ = 1e-8))
mowi = MoWi(parareal; τ = 2.0, Δτ = 1.0, adaptive = StretchParameters())
solution = solve(problem, mowi)

solution.restarts                        # restarts per window
X̄ = ensemblemean(solution, u -> u[1])    # ensemble statistics over the windows
σ² = ensemblevariance(solution, u -> u[1])
sem = ensemblesem(solution, u -> u[1])   # precision of the mean: shrinks as 1/√M
```

Extra keyword arguments to `solve` (e.g. `mode = "THREADS"`, `saveiterates = true`) go straight to the time-parallel subroutine. `solution(t)` evaluates the latest window covering `t` — a plotting and diagnostic helper, not a single continuous trajectory.

- The [Strategies](strategies.md) page covers fixed MoWi and the three adaptive variants.
- The [API](api.md) holds the full reference.
