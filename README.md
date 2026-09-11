# NSDEMovingWindow.jl

A Julia package implementing the moving-window (MoWi) algorithm.

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://giancarloantonucci.github.io/NSDEMovingWindow.jl/dev) ![Build Status](https://img.shields.io/github/actions/workflow/status/giancarloantonucci/NSDEMovingWindow.jl/CI.yml) ![Coverage Status](https://img.shields.io/codecov/c/github/giancarloantonucci/NSDEMovingWindow.jl)

## Installation

<!-- This package is a [registered package](https://juliahub.com/ui/Search?q=NSDEMovingWindow&type=packages) compatible with Julia v1.6 and above. From the Julia REPL,

```
]add NSDEMovingWindow
``` -->

This package is compatible with Julia v1.6 and above. From the Julia REPL,

```
]add https://github.com/giancarloantonucci/NSDEMovingWindow.jl
```

Read the [documentation](https://giancarloantonucci.github.io/NSDEMovingWindow.jl/dev) for a complete overview of this package.

## Usage

MoWi covers a long time span with a sequence of overlapping windows, each solved by a time-parallel subroutine seeded from the previous window over the overlap. The output is an **ensemble** of short window solutions, not one stitched trajectory — in a chaotic system every shift jumps orbit, by design — and the package ships the matching ensemble statistics.

```julia
using NSDERungeKutta, NSDETimeParallel, NSDEMovingWindow

problem = Lorenz([2.0, 3.0, -14.0], (0.0, 100.0))
parareal = Parareal(RK4(h = 1e-4), RK4(h = 1e-2);
                    parameters = PararealParameters(N = 8),
                    tolerance = Tolerance(ϵ = 1e-8))
mowi = MoWi(parareal; τ = 2.0, Δτ = 1.0, adaptive = StretchParameters())
solution = solve(problem, mowi)

solution.restarts                       # restarts per window
X̄ = ensemblemean(solution, u -> u[1])   # ensemble statistics over the windows
σ² = ensemblevariance(solution, u -> u[1])
```

## Strategies

The `adaptive` parameter type selects the strategy — there are no strategy integers:

- **Fixed** (`adaptive = nothing`, default) — fixed window length τ and shift Δτ, no restarts.
- **`StretchParameters`** (AMoWi-1) — adapt the window length: shrink τ (and the coarse step with it) on a failed window, stretch on success. The fine step stays fixed unless `fineupdate = true`.
- **`LeapParameters`** (AMoWi-2) — adapt the shift: hop back on failure, leap forward on success, capped at the window length. The first window is exempt from restarts — pick τ so it converges.
- **`ZoomParameters`** (AMoWi-3) — adapt the base weight of the weighted proximity function `ψ₂`, multiplicatively per restart: loosen on failure, tighten on success. Two constructor guards: the tolerance must use `ψ₂` (with `ψ₁`, zooming is a no-op), and the base weight must be **frozen** (`updatew = true` would silently undo every zoom-in) — measure `w` once with a probe window, then pass `Weights(w = ŵ)`. See the Strategies page of the docs for the recipe and tuning guidance.

Extra keyword arguments to `solve` (e.g. `mode = "THREADS"`, `saveiterates = true`) go straight to the time-parallel subroutine.
