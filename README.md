# NSDEMovingWindow.jl

Overlapping time windows with Parareal as the local solver.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://giancarloantonucci.github.io/NSDEMovingWindow.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://giancarloantonucci.github.io/NSDEMovingWindow.jl/dev) ![Build Status](https://img.shields.io/github/actions/workflow/status/giancarloantonucci/NSDEMovingWindow.jl/CI.yml) ![Coverage Status](https://img.shields.io/codecov/c/github/giancarloantonucci/NSDEMovingWindow.jl)

## Start here

```julia
using NSDERungeKutta, NSDETimeParallel, NSDEMovingWindow
problem = Logistic(0.3, (0.0, 2.3))
par = Parareal(RK4(h=0.01), RK4(h=0.05); N=4)
solution = solve(problem, MoWi(par; τ=1.0, Δτ=0.5))
@assert solution.status == :converged
ensemblemean(solution, first; samples=20)
```

Each accepted window seeds the next one over the overlap. The result retains
all window solutions and diagnostics; `solution(t)` selects the latest covering
window and is not guaranteed to be continuous. Window shifts do not necessarily
change orbit, and windows are not independent statistical samples.

The final time is a hard boundary. Genuinely shortened last windows are marked
in `solution.clipped` and excluded from ensemble statistics by default.
`include_clipped=true` opts in. Statistics compare equal **relative positions**
inside retained windows; they do not establish invariant-measure accuracy.

Fixed windows are the default. [Strategies](docs/src/strategies.md) explains
Stretch (length), Leap (shift), and Zoom (fixed discount base). A finite window
that exhausts its budget advances with a warning; an unrecovered nonfinite
window stops the run and makes ensemble statistics unavailable.

## Learn by running

The [guided tour](docs/src/guide.md) explains each step. Its companion
[`examples/guided_tour.jl`](examples/guided_tour.jl) runs the numerical examples
with assertions. Run it in a project containing this package and its example
dependencies; use Julia 1.10 or later.

The ecosystem is NSDEBase → NSDERungeKutta → NSDETimeParallel → NSDEMovingWindow.
ArrowMacros provides supporting macros; they are not required in user examples.

## Installation and versions

For a registered release, use `Pkg.add("NSDEMovingWindow")`. These docs describe the
current development API; a registered version may differ. To use a local
checkout, activate your working environment and `Pkg.develop(path="/path/to/NSDEMovingWindow")`.
Develop matching local dependency checkouts in that environment too, then check
`Pkg.status()`. The serial and parallel guides also require NSDERungeKutta;
plotting requires Plots.

See the [documentation](https://giancarloantonucci.github.io/NSDEMovingWindow.jl/dev/)
and [API reference](docs/src/api.md) for details.

## Build the manual locally

From the package root, with matching dependency checkouts already available:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs --threads=4 docs/make.jl
```

Documenter executes the guided-tour examples during the build. For local
unreleased dependencies, develop their paths into `docs` before instantiating.
The deployment step only publishes in an appropriately configured CI environment.
