# Strategies

The window loop is written once; the strategies differ only in what they adjust
when a window fails to converge (a restart) or is accepted. The `adaptive`
parameter type selects the strategy — there are no strategy integers. Every
strategy caps restarts per window at its `R`; a finite window that exhausts the
budget advances with a warning, while an unrecovered nonfinite window stops the
run.

## Fixed MoWi (`adaptive = nothing`)

Fixed window length τ and shift Δτ, no restarts. The baseline every adaptive
variant reduces to.

## Stretching and shrinking: `StretchParameters` (AMoWi-1)

Adapt the window length. On a failed window the restart shrinks τ, Δτ and the
coarse step by `δτ⁻`; on success they grow by `δτ⁺`. The number of chunks is
fixed, so the coarse step scales with τ. The fine step stays **fixed** by
default — changing it can alter the resolved physics — unless
`fineupdate = true`, which reproduces the idealised model of the convergence
analysis rather than the practical scheme. The fine step then scales with τ, so
it stops shrinking when τ reaches its floor.

## Leaping and hopping: `LeapParameters` (AMoWi-2)

Adapt the shift with the window length fixed. On success the next shift grows,
capped by the window length — `Δτ ← min(δΔτ⁺ Δτ, τ)`, the overlap constraint
that stops the algorithm out-running the speed at which the subroutine can carry
information. On failure the shift shrinks and the window is re-attempted at the
reduced shift. The first window is exempt from restarts: the strategy leaves τ
alone, so if window 1 cannot converge, no shift policy can save it — pick τ with
care.

## Zooming in and out: `ZoomParameters` (AMoWi-3)

Adapt the base of `WeightedDefect` multiplicatively between attempts. Zoom
requires `FixedWeight`: a same-chunk or neighbouring-chunk updater would choose
the weight independently of the window strategy. The adapted settings are local
to one `solve`; reusing the `MoWi` object starts from its configured settings.

A probe can estimate the initial base explicitly:

```julia
probe = Parareal(fine, coarse; N=32,
    criterion=WeightedDefect(ϵ=1e-9, weight=SecantWeight()))
measured = solve(remake(problem; tspan=(0.0, τ)), probe)
w0 = last(measured.diagnostics.weights)
fixed = remake(probe; criterion=WeightedDefect(ϵ=1e-9, weight=FixedWeight(w0)))
mowi = MoWi(fixed; τ, Δτ, adaptive=ZoomParameters())
```

The state carried between windows is an explicit replacement configuration;
neither the probe nor the caller's solver is mutated. Stretch similarly replaces
step-size settings, and chunk schedules are rebuilt when the window changes.

## Tuning

The adaptation factors act on each accepted or retried window. Strong shrinkage
can make subsequent windows unnecessarily short; strong growth can undo a useful
reduction. Large Zoom weights discount the tail more heavily, so early acceptance
alone does not establish tail accuracy. Choose factors with the required accuracy
and work costs in mind; there is no universal best range.

The shift is floored at `Δτ ≥ τ₀/N`, using the original window length, and
Stretch similarly floors the window length. Each window has a finite restart
budget. At its limit, a finite unconverged window advances with a warning and
leaves `solution.status == :budget_exhausted`; a nonfinite window stops the run
with `:nonfinite`. Restart counts and window diagnostics retain the evidence.

Two placement rules do more than any factor tuning. Pick τ so that **window 1
converges** — Leap exempts it from restarts because no shift policy can save a
first window that cannot converge, and the other strategies only shrink their
way out of it slowly. And keep the iteration budget honest: the per-window
fine-batch count `K` your criterion demands must sit under
`iteration_budget(S, N, ζ)` for your target speed-up, using the same cost model
for both. This is a work budget, not a convergence guarantee.

## Cache reuse and MPI windows

MoWi reuses caches when the new RHS has the same concrete type, refreshing the
window problem and allowing Parareal to rebuild its RHS-dependent caches. A
different RHS type requires a new MoWi cache.

The bulk backends use overlap values as initial guesses at interior chunk starts.
`PipelinedMPIBackend` uses only the seeded window initial value: its sequential
coarse start-up overwrites interior overlap guesses. Thus its next window has the
correct initial condition but does not benefit from those interior warm starts.
Use a bulk backend when overlap warm starts are needed. Simply skipping a coarse
receive would leave an unmatched MPI send; that protocol change is not implemented.

## Endpoints and ensemble uncertainty

The final time is a hard boundary. Compensated shifts prevent accumulated start
drift. Endpoints within the [scheduling tolerance](index.md) are snapped to it;
genuine overshoots are clipped and recorded in `solution.clipped`. This applies
also after Stretch growth. The last window can therefore be shorter than the
chosen window length. Fixed and Leap/Zoom windows do not require a fine solver
with a `stepsize` field. Stretch requires a resizable coarse solver;
`fineupdate=true` also requires a resizable fine solver.

Ensemble observables may be real or complex. `ensemblevariance` uses the
population mean squared modulus; `ensemblesem` instead uses sample variance
(divisor M−1) and requires at least two windows. Its standard-error formula
assumes independence: seeding and overlapping windows introduce correlations,
so it is not a correlation-adjusted uncertainty estimate. Positive correlations
make the formula optimistic; there is no general lower-bound guarantee.
