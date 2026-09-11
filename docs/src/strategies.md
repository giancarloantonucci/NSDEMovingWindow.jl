# Strategies

The window loop is written once; the strategies differ only in what they adjust when a window fails to converge (a restart) or is accepted. The `adaptive` parameter type selects the strategy — there are no strategy integers. Every strategy caps restarts per window at its `R`; a window that exhausts the budget is accepted with a warning.

## Fixed MoWi (`adaptive = nothing`)

Fixed window length τ and shift Δτ, no restarts. The baseline every adaptive variant reduces to.

## Stretching and shrinking: `StretchParameters` (AMoWi-1)

Adapt the window length. On a failed window the restart shrinks τ, Δτ and the coarse step by `δτ⁻`; on success they grow by `δτ⁺`. The number of chunks is fixed, so the coarse step scales with τ. The fine step stays **fixed** by default — changing it can alter the resolved physics — unless `fineupdate = true`, which reproduces the idealized model of the convergence analysis rather than the practical scheme.

## Leaping and hopping: `LeapParameters` (AMoWi-2)

Adapt the shift with the window length fixed. On success the next shift grows, capped by the window length — `Δτ ← min(δΔτ⁺ Δτ, τ)`, the overlap constraint that stops the algorithm out-running the speed at which the subroutine can carry information. On failure the shift shrinks and the window is re-attempted at the reduced shift. The first window is exempt from restarts: the strategy leaves τ alone, so if window 1 cannot converge, no shift policy can save it — pick τ with care.

## Zooming in and out: `ZoomParameters` (AMoWi-3)

Adapt the base weight of the **weighted** proximity function `ψ₂`, multiplicatively per restart: a failed window loosens the measure (`w ← δw⁻ · w`, discounting discontinuities further along the window), a success tightens it. Two things are required and **refused at construction** when absent:

1. the parallel solver's `Tolerance` must use `ψ₂` — with the unweighted `ψ₁` zooming is a silent no-op;
2. the base weight must be **frozen** — `Weights(updatew = true)` re-measures `w` every iteration and floors it at the measured rate, silently undoing every zoom-in.

The recipe: measure once with a probe, then freeze —

```julia
probe = Parareal(fine, coarse; parameters = params,
                 tolerance = Tolerance(ϵ = ϵ, ψ = ψ₂, weights = Weights(updatew = true)))
solve(copy(problem, u0, t0, t0 + τ), probe)          # one representative window
mowi = MoWi(Parareal(fine, coarse; parameters = params,
                     tolerance = Tolerance(ϵ = ϵ, ψ = ψ₂,
                                           weights = Weights(w = probe.tolerance.weights.w)));
            τ, Δτ, adaptive = ZoomParameters())
```

## Tuning

The adaptation factors are multiplicative per event, and the sweeps behind the thesis figures say the useful range is narrower than it looks. **Gentle factors recover; aggressive ones spiral.** A shrink factor around 2 (halve on failure, double on success) sits before the cliff; well past it, a failed window that shrinks hard fails *again* at the reduced size, shrinks again, and the run pays a long tail of tiny windows for one bad stretch. Loosening factors for Zoom behave the same way in `w`: `δw⁻` far above 2 forgives so much that the accepted window is barely constrained.

Three guard-rails are built in, so mis-tuning degrades rather than hangs. The shift is floored **absolutely** at `Δτ ≥ τ₀/N` (the original window over the chunk count — a floor relative to the *current* state would decay geometrically with it), window lengths are floored analogously under Stretch, and every window's restarts are capped at the strategy's `R`; a window that exhausts its budget is accepted unconverged, with a warning. That acceptance is deliberate: the alternative is non-termination, and an accepted-unconverged window is visible in `solution.restarts` and in the seams, where a hang is visible in nothing.

Two placement rules do more than any factor tuning. Pick τ so that **window 1 converges** — Leap exempts it from restarts because no shift policy can save a first window that cannot converge, and the other strategies only shrink their way out of it slowly. And keep the iteration budget honest: the per-window `K` your criterion demands must sit under `iteration_budget(S, N, ζ)` for your target speed-up, or the strategy is optimising a run that has already lost.
