# src/mowi/adaptive.jl
#
# The adaptive strategies of AMoWi (thesis, Ch. 5). The parameter type IS the
# strategy: pass one of these as `adaptive` to `MoWi` and the matching policy
# runs — no strategy integers. Each acts through two hooks, `adapt_failure!`
# (on a window restart) and `adapt_success!` (on acceptance), defined next to
# the window skeleton in `mowi.jl`. The unexplored combined strategy (thesis
# §5.4, future work) has no type here on purpose.

"""
    StretchParameters <: AbstractAdaptiveMoWiParameters

AMoWi-1: adapt the window LENGTH τ (stretching/shrinking). On a failed window
the restart shrinks τ, Δτ and the coarse step by `δτ⁻` (τ floored at `N·h_F`, Δτ at one chunk — termination guards); on success they grow
by `δτ⁺`. The number of chunks N is fixed, so the coarse step scales with τ.
The fine step stays FIXED by default — changing it can alter the resolved
physics — unless `fineupdate = true`, which reproduces the idealized model
used in the convergence analysis (thesis Thm. 5.1) rather than the practical
scheme.

# Constructors
```julia
StretchParameters(; δτ⁻=0.5, δτ⁺=1.25, R=10, fineupdate=false)
```
Thesis experiments put the sweet spot at `δτ⁺ ∈ [1.25, 1.75]`, with a
performance cliff past 2; sensitivity to `δτ⁻` is mild. `R` caps restarts per
window.
"""
struct StretchParameters{δτ⁻_T<:Real, δτ⁺_T<:Real, R_T<:Integer, fineupdate_T<:Bool} <: AbstractAdaptiveMoWiParameters
    δτ⁻::δτ⁻_T
    δτ⁺::δτ⁺_T
    R::R_T
    fineupdate::fineupdate_T
end

StretchParameters(; δτ⁻::Real=0.5, δτ⁺::Real=1.25, R::Integer=10, fineupdate::Bool=false) = StretchParameters(δτ⁻, δτ⁺, R, fineupdate)

"""
    LeapParameters <: AbstractAdaptiveMoWiParameters

AMoWi-2: adapt the window SHIFT Δτ (leaping/hopping) with the window length τ
fixed. On success the next shift grows, capped by the window length —
`Δτ ← min(δΔτ⁺ Δτ, τ)` — the overlap constraint that stops the algorithm
out-running the speed at which the time-parallel subroutine can carry
information. On failure the shift shrinks by `δΔτ⁻` — floored at one chunk
length `τ/N`, a termination guard — and the window is re-attempted at the
reduced shift. The FIRST window is exempt from restarts:
the strategy leaves τ alone, so if window 1 cannot converge, no shift policy
can save it (thesis §5.2) — pick τ with care.

# Constructors
```julia
LeapParameters(; δΔτ⁻=0.25, δΔτ⁺=2.0, R=10)
```
"""
struct LeapParameters{δΔτ⁻_T<:Real, δΔτ⁺_T<:Real, R_T<:Integer} <: AbstractAdaptiveMoWiParameters
    δΔτ⁻::δΔτ⁻_T
    δΔτ⁺::δΔτ⁺_T
    R::R_T
end

LeapParameters(; δΔτ⁻::Real=0.25, δΔτ⁺::Real=2.0, R::Integer=10) = LeapParameters(δΔτ⁻, δΔτ⁺, R)

"""
    ZoomParameters <: AbstractAdaptiveMoWiParameters

AMoWi-3: adapt the base weight `w` of the WEIGHTED proximity function with τ
and Δτ fixed. The action is multiplicative and per restart (thesis Thm. 5.3:
after ι restarts the weight is `δ^ι w`): a failed window loosens the measure
(`w ← δw⁻·w`, `δw⁻ > 1` — larger `w` makes the `w^{-Δt}` weights decay
faster, so discontinuities further along the window count for less — "zoom
out"); a success tightens it (`w ← δw⁺·w`, `δw⁺ < 1`, clamped at `w ≥ 1`
inside `ψ₂` — "zoom in"). Requires the parallel solver's tolerance to use a
WEIGHTED error function (`ψ₂`); with the unweighted `ψ₁` zooming would be a
silent no-op, so `MoWi` refuses that combination at construction. It likewise
requires a FROZEN base weight: `Weights(updatew = true)` re-measures `w`
every iteration and floors it at the measured rate, silently undoing every
zoom-in, so `MoWi` refuses that pairing too — measure Λ once with a probe
solve and pass `Weights(w = exp(Λ))`.

# Constructors
```julia
ZoomParameters(; δw⁻=2.0, δw⁺=0.25, R=10)
```
"""
struct ZoomParameters{δw⁻_T<:Real, δw⁺_T<:Real, R_T<:Integer} <: AbstractAdaptiveMoWiParameters
    δw⁻::δw⁻_T
    δw⁺::δw⁺_T
    R::R_T
end

# @! Adaptive per-window iteration BUDGETS (the AMoWi-3 ↔ budget coupling of
# @! rem:amowi3_limitation) are not expressible today: PararealParameters.K is
# @! fixed for the whole run. The MoWA-era precedent (`Budget(Kmax)` on the
# @! window driver) says where the feature belongs if built — MoWi-LEVEL
# @! state adapted by the strategy hooks, never a solver parameter — so the
# @! budget can move with the window the way Δτ and w already do.

ZoomParameters(; δw⁻::Real=2.0, δw⁺::Real=0.25, R::Integer=10) = ZoomParameters(δw⁻, δw⁺, R)
