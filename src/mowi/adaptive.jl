# The adaptive strategies of AMoWi. The parameter type is the strategy: pass
# one of these as `adaptive` to `MoWi` and the matching policy runs. Each acts
# through two hooks, `adapt_failure!` (on a window restart) and
# `adapt_success!` (on acceptance), defined next to the window loop in
# `mowi.jl`. Combined strategies are not implemented.

"""
    StretchParameters <: AbstractAdaptiveMoWiParameters

AMoWi-1: adapt the window length τ (stretching and shrinking). On a failed
window the restart shrinks τ, Δτ and the coarse step by `δτ⁻`, with τ floored at
`N·h_F` and Δτ at one chunk as termination guards; on success they grow by
`δτ⁺`. The number of chunks N is fixed, so the coarse step scales with τ. The
fine step stays fixed by default, since changing it can alter the resolved
physics, unless `fineupdate = true`, which reproduces the idealised model used
in the convergence analysis rather than the practical scheme.

# Constructors
```julia
StretchParameters(; δτ⁻=0.5, δτ⁺=1.25, R=10, fineupdate=false)
```
`0 < δτ⁻ < 1`, `δτ⁺ ≥ 1`, and `R ≥ 0` are required. All factors must be finite.
"""
struct StretchParameters{δτ⁻_T<:Real,δτ⁺_T<:Real,R_T<:Integer,fineupdate_T<:Bool} <:
       AbstractAdaptiveMoWiParameters
    δτ⁻::δτ⁻_T
    δτ⁺::δτ⁺_T
    R::R_T
    fineupdate::fineupdate_T
    function StretchParameters(δτ⁻::Real, δτ⁺::Real, R::Integer, fineupdate::Bool)
        isfinite(δτ⁻) && 0 < δτ⁻ < 1 || throw(ArgumentError("δτ⁻ must lie in (0, 1)"))
        isfinite(δτ⁺) && δτ⁺ >= 1 || throw(ArgumentError("δτ⁺ must be finite and at least 1"))
        R >= 0 || throw(ArgumentError("R must be nonnegative"))
        return new{typeof(δτ⁻),typeof(δτ⁺),typeof(R),typeof(fineupdate)}(δτ⁻, δτ⁺, R, fineupdate)
    end
end

function StretchParameters(; δτ⁻::Real=0.5, δτ⁺::Real=1.25, R::Integer=10, fineupdate::Bool=false)
    return StretchParameters(δτ⁻, δτ⁺, R, fineupdate)
end

"""
    LeapParameters <: AbstractAdaptiveMoWiParameters

AMoWi-2: adapt the window shift Δτ (leaping and hopping) with the window length
τ fixed. On success the next shift grows, capped by the window length:
`Δτ ← min(δΔτ⁺ Δτ, τ)`, the overlap constraint that stops the algorithm
outrunning the speed at which the time-parallel subroutine can carry
information. On failure the shift shrinks by `δΔτ⁻`, floored at one chunk
length `τ/N` as a termination guard, and the window is re-attempted at the
reduced shift. The first window is exempt from restarts: the strategy leaves τ
alone, so if window 1 cannot converge, no shift policy can save it. Pick τ with
care.

# Constructors
```julia
LeapParameters(; δΔτ⁻=0.25, δΔτ⁺=2.0, R=10)
```
`0 < δΔτ⁻ < 1`, `δΔτ⁺ ≥ 1`, and `R ≥ 0` are required. All factors must be finite.
"""
struct LeapParameters{δΔτ⁻_T<:Real,δΔτ⁺_T<:Real,R_T<:Integer} <: AbstractAdaptiveMoWiParameters
    δΔτ⁻::δΔτ⁻_T
    δΔτ⁺::δΔτ⁺_T
    R::R_T
    function LeapParameters(δΔτ⁻::Real, δΔτ⁺::Real, R::Integer)
        isfinite(δΔτ⁻) && 0 < δΔτ⁻ < 1 || throw(ArgumentError("invalid failure factor δΔτ⁻"))
        isfinite(δΔτ⁺) && δΔτ⁺ >= 1 || throw(ArgumentError("invalid success factor δΔτ⁺"))
        R >= 0 || throw(ArgumentError("R must be nonnegative"))
        return new{typeof(δΔτ⁻),typeof(δΔτ⁺),typeof(R)}(δΔτ⁻, δΔτ⁺, R)
    end
end

LeapParameters(; δΔτ⁻::Real=0.25, δΔτ⁺::Real=2.0, R::Integer=10) = LeapParameters(δΔτ⁻, δΔτ⁺, R)

"""
    ZoomParameters <: AbstractAdaptiveMoWiParameters

AMoWi-3: adapt the base weight `w` of the weighted proximity function with τ
and Δτ fixed. The action is multiplicative and per restart: after ι restarts
the weight is `δ^ι w`. A failed window loosens the measure (`w ← δw⁻·w` with
`δw⁻ > 1`: a larger `w` makes the `w^{-Δt}` weights decay faster, so
discontinuities further along the window count for less, "zooming out"); a
success tightens it (`w ← δw⁺·w` with `δw⁺ ≤ 1`, clamped at `w ≥ 1` in the
stored configuration). Requires `WeightedDefect(weight=FixedWeight(w))`:
changing weights within a window would overwrite the zoom adjustment. Use a
probe run to measure a suitable base, then construct a fixed-weight criterion
for Zoom.

# Constructors
```julia
ZoomParameters(; δw⁻=2.0, δw⁺=0.25, R=10)
```
`δw⁻ > 1`, `0 < δw⁺ ≤ 1`, and `R ≥ 0` are required. All factors must be finite.
"""
struct ZoomParameters{δw⁻_T<:Real,δw⁺_T<:Real,R_T<:Integer} <: AbstractAdaptiveMoWiParameters
    δw⁻::δw⁻_T
    δw⁺::δw⁺_T
    R::R_T
    function ZoomParameters(δw⁻::Real, δw⁺::Real, R::Integer)
        isfinite(δw⁻) && δw⁻ > 1 || throw(ArgumentError("invalid failure factor δw⁻"))
        isfinite(δw⁺) && 0 < δw⁺ <= 1 || throw(ArgumentError("invalid success factor δw⁺"))
        R >= 0 || throw(ArgumentError("R must be nonnegative"))
        return new{typeof(δw⁻),typeof(δw⁺),typeof(R)}(δw⁻, δw⁺, R)
    end
end

ZoomParameters(; δw⁻::Real=2.0, δw⁺::Real=0.25, R::Integer=10) = ZoomParameters(δw⁻, δw⁺, R)
