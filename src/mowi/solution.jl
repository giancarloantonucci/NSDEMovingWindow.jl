# src/mowi/solution.jl

"""
    MoWiSolution <: AbstractMovingWindowSolution

The output of [`MoWi`](@ref): a vector of window solutions plus the restart
count per window. `clipped` flags windows shortened by the global endpoint
(tolerance-sized snaps are not marked). `status` is `:not_started`, `:converged` (all windows accepted
or completed exactly), `:budget_exhausted` (at least one finite window used up its
budget), or `:nonfinite`. A nonfinite run stops after its last retry and retains
the failed window for inspection; ensemble statistics then raise `ArgumentError`.
Read it as an ENSEMBLE of short-time solutions, not as one
guaranteed continuous trajectory. Overlapping approximations can disagree;
a shift need not change orbit when the solves agree. The statistics
([`ensemblemean`](@ref), [`ensemblevariance`](@ref), [`ensemblesem`](@ref))
compare equal relative positions across windows, which need not be independent.

# Constructors
```julia
MoWiSolution(windows, restarts)
MoWiSolution(problem::AbstractInitialValueProblem, mowi::MoWi)
```
"""
mutable struct MoWiSolution{
    windows_T<:(AbstractVector{𝕊} where 𝕊<:AbstractTimeParallelSolution),
    restarts_T<:AbstractVector{Int},
} <: AbstractMovingWindowSolution
    windows::windows_T
    restarts::restarts_T
    status::Symbol
    clipped::Vector{Bool}
end

MoWiSolution(windows, restarts) = MoWiSolution(windows, restarts, :not_started)
function MoWiSolution(windows, restarts, status::Symbol)
    return MoWiSolution(windows, restarts, status, fill(false, length(windows)))
end

function MoWiSolution(problem::AbstractInitialValueProblem, mowi::MoWi)
    @↓ (t0, tN) ← tspan = problem
    @↓ τ, Δτ = mowi
    M = max(ceil(Int, (tN - t0 - τ) / Δτ) + 1, 1)
    windows = Vector{AbstractTimeParallelSolution}(undef, M)
    restarts = zeros(Int, M)
    return MoWiSolution(windows, restarts)
end

#----------------------------------- METHODS -----------------------------------

function windowspan(window::AbstractTimeParallelSolution)
    return (window.lastiterate[1].t[1], window.lastiterate[end].t[end])
end

"""
    (solution::MoWiSolution)(tₚ::Real)

evaluates the LATEST window covering `tₚ` (clamped at the ends). A diagnostic
and plotting helper: windows overlap and may disagree, so this is NOT guaranteed to be a single continuous trajectory — for
statistics use [`ensemblemean`](@ref) and friends instead.
"""
function (solution::MoWiSolution)(tₚ::Real)
    for m in length(solution):-1:1
        ta, tb = windowspan(solution[m])
        if ta ≤ tₚ ≤ tb
            return solution[m](tₚ)
        end
    end
    ta1, _ = windowspan(solution[1])
    return tₚ < ta1 ? solution[1](tₚ) : solution[end](tₚ)
end

#---------------------------------- FUNCTIONS ----------------------------------

"""
    length(solution::MoWiSolution)

returns the number of windows of `solution`.
"""
Base.length(solution::MoWiSolution) = length(solution.windows)

"""
    getindex(solution::MoWiSolution, m::Integer)

returns the `m`-th window of `solution`.
"""
Base.getindex(solution::MoWiSolution, m::Integer) = solution.windows[m]

"""
    setindex!(solution::MoWiSolution, window::AbstractTimeParallelSolution, m::Integer)

stores `window` as the `m`-th window of `solution`.
"""
function Base.setindex!(solution::MoWiSolution, window::AbstractTimeParallelSolution, m::Integer)
    return solution.windows[m] = window
end

"""
    firstindex(solution::MoWiSolution)

returns the first index of `solution`.
"""
Base.firstindex(solution::MoWiSolution) = firstindex(solution.windows)

"""
    lastindex(solution::MoWiSolution)

returns the last index of `solution`.
"""
Base.lastindex(solution::MoWiSolution) = lastindex(solution.windows)

"""
    MovingWindowSolution(problem::AbstractInitialValueProblem, mowi::MoWi)

returns a [`MoWiSolution`](@ref) for solving `problem` with `mowi`.
"""
MovingWindowSolution(problem::AbstractInitialValueProblem, mowi::MoWi) = MoWiSolution(problem, mowi)

#--------------------------------- STATISTICS ----------------------------------
# The windows form an ensemble of time series. Windows may have
# unequal lengths under AMoWi-1, so each is sampled on ITS OWN relative time
# grid s ∈ [0, 1]; the statistics are taken across windows at equal s.

"""
    ensemblemean(solution::MoWiSolution, g::Function=first; samples::Integer=100, include_clipped::Bool=false)

the ensemble mean of the functional `g(u)` across windows, sampled at
`samples` equally spaced relative times over each window's own span. Returns
a vector of length `samples`. `g` may return a number or a fixed-length numeric
vector; vector-valued results are componentwise. Windows shortened at the global
endpoint are excluded unless `include_clipped=true`. A tolerance-sized endpoint
snap does not exclude a full window. Stretch windows still use their own relative
time grids; these are comparisons at equal relative position, not equal elapsed
time. `solution.clipped` records the endpoint truncations.
"""
function ensemblemean end

function ensemblesamples(solution::MoWiSolution, g, samples::Integer; include_clipped::Bool=false)
    solution.status == :nonfinite && throw(
        ArgumentError(
            "ensemble statistics are unavailable after a nonfinite window; inspect solution.windows",
        ),
    )
    samples >= 1 || throw(ArgumentError("samples must be positive"))
    isempty(solution.windows) && throw(ArgumentError("ensemble contains no windows"))
    windows = include_clipped ? solution.windows : solution.windows[.!solution.clipped]
    isempty(windows) && throw(
        ArgumentError(
            "ensemble contains no full windows; pass include_clipped=true to include endpoint-truncated windows",
        ),
    )
    grid = samples == 1 ? (0.0,) : range(0.0, 1.0; length=samples)
    return [
        begin
            ta, tb = windowspan(w)
            g(w(ta + s * (tb - ta)))
        end for s in grid, w in windows
    ]
end

function ensemblemean(
    solution::MoWiSolution, g::Function=first; samples::Integer=100, include_clipped::Bool=false
)
    X = ensemblesamples(solution, g, samples; include_clipped)
    return vec(sum(X; dims=2)) ./ size(X, 2)
end

"""
    ensemblevariance(solution::MoWiSolution, g::Function=first; samples::Integer=100, include_clipped::Bool=false)

Population variance across windows at equal relative times. For complex-valued
observables this is the real mean squared modulus of the centred samples,
componentwise for vector observables. Endpoint-truncated windows are excluded
unless `include_clipped=true`, as in `ensemblemean`.
This describes spread, not uncertainty of the mean.
"""
function ensemblevariance(
    solution::MoWiSolution, g::Function=first; samples::Integer=100, include_clipped::Bool=false
)
    X = ensemblesamples(solution, g, samples; include_clipped)
    centre = sum(X; dims=2) ./ size(X, 2)
    return vec(sum(x -> abs2.(x), X .- centre; dims=2)) ./ size(X, 2)
end

"""
    ensemblesem(solution::MoWiSolution, g::Function=first; samples::Integer=100, include_clipped::Bool=false)

Independence-based estimate of the mean's standard error, using sample variance
(divisor M-1) and at least two windows. Windows are seeded from earlier windows
and may overlap, so this is not a correlation-adjusted uncertainty estimate.
Positive correlations make it optimistic; it is not a guaranteed lower bound.
Vector observables are handled componentwise. Endpoint-truncated windows are
excluded unless `include_clipped=true`; M counts the retained windows.
"""
function ensemblesem(
    solution::MoWiSolution, g::Function=first; samples::Integer=100, include_clipped::Bool=false
)
    X = ensemblesamples(solution, g, samples; include_clipped)
    M = size(X, 2)
    M >= 2 || throw(ArgumentError("ensemblesem requires at least two windows"))
    centre = sum(X; dims=2) ./ M
    return map(x -> sqrt.(x ./ (M * (M - 1))), vec(sum(x -> abs2.(x), X .- centre; dims=2)))
end

"""
    flatten(solution::MoWiSolution)

concatenates every window's [`flatten`](@ref)ed trajectory into one `(u, t)`
pair. Windows OVERLAP and their approximations may disagree: `t`
is non-monotone at window seams, and overlap regions appear once per covering
window. Repeated overlap samples change sampling weights; concatenation alone does not
justify attractor statistics or a single continuous trajectory. For relative-time
ensemble summaries, use [`ensemblemean`](@ref) and friends.
"""
function NSDETimeParallel.flatten(solution::MoWiSolution)
    u, t = flatten(solution[1])
    for m in 2:length(solution)
        uₘ, tₘ = flatten(solution[m])
        append!(u, uₘ)
        append!(t, tₘ)
    end
    return u, t
end
