# src/mowi/solution.jl

"""
    MoWiSolution <: AbstractMovingWindowSolution

The output of [`MoWi`](@ref): a vector of window solutions plus the restart
count per window. Read it as an ENSEMBLE of short-time solutions, not as one
stitched trajectory — in a chaotic system every window shift jumps orbit, and
that is the feature the statistics ([`ensemblemean`](@ref),
[`ensemblevariance`](@ref), [`ensemblesem`](@ref)) build on.

`windows` deliberately keeps an abstract element type: it is filled once per
window (cold path, tens of entries), so boxing here costs nothing — the hot
seams are the CHUNKS inside each window, which are concrete by construction.

# Constructors
```julia
MoWiSolution(windows, restarts)
MoWiSolution(problem::AbstractInitialValueProblem, mowi::MoWi)
```
"""
mutable struct MoWiSolution{windows_T<:(AbstractVector{𝕊} where 𝕊<:AbstractTimeParallelSolution), restarts_T<:AbstractVector{Int}} <: AbstractMovingWindowSolution
    windows::windows_T
    restarts::restarts_T
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

windowspan(window::AbstractTimeParallelSolution) = (window.lastiterate[1].t[1], window.lastiterate[end].t[end])

"""
    (solution::MoWiSolution)(tₚ::Real)

evaluates the LATEST window covering `tₚ` (clamped at the ends). A diagnostic
and plotting helper: windows overlap and, for chaotic systems, sit on
different orbits, so this is NOT a single continuous trajectory — for
statistics use [`ensemblemean`](@ref) and friends instead.
"""
function (solution::MoWiSolution)(tₚ::Real)
    for m = length(solution):-1:1
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
Base.setindex!(solution::MoWiSolution, window::AbstractTimeParallelSolution, m::Integer) = solution.windows[m] = window

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
# Thesis §4.2: the windows form an ensemble of time series. Windows may have
# unequal lengths under AMoWi-1, so each is sampled on ITS OWN relative time
# grid s ∈ [0, 1]; the statistics are taken across windows at equal s.

"""
    ensemblemean(solution::MoWiSolution, g::Function=first; samples::Integer=100)

the ensemble mean of the functional `g(u)` across windows, sampled at
`samples` equally spaced relative times over each window's own span. Returns
a vector of length `samples`.
"""
function ensemblemean(solution::MoWiSolution, g::Function=first; samples::Integer=100)
    M = length(solution)
    grid = range(0.0, 1.0; length=samples)
    X̄ = zeros(samples)
    for m = 1:M
        ta, tb = windowspan(solution[m])
        for (j, s) in enumerate(grid)
            X̄[j] += g(solution[m](ta + s * (tb - ta)))
        end
    end
    X̄ ./= M
    return X̄
end

"""
    ensemblevariance(solution::MoWiSolution, g::Function=first; samples::Integer=100)

the ensemble variance of `g(u)` across windows on the relative time grid
(population normalization, matching the thesis). Quantifies the physical
spread of the ensemble; it does NOT shrink with more windows — for the
precision of the mean, see [`ensemblesem`](@ref).
"""
function ensemblevariance(solution::MoWiSolution, g::Function=first; samples::Integer=100)
    M = length(solution)
    grid = range(0.0, 1.0; length=samples)
    X̄ = ensemblemean(solution, g; samples)
    σ² = zeros(samples)
    for m = 1:M
        ta, tb = windowspan(solution[m])
        for (j, s) in enumerate(grid)
            σ²[j] += abs2(g(solution[m](ta + s * (tb - ta))) - X̄[j])
        end
    end
    σ² ./= M
    return σ²
end

"""
    ensemblesem(solution::MoWiSolution, g::Function=first; samples::Integer=100)

the standard error of the ensemble mean, `σ/√M`: the precision of
[`ensemblemean`](@ref), which — unlike the variance — shrinks as the window
count grows.
"""
function ensemblesem(solution::MoWiSolution, g::Function=first; samples::Integer=100)
    M = length(solution)
    return sqrt.(ensemblevariance(solution, g; samples) ./ M)
end

"""
    flatten(solution::MoWiSolution)

concatenates every window's [`flatten`](@ref)ed trajectory into one `(u, t)`
pair. Windows OVERLAP and, on chaotic problems, sit on different orbits: `t`
is non-monotone at window seams, and overlap regions appear once per covering
window. Fit for attractor statistics (Poincaré sections, histograms); wrong
for anything assuming a single continuous trajectory — for statistics on the
ensemble as an ensemble, use [`ensemblemean`](@ref) and friends.
"""
function NSDETimeParallel.flatten(solution::MoWiSolution)
    u, t = flatten(solution[1])
    for m = 2:length(solution)
        uₘ, tₘ = flatten(solution[m])
        append!(u, uₘ)
        append!(t, tₘ)
    end
    return u, t
end
