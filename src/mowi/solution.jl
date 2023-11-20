"""
    MoWiSolution <: AbstractMovingWindowSolution

A composite type for an [`AbstractMovingWindowSolution`](@ref) obtained using an [`AbstractMovingWindowSolver`](@ref).

# Constructors
```julia
MoWiSolution(windows, restarts)
MoWiSolution(problem::AbstractInitialValueProblem, mowi::MoWi)
```

## Arguments
- `windows :: AbstractVector{𝕊} where 𝕊<:AbstractTimeParallelSolution`
- `restarts :: AbstractVector{ℤ} where ℤ<:Integer` : collects the number of restarts executed in all windows.

# Functions
- [`firstindex`](@ref) : first index.
- [`getindex`](@ref) : get window.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of windows.
- [`setindex!`](@ref) : set window.
"""
mutable struct MoWiSolution{windows_T<:(AbstractVector{𝕊} where 𝕊<:AbstractTimeParallelSolution), restarts_T<:(AbstractVector{ℤ} where ℤ<:Integer)} <: AbstractMovingWindowSolution
    windows::windows_T
    restarts::restarts_T
end

function MoWiSolution(problem::AbstractInitialValueProblem, mowi::MoWi)
    @↓ (t0, tN) ← tspan = problem
    @↓ τ, Δτ = mowi
    M = Δτ > 0 ? ceil(Int, (tN - τ) / Δτ) + 1 : 1 # Δτ == 0 => T = τ => M = 1
    windows = Vector{AbstractTimeParallelSolution}(undef, M)
    restarts = zeros(Integer, M)
    return MoWiSolution(windows, restarts)
end

#---------------------------------- FUNCTIONS ----------------------------------

"""
    length(solution::MoWiSolution)

returns the number of windows of `solution`.
"""
Base.length(solution::MoWiSolution) = length(solution.windows)

"""
    getindex(solution::MoWiSolution, m::Int)

returns the `m`-th window of `solution`.
"""
Base.getindex(solution::MoWiSolution, m::Integer) = solution.windows[m]

"""
    setindex!(solution::MoWiSolution, timeparallelsolution::AbstractTimeParallelSolution, m::Int)

stores a `timeparallelsolution` as the `m`-th window of `solution`.
"""
Base.setindex!(solution::MoWiSolution, timeparallelsolution::AbstractTimeParallelSolution, m::Integer) = solution.windows[m] = timeparallelsolution

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

returns a [`MoWiSolution`](@ref) constructor for the solution `problem` with `mowi`.
"""
MovingWindowSolution(problem::AbstractInitialValueProblem, mowi::MoWi) = MoWiSolution(problem, mowi)
