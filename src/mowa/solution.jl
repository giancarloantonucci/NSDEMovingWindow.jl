"""
    MovingWindowSolution <: AbstractMovingWindowSolution

A composite type for an [`AbstractMovingWindowSolution`](@ref) obtained using an [`AbstractMovingWindowSolver`](@ref).

# Constructors
```julia
MovingWindowSolution(windows, restarts)
MovingWindowSolution(problem::AbstractInitialValueProblem, mowa::MovingWindowSolver)
```

## Arguments
- `windows :: AbstractVector{𝕊} where 𝕊<:AbstractTimeParallelSolution`
- `restarts :: AbstractVector{ℤ} where ℤ<:Integer` : counts window restarts in adaptive loop.

# Functions
- [`firstindex`](@ref) : first index.
- [`getindex`](@ref) : get window.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of windows.
- [`setindex!`](@ref) : set window.
"""
mutable struct MovingWindowSolution{windows_T<:(AbstractVector{𝕊} where 𝕊<:AbstractTimeParallelSolution), restarts_T<:(AbstractVector{ℤ} where ℤ<:Integer)} <: AbstractMovingWindowSolution
    windows::windows_T
    restarts::restarts_T
end

function MovingWindowSolution(problem::AbstractInitialValueProblem, mowa::MovingWindowSolver)
    @↓ (t0, tN) ← tspan = problem
    @↓ τ, Δτ = mowa
    # Δτ == 0 => T = τ => M = 1
    M = Δτ > 0 ? trunc(Integer, (tN - τ) / Δτ) + 1 : 1
    windows = Vector{AbstractTimeParallelSolution}(undef, M)
    restarts = zeros(Integer, M)
    return MovingWindowSolution(windows, restarts)
end

#####
##### Functions
#####

"""
    length(solution::MovingWindowSolution)

returns the number of windows of `solution`.
"""
Base.length(solution::MovingWindowSolution) = length(solution.windows)

"""
    getindex(solution::MovingWindowSolution, m::Int)

returns the `m`-th window of `solution`.
"""
Base.getindex(solution::MovingWindowSolution, m::Integer) = solution.windows[m]

"""
    setindex!(solution::MovingWindowSolution, timeparallelsolution::AbstractTimeParallelSolution, m::Int)

stores a `timeparallelsolution` as the `m`-th window of `solution`.
"""
Base.setindex!(solution::MovingWindowSolution, timeparallelsolution::AbstractTimeParallelSolution, m::Integer) = solution.windows[m] = timeparallelsolution

"""
    firstindex(solution::MovingWindowSolution)

returns the first index of `solution`.
"""
Base.firstindex(solution::MovingWindowSolution) = firstindex(solution.windows)

"""
    lastindex(solution::MovingWindowSolution)

returns the last index of `solution`.
"""
Base.lastindex(solution::MovingWindowSolution) = lastindex(solution.windows)
