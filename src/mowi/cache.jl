# src/mowi/cache.jl

"""
    MoWiCache <: AbstractMovingWindowCache

One window problem and one time-parallel cache, reused for every window: the
skeleton re-targets them with [`NSDETimeParallel.shiftwindow!`](@ref) as the
window moves.
"""
struct MoWiCache{windowproblem_T<:AbstractInitialValueProblem, windowcache_T<:NSDETimeParallel.AbstractTimeParallelCache} <: AbstractMovingWindowCache
    windowproblem::windowproblem_T
    windowcache::windowcache_T
end

function MoWiCache(problem::AbstractInitialValueProblem, mowi::MoWi)
    @↓ (t0, tN) ← tspan = problem
    @↓ parallelsolver, τ = mowi
    windowproblem = copy(problem, t0, t0 + τ)
    windowcache = NSDETimeParallel.TimeParallelCache(windowproblem, parallelsolver)
    return MoWiCache(windowproblem, windowcache)
end
