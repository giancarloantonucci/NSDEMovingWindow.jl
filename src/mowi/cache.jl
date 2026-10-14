"""
    MoWiCache <: AbstractMovingWindowCache

One window problem and one time-parallel cache, reused for every window: the
window loop retargets them with `NSDETimeParallel.shiftwindow!` as the window
moves.
"""
struct MoWiCache{
    windowproblem_T<:AbstractInitialValueProblem,
    windowcache_T<:NSDETimeParallel.AbstractTimeParallelCache,
} <: AbstractMovingWindowCache
    windowproblem::windowproblem_T
    windowcache::windowcache_T
end

# Window times work in at least Float64, like the chunk grid of
# NSDETimeParallel: Float32 spans and window lengths are widened, and BigFloat
# keeps its precision. The window problem, the loop state and the chunk grid
# then agree in one number type.
function windowtime(problem::AbstractInitialValueProblem, mowi::MoWi)
    return promote_type(Float64, map(t -> typeof(float(t)), (problem.tspan..., mowi.τ, mowi.Δτ))...)
end

function MoWiCache(problem::AbstractInitialValueProblem, mowi::MoWi)
    @↓ (t0, tN) ← tspan = problem
    @↓ parallelsolver, τ = mowi
    W = windowtime(problem, mowi)
    windowproblem = remake(problem; tspan=(W(t0), W(t0) + W(τ)))
    windowcache = NSDETimeParallel.TimeParallelCache(
        windowproblem, parallelsolver; duration=τ / parallelsolver.parameters.N
    )
    return MoWiCache(windowproblem, windowcache)
end
