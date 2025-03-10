struct MoWiCache{windowboundaries_T<:(Tuple{ℝ,ℝ,ℝ} where ℝ<:Real), windowproblem_T<:AbstractInitialValueProblem, windowcache_T<:NSDETimeParallel.AbstractTimeParallelCache} <: AbstractMovingWindowCache
    windowboundaries::windowboundaries_T
    windowproblem::windowproblem_T
    windowcache::windowcache_T
    # τ0::AbstractFloat
    # τN::AbstractFloat
    # τJ::AbstractFloat
end

function MoWiCache(problem::AbstractInitialValueProblem, mowi::MoWi)
    @↓ (t0, tN) ← tspan = problem
    @↓ parallelsolver, τ = mowi
    # τJ = T_{N-ΔN}^{m-1} = T_{ΔN}^m
    τJ = τ0 = t0
    τN = τ0 + τ
    windowboundaries = (τ0, τJ, τN)
    windowproblem = copy(problem, τ0, τN)
    windowcache = NSDETimeParallel.TimeParallelCache(windowproblem, parallelsolver)
    return MoWiCache(windowboundaries, windowproblem, windowcache)
end
