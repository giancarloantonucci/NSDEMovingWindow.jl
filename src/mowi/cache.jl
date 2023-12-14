struct MoWiCache <: AbstractMovingWindowCache
    τ0::AbstractFloat
    τN::AbstractFloat
    τJ::AbstractFloat
end

function MoWiCache(problem::AbstractInitialValueProblem, mowi::MoWi)
    @↓ t0 ← tspan[1] = problem
    τ0 = t0
    @↓ τ = mowi
    τN = τ0 + τ
    τJ = τ0 # τJ = T_{N-ΔN}^{m-1} = T_{ΔN}^m
    return MoWiCache(τ0, τN, τJ)
end