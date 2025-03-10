"""
    MoWi <: AbstractMovingWindowSolver

A composite type for the moving-window algorithm.

# Constructors
```julia
MoWi(parallelsolver, adaptive, τ, Δτ)
MoWi(parallelsolver; adaptive, τ, Δτ)
```

# Arguments
- `parallelsolver :: AbstractTimeParallelSolver`
- `adaptive :: Union{AbstractAdaptiveMoWiParameters, Nothing}`
- `τ :: Real` : time-window length.
- `Δτ :: Real` : time-window shift.
"""
mutable struct MoWi{parallelsolver_T<:AbstractTimeParallelSolver, adaptive_T<:Union{AbstractAdaptiveMoWiParameters, Nothing}, τ_T<:Real, Δτ_T<:Real} <: AbstractMovingWindowSolver
    parallelsolver::parallelsolver_T
    adaptive::adaptive_T
    τ::τ_T
    Δτ::Δτ_T
    function MoWi(parallelsolver::parallelsolver_T, adaptive::adaptive_T, τ::τ_T, Δτ::Δτ_T) where {parallelsolver_T<:AbstractTimeParallelSolver, adaptive_T<:Union{AbstractAdaptiveMoWiParameters, Nothing}, τ_T<:Real, Δτ_T<:Real}
        parallelsolver2 = deepcopy(parallelsolver)
        return new{parallelsolver_T, adaptive_T, τ_T, Δτ_T}(parallelsolver2, adaptive, τ, Δτ)
    end
end

function MoWi(parallelsolver::AbstractTimeParallelSolver; adaptive::Union{AbstractAdaptiveMoWiParameters, Nothing}=nothing, τ::Real, Δτ::Real=τ)
    @↓ N = parallelsolver.parameters
    if Δτ > τ
        error("τ = $τ and Δτ = $(Δτ). Please, select Δτ ≤ τ.")
    end
    return MoWi(parallelsolver, adaptive, τ, Δτ)
end
