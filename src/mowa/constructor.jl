"""
    MovingWindowSolver <: AbstractMovingWindowSolver

A composite type for the moving-window algorithm.

# Constructors
```julia
MovingWindowSolver(parallelsolver, budget, adaptive, τ, Δτ)
MovingWindowSolver(parallelsolver; budget, adaptive, τ, Δτ)
MoWA(args...; kwargs...)
```

# Arguments
- `parallelsolver :: AbstractTimeParallelSolver`
- `adaptive :: Union{AbstractAdaptiveParameters, Nothing}`
- `budget :: Union{AbstractBudget, Nothing}` : budget available.
- `τ :: Real` : time-window length.
- `Δτ :: Real` : time-window shifting.
"""
mutable struct MovingWindowSolver{parallelsolver_T<:AbstractTimeParallelSolver, adaptive_T<:Union{AbstractAdaptiveParameters, Nothing}, budget_T<:Union{AbstractBudget, Nothing}, τ_T<:Real, Δτ_T<:Real} <: AbstractMovingWindowSolver
    parallelsolver::parallelsolver_T
    adaptive::adaptive_T
    budget::budget_T
    τ::τ_T
    Δτ::Δτ_T
end
function MovingWindowSolver(parallelsolver::AbstractTimeParallelSolver; adaptive::Union{AbstractAdaptiveParameters, Nothing}=nothing, budget::Union{AbstractBudget, Nothing}=nothing, τ::Real, Δτ::Real=τ)
    @↓ P = parallelsolver
    if Δτ < τ/P || Δτ > τ
        error("Select τ/P ≤ Δτ ≤ τ.")
    end
    return MovingWindowSolver(parallelsolver, adaptive, budget, τ, Δτ)
end
@doc (@doc MovingWindowSolver) MoWA(args...; kwargs...) = MovingWindowSolver(args...; kwargs...)
