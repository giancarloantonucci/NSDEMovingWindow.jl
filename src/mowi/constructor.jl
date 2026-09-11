# src/mowi/constructor.jl

"""
    MoWi <: AbstractMovingWindowSolver

The moving-window algorithm (thesis, Ch. 4): long-time integration through a
sequence of overlapping windows, each solved by a time-parallel subroutine
seeded from the previous window's solution over the overlap. The output is an
ENSEMBLE of short window solutions, not one long trajectory — in chaotic
systems each shift jumps orbit by design; see [`ensemblemean`](@ref) and
friends for the statistics this enables.

# Constructors
```julia
MoWi(parallelsolver; adaptive=nothing, τ, Δτ=τ)
```

# Arguments
- `parallelsolver :: AbstractTimeParallelSolver` (deep-copied: the adaptive
  strategies mutate solver step sizes and weights during a run)
- `adaptive` : `nothing` for fixed windows, or one of [`StretchParameters`](@ref)
  (AMoWi-1), [`LeapParameters`](@ref) (AMoWi-2), [`ZoomParameters`](@ref)
  (AMoWi-3) — the parameter type selects the strategy
- `τ :: Real` : window length, `τ > 0`
- `Δτ :: Real` : window shift, `0 < Δτ ≤ τ`
"""
mutable struct MoWi{parallelsolver_T<:AbstractTimeParallelSolver, adaptive_T<:Union{AbstractAdaptiveMoWiParameters, Nothing}, τ_T<:Real, Δτ_T<:Real} <: AbstractMovingWindowSolver
    parallelsolver::parallelsolver_T
    adaptive::adaptive_T
    τ::τ_T
    Δτ::Δτ_T
    function MoWi(parallelsolver::parallelsolver_T, adaptive::adaptive_T, τ::τ_T, Δτ::Δτ_T) where {parallelsolver_T<:AbstractTimeParallelSolver, adaptive_T<:Union{AbstractAdaptiveMoWiParameters, Nothing}, τ_T<:Real, Δτ_T<:Real}
        τ > 0 || throw(ArgumentError("MoWi needs a window length τ > 0, got τ = $τ."))
        0 < Δτ ≤ τ || throw(ArgumentError("MoWi needs a window shift 0 < Δτ ≤ τ, got Δτ = $Δτ with τ = $τ."))
        if adaptive isa ZoomParameters && parallelsolver.tolerance.ψ === NSDETimeParallel.ψ₁
            throw(ArgumentError("ZoomParameters (AMoWi-3) adapts the ψ weights, but the parallel solver's tolerance uses the unweighted ψ₁ — zooming would be a silent no-op. Build the solver's Tolerance with the weighted error function ψ₂."))
        end
        if adaptive isa ZoomParameters && parallelsolver.tolerance.weights.updatew
            throw(ArgumentError("ZoomParameters (AMoWi-3) adapts the base weight w, but Weights(updatew = true) re-measures w every iteration and floors it at the measured rate — the zoom-in direction (δw⁺) becomes a silent no-op. Freeze the weight: Weights(w = exp(Λ)) with Λ known or measured by a probe solve."))
        end
        return new{parallelsolver_T, adaptive_T, τ_T, Δτ_T}(deepcopy(parallelsolver), adaptive, τ, Δτ)
    end
end

MoWi(parallelsolver::AbstractTimeParallelSolver; adaptive::Union{AbstractAdaptiveMoWiParameters, Nothing}=nothing, τ::Real, Δτ::Real=τ) = MoWi(parallelsolver, adaptive, τ, Δτ)
