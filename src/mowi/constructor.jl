"""
    MoWi <: AbstractMovingWindowSolver

The moving-window algorithm: long-time integration through a sequence of
overlapping windows, each solved by a time-parallel subroutine seeded from the
previous window's solution over the overlap. The output is an ensemble of short
window solutions, not one long trajectory: in chaotic systems each shift may
jump orbit. See [`ensemblemean`](@ref) and friends for the statistics this
enables.

The requested final time is a hard boundary: no window is integrated beyond it.
A genuine final truncation is recorded in `solution.clipped` and excluded from
ensemble statistics unless `include_clipped=true`. Gaps within the endpoint
snap tolerance are treated as rounding and do not mark a window as clipped.

# Constructors
```julia
MoWi(parallelsolver; adaptive=nothing, τ, Δτ=τ)
```

# Arguments
- `parallelsolver :: AbstractTimeParallelSolver` : copied for isolation;
  adaptive strategies replace settings in run-local state
- `adaptive` : `nothing` for fixed windows, or one of [`StretchParameters`](@ref)
  (AMoWi-1), [`LeapParameters`](@ref) (AMoWi-2), [`ZoomParameters`](@ref)
  (AMoWi-3); the parameter type selects the strategy
- `τ :: Real` : window length, `τ > 0`
- `Δτ :: Real` : window shift, `0 < Δτ ≤ τ`
"""
struct MoWi{
    parallelsolver_T<:AbstractTimeParallelSolver,
    adaptive_T<:Union{AbstractAdaptiveMoWiParameters,Nothing},
    τ_T<:Real,
    Δτ_T<:Real,
} <: AbstractMovingWindowSolver
    parallelsolver::parallelsolver_T
    adaptive::adaptive_T
    τ::τ_T
    Δτ::Δτ_T
    function MoWi(
        parallelsolver::parallelsolver_T, adaptive::adaptive_T, τ::τ_T, Δτ::Δτ_T
    ) where {
        parallelsolver_T<:AbstractTimeParallelSolver,
        adaptive_T<:Union{AbstractAdaptiveMoWiParameters,Nothing},
        τ_T<:Real,
        Δτ_T<:Real,
    }
        τ > 0 || throw(ArgumentError("MoWi needs a window length τ > 0, got τ = $τ."))
        0 < Δτ ≤ τ ||
            throw(ArgumentError("MoWi needs a window shift 0 < Δτ ≤ τ, got Δτ = $Δτ with τ = $τ."))
        if adaptive isa StretchParameters
            g = parallelsolver.coarsesolver
            hasproperty(g, :stepsize) && applicable(withstepsize, g, g.stepsize.h) || throw(
                ArgumentError("Stretch requires a coarse solver with stepsize and withstepsize")
            )
            f = parallelsolver.finesolver
            if adaptive.fineupdate
                hasproperty(f, :stepsize) && applicable(withstepsize, f, f.stepsize.h) || throw(
                    ArgumentError(
                        "Stretch(fineupdate=true) requires a fine solver with stepsize and withstepsize",
                    ),
                )
            end
        end
        if adaptive isa ZoomParameters
            c = parallelsolver.criterion
            c isa NSDETimeParallel.WeightedDefect ||
                throw(ArgumentError("Zoom requires WeightedDefect"))
            c.weight isa NSDETimeParallel.FixedWeight || throw(
                ArgumentError(
                    "Zoom controls a fixed weight; secant and neighbour updates are incompatible",
                ),
            )
        end
        return new{parallelsolver_T,adaptive_T,τ_T,Δτ_T}(deepcopy(parallelsolver), adaptive, τ, Δτ)
    end
end

function MoWi(
    parallelsolver::AbstractTimeParallelSolver;
    adaptive::Union{AbstractAdaptiveMoWiParameters,Nothing}=nothing,
    τ::Real,
    Δτ::Real=τ,
)
    return MoWi(parallelsolver, adaptive, τ, Δτ)
end
