mutable struct MovingWindowSolver{𝒫_T, τ_T, Δτ_T} <: InitialValueSolver
    𝒫::𝒫_T
    τ::τ_T
    Δτ::Δτ_T
end

# function MovingWindowSolver(parallelsolver::TimeParallelSolver, τ, Δτ)
#     if Δτ < τ / parallelsolver.P
#         error("Select Δτ ≥ τ / P!")
#     elseif Δτ > τ
#         error("Select Δτ ≤ τ!")
#     end
#     new(parallelsolver, τ, Δτ)
# end

MovingWindowSolver(𝒫; τ, Δτ) = MovingWindowSolver(𝒫, τ, Δτ)
@doc (@doc MovingWindowSolver) MoWA(args...; kwargs...) = MovingWindowSolver(args...; kwargs...)
