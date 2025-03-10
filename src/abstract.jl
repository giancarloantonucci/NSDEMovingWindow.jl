"An abstract type for the [`MoWi`](@ref) algorithm."
abstract type AbstractMovingWindowSolver <: AbstractInitialValueSolver end

"An abstract type for solutions computed with [`MoWi`](@ref)."
abstract type AbstractMovingWindowSolution <: AbstractInitialValueSolution end

"An abstract type for caching intermediate computations in [`MoWi`](@ref)."
abstract type AbstractMovingWindowCache <: AbstractInitialValueCache end

"An abstract type for parameters in [`MoWi`](@ref)."
abstract type AbstractMovingWindowParameters <: AbstractInitialValueParameters end

"An abstract type for adaptive parameters in [`MoWi`](@ref)."
abstract type AbstractAdaptiveMoWiParameters <: AbstractMovingWindowParameters end
