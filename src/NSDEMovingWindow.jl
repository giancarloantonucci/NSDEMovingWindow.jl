module NSDEMovingWindow

using Reexport
using NSDEBase
using ArrowMacros
@reexport using NSDETimeParallel
using RecipesBase

include("abstract.jl")
include("mowi/adaptive.jl")
include("mowi/constructor.jl")
include("mowi/cache.jl")
include("mowi/solution.jl")
include("mowi/mowi.jl")
include("mowi/solve.jl")
include("solve.jl")
include("plots_recipes.jl")

export AbstractMovingWindowSolver
export AbstractMovingWindowSolution
export AbstractMovingWindowCache
export AbstractMovingWindowParameters

export MoWi
export MoWiSolution
export StretchParameters, LeapParameters, ZoomParameters

export MovingWindowSolution
export ensemblemean, ensemblevariance, ensemblesem

end
