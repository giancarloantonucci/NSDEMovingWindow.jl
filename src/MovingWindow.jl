module MovingWindow

export MovingWindowSolver
export MovingWindowSolution
export MoWA

using Reexport
using ArrowMacros
@reexport using TimeParallel
using RecipesBase

include("solver.jl")
include("solution.jl")
include("solve.jl")
include("plot.jl")

end
