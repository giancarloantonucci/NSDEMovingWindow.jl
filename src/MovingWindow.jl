module MovingWindow

using Reexport
using ArrowMacros
@reexport using TimeParallel
using RecipesBase

include("abstract.jl")
include("mowa/budget.jl")
include("mowa/adaptive.jl")
include("mowa/constructor.jl")
include("mowa/solution.jl")
include("mowa/solve.jl")
include("solve.jl")
include("plot.jl")

export MovingWindowSolver
export MovingWindowSolution
export MoWA

export Budget

end
