# function NSDEBase.cost(problem::AbstractInitialValueProblem, solver::MovingWindow.AbstractMovingWindowSolver, solution::MovingWindow.AbstractMovingWindowSolution)
#     @↓ (t0, tN) ← tspan = problem
#     @↓ τ, Δτ, parallelsolver, budget = solver
#     @↓ finesolver, coarsolver,  P = parallelsolver
#     @↓ hF ← h = finesolver.stepsize
#     @↓ hG ← h = coarsolver.stepsize
#     ξ = hG / hF
#     @↓ B ← budget, a₁, a₂, a₃ = budget
#     M = length(solution)
#     r = solution.restarts
#     α = τ / Δτ
#     return τ / abs(tN - t0) * (M*(M-1)/(2α*ξ) + (P+ξ)/(2P*ξ)*sum(1 + (B+1)*r[m] + length(solution[m].errors) for m=1:M))
# end
# cost(problem, mowa, mowasolution) * Cₛ
