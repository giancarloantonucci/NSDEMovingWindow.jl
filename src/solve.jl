"""
    solve!(solution::AbstractMovingWindowSolution, problem, solver::AbstractMovingWindowSolver; kwargs...) :: AbstractMovingWindowSolution

returns the [`AbstractMovingWindowSolution`](@ref) of an [`NSDEBase.AbstractInitialValueProblem`](@extref).
"""
function NSDEBase.solve!(solution::AbstractMovingWindowSolution, problem::AbstractInitialValueProblem, solver::AbstractMovingWindowSolver; kwargs...)
    return solver(solution, problem; kwargs...)
end

"""
    solve(problem, solver::AbstractMovingWindowSolver; kwargs...) :: AbstractMovingWindowSolution

returns the [`AbstractMovingWindowSolution`](@ref) of an [`NSDEBase.AbstractInitialValueProblem`](@extref).
"""
function NSDEBase.solve(problem::AbstractInitialValueProblem, solver::AbstractMovingWindowSolver; kwargs...)
    return solver(problem; kwargs...)
end
