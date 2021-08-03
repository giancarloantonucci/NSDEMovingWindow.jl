mutable struct MovingWindowSolution{windows_T} <: InitialValueSolution
    windows::windows_T
end

function MovingWindowSolution(problem, solver::MovingWindowSolver)
    @↓ (t0, tN) ← tspan = problem
    @↓ τ, Δτ = solver
    M = 1 # if solver.Δτ == 0 then T = τ and M = 1
    if Δτ > 0
        M = trunc(Int, (tN - τ) / Δτ) + 1
    end
    windows = Vector{TimeParallelSolution}(undef, M)
    return MovingWindowSolution(windows)
end

Base.length(solution::MovingWindowSolution) = length(solution.windows)
Base.getindex(solution::MovingWindowSolution, m::Int) = solution.windows[m] # solution[m] ≡ solution.windows[m]
Base.setindex!(solution::MovingWindowSolution, value::TimeParallelSolution, m::Int) = solution.windows[m] = value
Base.lastindex(solution::MovingWindowSolution) = lastindex(solution.windows)
