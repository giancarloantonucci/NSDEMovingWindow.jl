# src/mowi/solve.jl

function (mowi::MoWi)(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem; kwargs...)
    return mowi!(cache, solution, problem, mowi; kwargs...)
end

function (mowi::MoWi)(solution::MoWiSolution, problem::AbstractInitialValueProblem; kwargs...)
    cache = MoWiCache(problem, mowi)
    return mowi(cache, solution, problem; kwargs...)
end

function (mowi::MoWi)(problem::AbstractInitialValueProblem; kwargs...)
    solution = MoWiSolution(problem, mowi)
    return mowi(solution, problem; kwargs...)
end
