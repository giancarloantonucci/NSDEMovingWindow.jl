function (mowi::MoWi)(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem; isstrategy=0, kwargs...)
    if isstrategy == 0
        mowi_0(cache, solution, problem, mowi)
    elseif isstrategy == 1
        mowi_1(cache, solution, problem, mowi)
    elseif isstrategy == 2
        mowi_2(cache, solution, problem, mowi)
    # elseif isstrategy == 3
    #     mowi_3(cache, solution, problem, mowi)
    end
    return solution
end

function (mowi::MoWi)(solution::MoWiSolution, problem::AbstractInitialValueProblem; kwargs...)
    cache = MoWiCache(problem, mowi)
    mowi(cache, solution, problem; kwargs...)
    return solution
end

function (mowi::MoWi)(problem::AbstractInitialValueProblem; kwargs...)
    solution = MoWiSolution(problem, mowi)
    mowi(solution, problem; kwargs...)
    return solution
end
