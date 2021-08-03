@recipe function f(solution::MovingWindowSolution; vars = nothing, label = "")
    M = length(solution)
    for m in 1:M
        vars    --> vars
        @series solution[m]
    end
end
