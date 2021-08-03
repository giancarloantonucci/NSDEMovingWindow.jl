@recipe function f(solution::MovingWindowSolution; vars = nothing, label = "")
    fontfamily --> "Computer Modern"
    framestyle --> :box
    gridalpha --> 0.2
    linewidth --> 1.5
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path
    tick_direction --> :out
    xwiden --> false
    M = length(solution)
    for m in 1:M
        vars --> vars
        @series begin
            seriescolor --> m
            solution[m]
        end
        @series begin
            seriescolor --> m
            seriestype := :vline
            linestyle --> :dash
            [solution[m][end][end].t[end]]
        end
    end
end
