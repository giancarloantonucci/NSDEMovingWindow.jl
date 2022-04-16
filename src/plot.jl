@recipe function f(solution::MovingWindowSolution)
    @↓ windows = solution
    for m = 1:length(windows)
        @series begin
            if m != 1
                label := ""
            else
                if haskey(plotattributes, :label)
                    label := plotattributes[:label]
                end
            end
            solution[m]
        end
    end
end

@recipe function f(wrappedobject::NSDEBase._PhasePlot{<:MovingWindowSolution})
    @↓ solution ← object = wrappedobject
    @↓ windows = solution
    for m = 1:length(windows)
        @series NSDEBase._PhasePlot(solution[m])
    end
end

@recipe function f(wrappedobject::NSDEBase._Convergence{<:MovingWindowSolution})
    gridalpha         --> 0.2
    markershape       --> :circle
    markerstrokewidth --> 0
    seriestype        --> :path
    xticks            --> 0:1000
    yticks            --> 0:1000
    @↓ solution ← object = wrappedobject
    @↓ windows, restarts = solution
    @series begin
        if haskey(plotattributes, :label)
            label := plotattributes[:label][1]
        end
        [length(solution[m].errors) for m = 1:length(windows)]
    end
    @series begin
        if haskey(plotattributes, :label)
            label := plotattributes[:label][2]
        end
        restarts
    end
end
