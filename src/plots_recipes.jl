@recipe function f(solution::MoWiSolution)
    @↓ windows = solution
    for m in 1:length(windows)
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

@recipe function f(wrapper::NSDEBase._PhasePlot{<:MoWiSolution})
    @↓ solution ← plottable = wrapper
    @↓ windows = solution
    for m in 1:length(windows)
        @series begin
            m == 1 || (label := "")
            NSDEBase._PhasePlot(solution[m])
        end
    end
end

@recipe function f(wrapper::NSDEBase._Convergence{<:MoWiSolution})
    gridalpha --> 0.2
    markershape --> :circle
    markerstrokewidth --> 0
    seriestype --> :path
    xticks --> 0:1000
    yticks --> 0:1000
    @↓ solution ← plottable = wrapper
    @↓ windows, restarts = solution
    given = get(plotattributes, :label, nothing)
    labelfor(i) =
        if given === nothing
            (i == 1 ? "iterations" : "restarts")
        elseif given isa Union{AbstractString,Symbol}
            (i == 1 ? string(given) : "")
        else
            (i <= length(given) ? given[i] : "")
        end
    @series begin
        label := labelfor(1)
        [length(solution[m].errors) for m in 1:length(windows)]
    end
    @series begin
        label := labelfor(2)
        restarts
    end
end
