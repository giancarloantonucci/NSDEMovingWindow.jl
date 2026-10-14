# Shared fixtures and helpers for this package test suite.
# Truth for accuracy assertions: the serial fine solve over the whole span.
const problem = Logistic(0.3, (0.0, 3.0))
const finesolver = RK4(; h=1e-3)
const finetruth = solve(problem, finesolver)

function good_parareal()
    return Parareal(
        finesolver,
        RungeKutta4(; h=5e-2);
        parameters=PararealParameters(; N=4, K=4),
        criterion=MeanRelativeUpdate(; ϵ=1e-9),
    )
end
# A setup built to FAIL windows: hopeless coarse solver, tiny budget, tight ϵ.
# ϵ sits above what a shrunken/loosened window reaches after a restart or two,
# so the adaptive machinery both FIRES and RECOVERS — a tolerance far below
# reach would send Stretch/Leap into a shrink spiral (Δτ collapsing, window
# count exploding), which is the algorithm as designed, not a test's place.
function bad_parareal(; criterion=MeanRelativeUpdate(; ϵ=1e-5))
    return Parareal(
        finesolver, Euler(; h=0.5); parameters=PararealParameters(; N=4, K=2), criterion
    )
end

windowstarts(sol) = [NSDEMovingWindow.windowspan(sol[m])[1] for m in 1:length(sol)]
function windowlengths(sol)
    return [
        let (ta, tb) = NSDEMovingWindow.windowspan(sol[m])
            tb - ta
        end for m in 1:length(sol)
    ]
end
