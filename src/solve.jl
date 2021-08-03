
function NSDEBase.solve!(solution::MovingWindowSolution, problem, solver::MovingWindowSolver)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ 𝒫, τ, Δτ = solver
    @↓ 𝒢, P = 𝒫
    for m = 1:length(solution)
        solution[m] = TimeParallelSolution(problem, 𝒫)
        tmp = solution[m]
        @↓ U, T = tmp
        if m == 1
            coarseguess!(solution[m], problem, u0, t0, t0 + τ, 𝒫)
        else
            ΔP = trunc(Int, P * Δτ / τ)
            N = P - ΔP + 1
            for n = 1:length(T)
                T[n] = solution[m-1].T[n] + Δτ
            end
            for n = 1:N
                U[n] = solution[m-1].U[ΔP+n]
            end
            for n = N:P
                chunk = 𝒢(problem, U[n], T[n], T[n+1])
                U[n+1] = chunk.u[end]
            end
        end
        solve_serial!(solution[m], problem, 𝒫)
    end
    solution
end

function NSDEBase.solve(problem, solver::MovingWindowSolver)
    solution = MovingWindowSolution(problem, solver)
    solve!(solution, problem, solver)
    solution
end
