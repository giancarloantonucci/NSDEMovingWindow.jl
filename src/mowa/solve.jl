function (mowa::MovingWindowSolver)(solution::MovingWindowSolution, problem::AbstractInitialValueProblem; saveiterates::Bool=false)
    @↓ parallelsolver, budget, adaptive, τ, Δτ = mowa
    @↓ finesolver, coarsolver, tolerance, P, K = parallelsolver
    K₀ = K
    @↓ finestepsize ← stepsize = finesolver
    @↓ hF ← h = finestepsize
    hF₀ = hF
    @↓ coarstepsize ← stepsize = coarsolver
    @↓ hG ← h = coarstepsize
    hG₀ = hG
    @↓ ϵ = tolerance
    if budget !== nothing
        @↓ Kmax = budget
        @↑ parallelsolver = K ← Kmax
    end
    if adaptive !== nothing
        @↓ δτ⁻, δτ⁺, δΔτ⁻, δΔτ⁺, Nᵣ, fineupdate = adaptive
    end
    @↓ windows, restarts = solution
    @↓ u0, (t0, tN) ← tspan = problem
    
    τ0 = t0
    τN = τ0 + τ
    τJ = zero(τ0)
    # τ0, τN, τJ = windowparams(problem)
    
    cache = TimeParallel.TimeParallelCache(problem, parallelsolver)
    @↓ U, G, T = cache
    
    m = 1
    M = M₀ = length(windows)
    counter = counter2 = 0
    while τ0 < tN
        if τN ≥ tN
            counter += 1
        end
        # if still not done, δτ⁺ease solution.windows length
        if m > M
            append!(windows, Vector{AbstractTimeParallelSolution}(undef, M₀))
            append!(restarts, zeros(M₀))
            M += M₀
        end
        # start with a coarse guess
        if m == 1
            # define windowproblem
            windowproblem = subproblemof(problem, u0, τ0, τN)
            # update cache
            coarseguess!(cache, windowproblem, parallelsolver)
        elseif m > 1
            𝑢 = restarts[m] == 0 ? windows[m-1] : windows[m]
            # recompute T
            for n = 1:P+1
                T[n] = (P - n + 1) / P * τ0 + (n - 1) / P * τN
            end
            # recompute U
            for n = 1:P
                # transfer initial conditions from (m-1)-th window (or m-th window at previous "restart")
                if T[n] ≤ τJ
                    U[n] = 𝑢(T[n])
                    # compute G for Parareal
                    if n == 1
                        G[n] = U[n]
                    else
                        windowchunkproblem = subproblemof(problem, U[n-1], T[n-1], T[n])
                        G[n] = coarsolver(windowchunkproblem).u[end]
                    end
                # compute the remaining initial conditions with coarse solver
                else
                    windowchunkproblem = subproblemof(problem, U[n-1], T[n-1], T[n])
                    G[n] = U[n] = coarsolver(windowchunkproblem).u[end]
                end
            end
            windowproblem = subproblemof(problem, U[1], τ0, τN)
            # @↑ cache = U, G, T
        end
        windows[m] = parallelsolver(cache, windowproblem; saveiterates=saveiterates)
        # println("$(restarts[m]) restarts, $(length(solution[m].errors)) Parareal iterations")
        # println("——————————————————————————————————")
        # save τJ for next window
        τJ = τN
        # adaptive check
        if adaptive !== nothing
            if (windows[m].errors[end] > ϵ) && (restarts[m] < Nᵣ) && (hG > hF)
                # δτ⁻₂ = max(δτ⁻, P * hG / τ) # P * hG is minimum possible
                # τ *= δτ⁻₂
                # Δτ *= δτ⁻₂
                # if fineupdate
                #     hF *= δτ⁻₂
                # end
                # hG = max(hG * δτ⁻₂, hF)
                Δτ = max(δΔτ⁻ * Δτ, τ / P)
                restarts[m] += 1
                counter2 += 1
            else
                counter2 = 0
                # τ *= δτ⁺
                # Δτ *= δτ⁺
                # if fineupdate
                #     hF *= δτ⁺
                # end
                # hG *= δτ⁺
                Δτ = min(δΔτ⁺ * Δτ, τ) # avoids overshoot
                τ0 = min(τ0 + Δτ, τJ)
                m += 1
            end
            if fineupdate
                @↑ finestepsize = h ← hF
            end
            @↑ coarstepsize = h ← hG
        else
            τ0 += Δτ
            m += 1
        end
        τN = τ0 + τ
        if counter ≥ 1 && counter2 == 0
            break
        end
    end
    m -= 1 # to correct m += 1 at end of loop
    @↑ parallelsolver = K ← K₀
    @↑ finestepsize = h ← hF₀
    @↑ coarstepsize = h ← hG₀
    resize!(windows, m)
    resize!(restarts, m)
    return solution
end

function (mowa::MovingWindowSolver)(problem::AbstractInitialValueProblem; saveiterates::Bool=false)
    solution = MovingWindowSolution(problem, mowa)
    mowa(solution, problem; saveiterates=saveiterates)
    return solution
end
