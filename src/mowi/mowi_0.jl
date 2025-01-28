function mowi_0(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem, mowi::MoWi; kwargs...)
    
    @↓ parallelsolver, τ, Δτ = mowi
    @↓ finesolver, coarsesolver, parameters, tolerance = parallelsolver
    @↓ N = parameters
    @↓ ϵ = tolerance
    @↓ τ0, τN, τJ = cache
    @↓ windows = solution
    @↓ (t0, tN) ← tspan = problem

    # main loop
    m = 1
    while τN < tN

        if m > 1
            τ0 += Δτ
            τJ = τN
            τN = τ0 + τ
        end

        println("window # $m: τ0 = $τ0, τJ = $τJ, τN = $τN")

        windowproblem = copy(problem, τ0, τN) # u0 here is the wrong one; careful

        if m > 1
            cache = NSDETimeParallel.TimeParallelCache(windowproblem, parallelsolver)
            @↓ makeGs, Uₘ ← U, Tₘ ← T = cache

            for n = 1:N
                # between τ0 and τJ copy parallelsolution
                if τ0 ≤ Tₘ[n] ≤ τJ
                    Uₘ[n] = windows[m-1](Tₘ[n])
                    makeGs[n] = false
                end

                # between τJ and τN coarse guess inside parareal
                if τJ < Tₘ[n] ≤ τN
                    makeGs[n] = true
                end
            end
        end

        windows[m] = solve(windowproblem, parallelsolver; kwargs...)

        # move on
        m += 1
        
    end

end
