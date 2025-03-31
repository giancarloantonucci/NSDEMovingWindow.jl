function mowi_3(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem, mowi::MoWi; kwargs...)
    # Extract parameters
    @↓ windowboundaries, windowproblem, windowcache = cache
    τ0, τJ, τN = windowboundaries
    @↓ makeGs, Uₘ ← U, Tₘ ← T = windowcache
    @↓ parallelsolver, adaptive, τ, Δτ = mowi
    @↓ finesolver, coarsesolver, parameters, tolerance = parallelsolver
    @↓ N = parameters
    @↓ ϵ, weights = tolerance

    @↓ δw⁻, δw⁺, R, fineupdate = adaptive
    @↓ windows, restarts = solution
    @↓ (t0, tN) ← tspan = problem

    # Main-loop variables
    needsrestart = false
    τ0_tmp = τ0

    # Main loop
    m = 1
    stillnan = false
    while τN < tN || stillnan

        # Update window boundaries
        if needsrestart
            τ0 = τ0_tmp
            τN = τ0 + τ
        elseif !needsrestart && m > 1
            τ0 = τ0_tmp + Δτ
            τJ = τN
            τN = τ0 + τ
        end
        println("window # $m: τ0 = $τ0, τJ = $τJ, τN = $τN")

        # Update window problem
        tspan = (τ0, τN)
        @↑ windowproblem = tspan

        if m > 1 || needsrestart
            # Compute Tₘ
            for n = 1:N+1
                Tₘ[n] = (N - n + 1) / N * τ0 + (n - 1) / N * τN # stable sum
            end

            # Compute Uₘ
            for n = 1:N
                tₙ = Tₘ[n]
                if τ0 ≤ tₙ ≤ τJ # between τ0 and τJ copy parallelsolution
                    Uₘ[n] = windows[m-1](tₙ)
                    makeGs[n] = false
                elseif τJ < Tₘ[n] ≤ τN # between τJ and τN coarse guess inside parareal
                    makeGs[n] = true
                end
            end
        end

        # Solve window problem
        if m == 1 && !needsrestart
            windows[m] = parallelsolver(windowproblem; kwargs...)
        else
            windows[m] = parallelsolver(windowcache, windowproblem; kwargs...)
        end

        # Check error and decide whether to restart
        τ0_tmp = τ0

        if (windows[m].errors[end] > ϵ || isnan(windows[m].errors[end])) && restarts[m] < R
            if isnan(windows[m].errors[end])
                stillnan = true
            else
                stillnan = false
            end

            needsrestart = true
            restarts[m] += 1

            # Update adaptive paramters
            weights.δ = δw⁻
        else
            stillnan = false
            needsrestart = false
            m += 1

            # Update adaptive paramters
            weights.δ = δw⁺
        end
    end

    # Resize final storage
    resize!(windows, m-1)
    resize!(restarts, m-1)
end
