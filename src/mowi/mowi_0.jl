function mowi_0(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem, mowi::MoWi; kwargs...)
    # Extract parameters
    @↓ windowboundaries, windowproblem, windowcache = cache
    τ0, τJ, τN = windowboundaries
    @↓ makeGs, Uₘ ← U, Tₘ ← T = windowcache
    @↓ parallelsolver, τ, Δτ = mowi
    @↓ finesolver, coarsesolver, parameters, tolerance = parallelsolver
    @↓ N = parameters
    @↓ ϵ = tolerance
    @↓ windows = solution
    @↓ u0, (t0, tN) ← tspan = problem

    # Solve the first window problem (m = 1)
    println("window # 1: τ0 = $τ0, τJ = $τJ, τN = $τN")
    windows[1] = parallelsolver(windowproblem; kwargs...)

    # Main loop (m > 1)
    m = 2
    while τN < tN
        # Update window boundaries
        τ0 += Δτ
        τJ = τN
        τN = τ0 + τ
        println("window # $m: τ0 = $τ0, τJ = $τJ, τN = $τN")

        # Update window problem
        tspan = (τ0, τN)
        @↑ windowproblem = tspan

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

        # Solve the window problem
        windows[m] = parallelsolver(windowcache, windowproblem; kwargs...)

        m += 1
    end
end
