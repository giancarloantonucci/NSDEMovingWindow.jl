function mowi_2(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem, mowi::MoWi; kwargs...)
    
    @↓ parallelsolver, adaptive, τ, Δτ = mowi
    @↓ finesolver, coarsesolver, parameters, tolerance = parallelsolver
    @↓ N = parameters
    @↓ ϵ = tolerance
    @↓ δΔτ⁻, δΔτ⁺ = adaptive
    @↓ τ0, τN, τJ = cache
    @↓ windows, restarts = solution
    @↓ (t0, tN) ← tspan = problem

    M = M₀ = length(windows)

    # main loop
    needsrestart = false
    τ0_old = τ0
    m = 1
    while τN < tN

        if m > M
            append!(windows, Vector{AbstractTimeParallelSolution}(undef, M₀)) # TODO: Vector{typeof(windows[m])}(undef, M₀)
            append!(restarts, zeros(M₀)) # TODO: Vector{typeof(restarts[m])}(undef, M₀)
            M += M₀
        end

        if needsrestart
            τ0 = τ0_old + Δτ
            τN = τ0 + τ
        end

        if m > 1 && !needsrestart
            τ0 = τ0_old + Δτ
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

        if windows[m].errors[end] > ϵ && restarts[m] < R
            needsrestart = true

            Δτ *= δΔτ⁻

            restarts[m] += 1
        else
            needsrestart = false
            τ0_old = τ0

            Δτ *= δΔτ⁺

            # move on
            m += 1
        end

    end

    resize!(windows, m-1)
    resize!(restarts, m-1)
    
end
