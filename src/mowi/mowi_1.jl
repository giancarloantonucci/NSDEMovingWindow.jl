function mowi_1(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem, mowi::MoWi; kwargs...)
    # Extract parameters
    @↓ windowboundaries, windowproblem, windowcache = cache
    τ0, τJ, τN = windowboundaries
    @↓ makeGs, Uₘ ← U, Tₘ ← T = windowcache
    @↓ parallelsolver, adaptive, τ, Δτ = mowi
    @↓ finesolver, coarsesolver, parameters, tolerance = parallelsolver
    @↓ N = parameters
    @↓ ϵ = tolerance
    @↓ δτ⁻, δτ⁺, R, fineupdate = adaptive
    @↓ windows, restarts = solution
    @↓ (t0, tN) ← tspan = problem

    # Initialize step sizes
    @↓ hF ← h = finesolver.stepsize
    # hF_tmp = hF
    @↓ hG ← h = coarsesolver.stepsize
    # hG_tmp = hG

    # Initialize window count
    M = M_tmp = length(windows)

    # Main-loop variables
    needsrestart = false
    τ0_tmp = τ0

    # Main loop
    m = 1
    while τN < tN
        # Resize storage if needed
        if m > M
            append!(windows, Vector{AbstractTimeParallelSolution}(undef, M_tmp)) # TODO: Vector{typeof(windows[m])}(undef, M_tmp)
            append!(restarts, zeros(M_tmp)) # TODO: Vector{typeof(restarts[m])}(undef, M_tmp)
            M += M_tmp
        end

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
                    if needsrestart
                        Uₘ[n] = windows[m](tₙ)
                    else
                        Uₘ[n] = windows[m-1](tₙ)
                    end
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
        if windows[m].errors[end] > ϵ && restarts[m] < R
            needsrestart = true
            restarts[m] += 1

            # Update adaptive paramters
            τ *= δτ⁻
            Δτ *= δτ⁻
            if fineupdate
                hF *= δτ⁻
                @↑ finesolver.stepsize = h ← hF
            end
            hG = max(hG * δτ⁻, hF)
            @↑ coarsesolver.stepsize = h ← hG
        else
            needsrestart = false
            m += 1

            # Update adaptive paramters
            τ *= δτ⁺
            Δτ *= δτ⁺
            if fineupdate
                hF *= δτ⁺
                @↑ finesolver.stepsize = h ← hF
            end
            hG *= δτ⁺
            @↑ coarsesolver.stepsize = h ← hG
        end
    end

    # Resize final storage
    resize!(windows, m-1)
    resize!(restarts, m-1)
end
