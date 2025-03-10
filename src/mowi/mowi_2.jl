function mowi_2(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem, mowi::MoWi; kwargs...)
    # Extract parameters
    @↓ windowboundaries, windowproblem, windowcache = cache
    τ0, τJ, τN = windowboundaries
    @↓ makeGs, Uₘ ← U, Tₘ ← T = windowcache
    @↓ parallelsolver, adaptive, τ, Δτ = mowi
    @↓ finesolver, coarsesolver, parameters, tolerance = parallelsolver
    @↓ N = parameters
    @↓ ϵ = tolerance
    @↓ δΔτ⁻, δΔτ⁺, R = adaptive
    @↓ windows, restarts = solution
    @↓ (t0, tN) ← tspan = problem

    # Solve the first window problem (m = 1)
    println("window # 1: τ0 = $τ0, τJ = $τJ, τN = $τN")
    windows[1] = parallelsolver(windowproblem; kwargs...)

    # Initialize window count
    M = M_tmp = length(windows)

    # Main-loop variables
    needsrestart = false
    τ0_tmp = τ0

    # Main loop
    m = 2
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
        else
            τ0 = τ0_tmp + Δτ
            τJ = τN
            τN = τ0 + τ
        end
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
                Uₘ[n] = windows[m-1](Tₘ[n])
                makeGs[n] = false
            elseif τJ < Tₘ[n] ≤ τN # between τJ and τN coarse guess inside parareal
                makeGs[n] = true
            end
        end

        windows[m] = parallelsolver(windowcache, windowproblem; kwargs...)

        if windows[m].errors[end] > ϵ && restarts[m] < R
            needsrestart = true
            restarts[m] += 1

            # Update adaptive paramters
            Δτ *= δΔτ⁻
        else
            needsrestart = false
            τ0_tmp = τ0
            m += 1

            # Update adaptive paramters
            Δτ = min(Δτ * δΔτ⁺, τ)
        end
    end

    # Resize final storage
    resize!(windows, m-1)
    resize!(restarts, m-1)
end
