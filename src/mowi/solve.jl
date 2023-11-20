function (mowi::MoWi)(solution::MoWiSolution, problem::AbstractInitialValueProblem)
    @↓ parallelsolver, adaptive, τ, Δτ = mowi
    if adaptive !== nothing
        @↓ δτ⁻, δτ⁺, δΔτ⁻, δΔτ⁺, R, fineupdate = adaptive
    end

    @↓ finesolver, coarsolver, tolerance = parallelsolver
    @↓ N = parallelsolver.parameters

    @↓ hF ← h = finesolver.stepsize
    hF₀ = hF
    @↓ hG ← h = coarsolver.stepsize
    hG₀ = hG

    @↓ ϵ = tolerance

    @↓ windows, restarts = solution
    @↓ u0, (t0, tN) ← tspan = problem

    τ0 = t0
    τN = τ0 + τ
    τJ = τ0 # τJ = T_{N-ΔN}^{m-1} = T_{ΔN}^m

    # main loop
    m = 1
    M = M₀ = length(windows)

    cache = TimeParallel.TimeParallelCache(problem, parallelsolver)
    @↓ skips, U, T = cache

    flag = true
    itsnan = false
    while (τ0 < tN && τJ < tN) || itsnan

        if m > M
            append!(windows, Vector{AbstractTimeParallelSolution}(undef, M₀)) # TO-DO: Vector{typeof(windows[m])}(undef, M₀)
            append!(restarts, zeros(M₀)) # TO-DO: Vector{typeof(restarts[m])}(undef, M₀)
            M += M₀
        end

        # coarse run (serial)
        for n = 1:N+1
            T[n] = (N - n + 1) / N * τ0 + (n - 1) / N * τN # stable sum
        end
        if m == 1
            # initial guess
            U[1] = u0
            skips[1] = true
            for n = 1:N-1
                windowchunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
                windowchunksolution = coarsolver(windowchunkproblem)
                U[n+1] = windowchunksolution(T[n+1])
                skips[n+1] = true
            end
        elseif m > 1
            Uₙs = if restarts[m] ≠ 0
                windows[m] # i.e. previous restart
            else
                windows[m-1]
            end
            # shifting procedure
            for n = 1:N
                if T[n] ≤ τJ
                    U[n] = Uₙs(T[n])
                    if n == 1
                        skips[n] = true
                    else
                        skips[n] = false
                    end
                else
                    windowchunkproblem = subproblemof(problem, U[n-1], T[n-1], T[n])
                    U[n] = coarsolver(windowchunkproblem).u[end]
                    skips[n] = true
                end
            end
        end
        υ0 = U[1]

        # time-parallel subroutine (parallel)
        windowproblem = subproblemof(problem, υ0, τ0, τN)
        windows[m] = parallelsolver(cache, windowproblem)

        τJ = τN # update τJ for next iteration (next window or restart)

        # adaptive check
        if adaptive == nothing
            τ0 += Δτ
            m += 1
            flag = true
        else
            if (lasterror(windows[m]) > ϵ || isnan(lasterror(windows[m]))) && (restarts[m] < R)

                itsnan = isnan(lasterror(windows[m])) && (τJ > tN) ? true : false

                # δτ⁻₂ = max(δτ⁻, N * hG / τ) # smallest allowed is ΔT = hG
                τ *= δτ⁻
                Δτ *= δτ⁻
                if fineupdate
                    hF *= δτ⁻
                    @↑ finesolver.stepsize = h ← hF
                end
                hG = max(hG * δτ⁻, hF)
                @↑ coarsolver.stepsize = h ← hG

                # Δτ = max(Δτ * δΔτ⁻, τ / N)

                restarts[m] += 1
                flag = false
            else
                itsnan = false

                τ *= δτ⁺
                Δτ *= δτ⁺
                if fineupdate
                    hF *= δτ⁺
                    @↑ finesolver.stepsize = h ← hF
                end
                hG *= δτ⁺
                @↑ coarsolver.stepsize = h ← hG

                # Δτ = min(Δτ * δΔτ⁺, τ) # avoids overshoot

                τ0 = min(τ0 + Δτ, τJ)

                m += 1
                flag = true
            end
            # print("m = $m, Rₘ = $(restarts[m]), τ0 = $τ0, old τN = $τN, ")
        end

        τN = τ0 + τ
        # println("new τN = $τN")

    end

    if flag
        m -= 1 # to correct m += 1 at end of loop
    else
        # println("m = $m, Rₘ = $(restarts[m]), τ0 = $τ0, old τN = $τN, new τN = $τN")
    end

    resize!(windows, m)
    resize!(restarts, m)
    # @↑ solution = windows, restarts

    @↑ finesolver.stepsize = h ← hF₀
    @↑ coarsolver.stepsize = h ← hG₀

    return solution
end

function (mowi::MoWi)(problem::AbstractInitialValueProblem)
    solution = MoWiSolution(problem, mowi)
    mowi(solution, problem)
    return solution
end
