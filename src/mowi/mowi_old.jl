function (mowi::MoWi)(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem; mode::String="SERIAL")
    @↓ τ0, τN, τJ = cache

    @↓ parallelsolver, adaptive, τ, Δτ = mowi
    if adaptive !== nothing
        @↓ δτ⁻, δτ⁺, δΔτ⁻, δΔτ⁺, R, fineupdate = adaptive
    end

    @↓ finesolver, coarsesolver, tolerance = parallelsolver
    @↓ N = parallelsolver.parameters

    @↓ hF ← h = finesolver.stepsize
    @↓ hG ← h = coarsesolver.stepsize
    hF₀ = hF
    hG₀ = hG

    @↓ ϵ = tolerance

    @↓ windows, restarts = solution
    @↓ u0, (t0, tN) ← tspan = problem

    # main loop
    m = 1
    M = M₀ = length(windows)

    parallelcache = NSDETimeParallel.TimeParallelCache(problem, parallelsolver)
    @↓ skips = parallelcache
    parallelsolution = NSDETimeParallel.TimeParallelSolution(problem, parallelsolver)
    @↓ U, T = parallelsolution.lastiterate

    flag = true
    itsnan = false
    
    while (τ0 < tN && τJ < tN) || itsnan

        if m > M
            append!(windows, Vector{AbstractTimeParallelSolution}(undef, M₀)) # TODO: Vector{typeof(windows[m])}(undef, M₀)
            append!(restarts, zeros(M₀)) # TODO: Vector{typeof(restarts[m])}(undef, M₀)
            M += M₀
        end

        # coarse run (serial)
        for n = 1:N+1
            T[n] = (N - n + 1) / N * τ0 + (n - 1) / N * τN # stable sum
        end
        if m == 1
            # initial guess
            U[1] = u0
            for n = 1:N-1
                windowchunkproblem = copy(problem, U[n], T[n], T[n+1])
                windowchunksolution = coarsesolver(windowchunkproblem)
                U[n+1] = windowchunksolution(T[n+1])
            end
            skips .= true
        elseif m > 1
            Uₙs = if restarts[m] ≠ 0
                windows[m] # i.e. previous restart
            else
                windows[m-1]
            end
            # shifting procedure
            for n = 1:N
                if T[n] ≤ τJ
                    # if mode == "SERIAL"
                    #     U[n] = Uₙs(T[n])
                    # elseif mode == "DISTRIBUTED"

                    # end

                    # println(Uₙs.lastiterate.T[n], " ", T[n], " ", findfirst(Uₙs.lastiterate.T .≈ T[n]))
                    # U[n] = Uₙs.lastiterate.U[n]

                    U[n] = Uₙs(T[n])
                    if n == 1
                        skips[n] = true
                    else
                        skips[n] = false
                    end
                else
                    windowchunkproblem = copy(problem, U[n-1], T[n-1], T[n])
                    windowchunksolution = coarsesolver(windowchunkproblem)
                    U[n] = windowchunksolution(T[n])
                    skips[n] = true
                end
            end
        end
        υ0 = U[1]

        # time-parallel subroutine (parallel)
        windowproblem = copy(problem, υ0, τ0, τN)
        windows[m] = parallelsolver(parallelcache, parallelsolution, windowproblem; mode)

        τJ = τN # update τJ for next iteration (next window or restart)
        parallelsolution = NSDETimeParallel.TimeParallelSolution(problem, parallelsolver)
        parallelsolution.lastiterate.U .= U
        parallelsolution.lastiterate.T .= T

        # adaptive check
        if isnothing(adaptive)
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
                @↑ coarsesolver.stepsize = h ← hG

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
                @↑ coarsesolver.stepsize = h ← hG

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
    @↑ coarsesolver.stepsize = h ← hG₀

    return solution
end
