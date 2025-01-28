function mowi_3(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem, mowi::MoWi;
    isstrategy=0, kwargs...)

    @↓ parallelsolver, adaptive, τ, Δτ = mowi
    @↓ finesolver, coarsesolver, parameters, tolerance = parallelsolver
    @↓ N = parameters
    @↓ ϵ = tolerance
    if adaptive !== nothing
        @↓ δτ⁻, δτ⁺, δΔτ⁻, δΔτ⁺, R, fineupdate = adaptive
    end
    @↓ τ0, τN, τJ = cache
    @↓ windows, restarts = solution
    @↓ (t0, tN) ← tspan = problem

    # main loop
    needsrestart = false
    m = 1
    while τN < tN

        if needsrestart && isstrategy == 2
            Δτ *= δΔτ⁻
        end

        # function update_window_parameters(τ0, τJ, τN, τ, Δτ, isstrategy, needsrestart, adaptive, n_restarts)
        #     if needsrestart
        #         if isstrategy == 1
        #             if adaptive !== nothing
        #                 @↓ δτ⁻, δτ⁺, R, fineupdate = adaptive
        #             end

        #             τ *= δτ⁻
        #             Δτ *= δτ⁻

        #             τN = τ0 + τ
        #         end
        #     else
        #         τ0 += Δτ
        #         τJ = τN
        #         τN = τ0 + τ
        #     end
        #     return τ0, τJ, τN
        # end

        # τ0, τJ, τN = update_window_parameters(τ0, τJ, τN, τ, Δτ, isstrategy, needsrestart, adaptive, restarts[m])

        if m > 1
            τ0 += Δτ
            τJ = τN
            τN = τ0 + τ
        end

        println("m = $m, τ0 = $τ0, τJ = $τJ, τN = $τN")

        windowproblem = copy(problem, τ0, τN) # u0 here is the wrong one; careful

        if m > 1
            cache = NSDETimeParallel.TimeParallelCache(windowproblem, parallelsolver)
            @↓ makeGs, Uₘ ← U, Tₘ ← T = cache

            for n = 1:N

                println("n = $n, Tₘ = $Tₘ, Tₘ[n] = $(Tₘ[n])")

                # between τ0 and τJ copy parallelsolution
                if τ0 ≤ Tₘ[n] ≤ τJ
                    if needsrestart && isstrategy == 1
                        Uₘ[n] = windows[m](Tₘ[n])
                        makeGs[n] = false
                    else
                        Uₘ[n] = windows[m-1](Tₘ[n])
                        makeGs[n] = false
                    end
                end

                if isstrategy == 2
                    if τJ < Tₘ[n] ≤ τN
                        Uₘ[n] = windows[m](Tₘ[n])
                        # makeGs[n] = true or false?
                    end
                else
                    # between τJ and τN coarse guess
                    if τJ < Tₘ[n] ≤ τN
                        # windowchunkproblem = copy(problem, Uₘ[n], Tₘ[n], Tₘ[n+1])
                        # windowchunksolution = coarsesolver(windowchunkproblem)
                        # Uₘ[n] = windowchunksolution(Tₘ[n])
                        makeGs[n] = true
                    end
                end

            end

            # solve with time-parallel algorithm
            if isstrategy == 3
                # parallelsolver.cache.a = ...
            end

        end

        windows[m] = solve(windowproblem, parallelsolver; kwargs...)

        if windows[m].errors[end] > ϵ
            needsrestart = true
            restarts[m] += 1
        else
            needsrestart = false
            if isstrategy == 1
                τ *= δτ⁺
            end
            # move on
            m += 1
        end

        #   update windowproblem / cache (that is, τ0 and τN)
    end

end
