# src/mowi/mowi.jl
#
# THE moving-window loop — written once. The adaptive strategies (thesis,
# Ch. 5) differ only in what they adjust on a restart or an acceptance, so
# they live in the small `adapt_failure!`/`adapt_success!` hooks below,
# dispatched on the parameter type. `adaptive === nothing` is fixed MoWi:
# never restart, never adapt.

"Mutable window-length/shift state threaded through the adaptive hooks."
mutable struct MoWiState
    τ::Float64
    Δτ::Float64
    Δτmin::Float64 # one ORIGINAL chunk length τ₀/N: the shift floor. Flooring
                   # at the CURRENT τ/N is not enough — Stretch shrinks τ, so
                   # that floor shrinks with it and the advance per window can
                   # fall to one fine step (~1e5 windows to cover a domain).
    τmin::Float64  # max(N·h_F, Δτmin): each chunk holds ≥ 1 fine step and the
                   # window still contains its shift.
end

# ------------------------------------------------------------ adaptive hooks --

adapt_failure!(::Nothing, state, parallelsolver) = nothing # unreachable: fixed MoWi never restarts
adapt_success!(::Nothing, state, parallelsolver) = nothing

# AMoWi-1 (Stretch): τ, Δτ and the coarse step scale together (N chunks fixed
# ⇒ h_G ∝ τ); the fine step moves only under `fineupdate` (the analysis-model
# variant — thesis fixes h_F in the practical scheme).
function adapt_failure!(p::StretchParameters, state, parallelsolver)
    @↓ δτ⁻, fineupdate = p
    hF = parallelsolver.finesolver.stepsize.h
    if fineupdate
        hF *= δτ⁻
        parallelsolver.finesolver.stepsize.h = hF
    end
    # Floors are ABSOLUTE (set at construction in MoWiState): under persistent
    # failure τ must not decay to zero and the advance per window must not
    # decay to one fine step.
    τnew = max(state.τ * δτ⁻, state.τmin)
    state.Δτ = clamp(state.Δτ * δτ⁻, state.Δτmin, τnew)
    state.τ = τnew
    parallelsolver.coarsesolver.stepsize.h = max(parallelsolver.coarsesolver.stepsize.h * δτ⁻, hF)
    return nothing
end

function adapt_success!(p::StretchParameters, state, parallelsolver)
    @↓ δτ⁺, fineupdate = p
    state.Δτ = min(state.Δτ * δτ⁺, state.τ)
    state.τ *= δτ⁺
    if fineupdate
        parallelsolver.finesolver.stepsize.h *= δτ⁺
    end
    parallelsolver.coarsesolver.stepsize.h *= δτ⁺
    return nothing
end

# AMoWi-2 (Leap): only the shift moves; growth is capped by the window length
# (the overlap constraint — thesis §5.2). Shrinkage is FLOORED at one ORIGINAL
# chunk length τ₀/N (`state.Δτmin`): a shift below one chunk re-solves N−1
# chunks of known data for under a chunk of progress, and without the floor
# Δτ decays geometrically across consecutively failing windows — the advance
# per window tends to zero (hitting exactly 0.0 after ~10³ halvings), so the
# loop never terminates. The thesis is silent on the floor; it is a
# termination guard, not policy.
adapt_failure!(p::LeapParameters, state, parallelsolver) =
    (state.Δτ = max(state.Δτ * p.δΔτ⁻, state.Δτmin); nothing)
adapt_success!(p::LeapParameters, state, parallelsolver) = (state.Δτ = min(state.Δτ * p.δΔτ⁺, state.τ); nothing)

# AMoWi-3 (Zoom): multiplicative, per restart, on the base weight of the
# weighted proximity function (thesis Thm. 5.3: after ι restarts the weight is
# δ^ι·w). Failure loosens (w grows), success tightens (w shrinks; ψ₂ clamps
# w ≥ 1).
function adapt_failure!(p::ZoomParameters, state, parallelsolver)
    weights = parallelsolver.tolerance.weights
    weights.w = weights.w * p.δw⁻
    return nothing
end

function adapt_success!(p::ZoomParameters, state, parallelsolver)
    weights = parallelsolver.tolerance.weights
    weights.w = weights.w * p.δw⁺
    return nothing
end

# On a restart, Stretch and Zoom re-attempt the SAME window start; Leap
# re-attempts the SHIFT at its reduced size (thesis: Δτ_{m,r} = δ⁻ Δτ_{m,r-1}).
restart_start(::AbstractAdaptiveMoWiParameters, τ0, τ0_accepted, Δτ) = τ0
restart_start(::LeapParameters, τ0, τ0_accepted, Δτ) = τ0_accepted + Δτ

# Leap leaves τ alone, so if the FIRST window cannot converge no shift policy
# can save it (thesis §5.2): window 1 is exempt from Leap restarts.
allowsrestart(adaptive, m) = !(adaptive isa Nothing) && !(adaptive isa LeapParameters && m == 1)

# ------------------------------------------------------------- the skeleton --

"""
    mowi!(cache, solution, problem, mowi; verbose=false, kwargs...)

runs the moving-window loop: re-target the shared window problem and parallel
cache, seed the chunk starts from the previous window over the overlap (coarse
elsewhere, through `makeGs`), solve the window with the time-parallel
subroutine, then either accept and shift or restart under the adaptive policy.
Extra `kwargs` (e.g. `mode`, `backend`, `saveiterates`) go straight to the
subroutine.
"""
function mowi!(cache::MoWiCache, solution::MoWiSolution, problem::AbstractInitialValueProblem, mowi::MoWi; verbose::Bool=false, kwargs...)
    @↓ windowproblem, windowcache = cache
    @↓ parallelsolver, adaptive = mowi
    @↓ ϵ = parallelsolver.tolerance
    @↓ windows, restarts = solution
    @↓ (t0, tN) ← tspan = problem
    N = parallelsolver.parameters.N

    Δτmin = float(mowi.τ) / N
    τmin = max(N * parallelsolver.finesolver.stepsize.h, Δτmin) # static: with
    # `fineupdate` the model variant shrinks h_F too, making this conservative
    state = MoWiState(float(mowi.τ), float(mowi.Δτ), Δτmin, τmin)
    τ0 = τ0_accepted = float(t0)
    τJ = float(t0) # end of the covered region so far
    τN = τ0 + state.τ
    isrestart = false

    M = M₀ = length(windows)
    m = 1
    while true
        if m > M # grow storage under adaptive shrinking
            append!(windows, Vector{AbstractTimeParallelSolution}(undef, M₀))
            append!(restarts, zeros(Int, M₀))
            M += M₀
        end
        verbose && @info "MoWi window" m τ0 τJ τN restarts = restarts[m]

        # Re-target the shared window problem and cache:
        tspan = (τ0, τN)
        @↑ windowproblem = tspan
        NSDETimeParallel.shiftwindow!(windowcache, τ0, τN)

        # Seed the chunk starts: previous data over the overlap (the failed
        # attempt on a same-start restart), the coarse solver elsewhere.
        if m > 1 || isrestart
            @↓ makeGs, U, T = windowcache
            source = (isrestart && !(adaptive isa LeapParameters)) ? windows[m] : windows[m-1]
            for n = 1:N
                if T[n] ≤ τJ
                    copyto!(U[n], source(T[n])) # copyto!, NEVER rebind: chunk problems reference these buffers
                    makeGs[n] = false
                else
                    makeGs[n] = true
                end
            end
        end

        windows[m] = parallelsolver(windowcache, windowproblem; kwargs...)

        err = windows[m].errors[end]
        failed = err > ϵ || isnan(err)
        if failed && allowsrestart(adaptive, m) && restarts[m] < adaptive.R
            restarts[m] += 1
            adapt_failure!(adaptive, state, parallelsolver)
            τ0 = restart_start(adaptive, τ0, τ0_accepted, state.Δτ)
            τN = τ0 + state.τ
            isrestart = true
        else
            if failed
                # This branch is OUTSIDE the thesis: §5.2 prescribes restarting
                # until convergence, and Thms 5.1–5.3 exist to guarantee finitely
                # many restarts suffice. The R cap is an implementation guard, and
                # accepting an unconverged window quietly voids every Ch. 5
                # guarantee downstream of it — so the warning is UNCONDITIONAL,
                # never verbose-gated.
                @warn "MoWi window accepted WITHOUT convergence — Ch. 5 guarantees do not cover the trajectory past this window" m error = err restarts = restarts[m]
            end
            # Success hooks fire on SUCCESS only. Firing them on a forced
            # acceptance would invert the strategies' semantics in a branch the
            # thesis never defines: Stretch would re-grow τ immediately after R
            # consecutive shrinks failed to converge the window, and Zoom would
            # TIGHTEN w in response to a failure. Leaving the state alone keeps
            # the next window at the conservative post-failure parameters.
            (adaptive isa Nothing || failed) || adapt_success!(adaptive, state, parallelsolver)
            τ0_accepted = τ0
            τJ = τN
            m += 1
            isrestart = false
            if τJ ≥ tN
                break
            end
            τ0 = τ0_accepted + state.Δτ
            τN = τ0 + state.τ
        end
    end

    resize!(windows, m - 1)
    resize!(restarts, m - 1)
    return solution
end
