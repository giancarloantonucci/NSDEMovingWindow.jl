# The moving-window strategies differ in what they adjust on a restart or
# acceptance, so they live in the small `adapt_failure!`/`adapt_success!` hooks
# below, dispatched on the parameter type. `adaptive === nothing` is fixed MoWi:
# never restart, never adapt.

nominal_step(solver) = hasproperty(solver, :stepsize) ? solver.stepsize.h : 0.0

"Mutable window-length/shift state threaded through the adaptive hooks."
mutable struct MoWiState{T<:Real}
    τ::T
    Δτ::T
    # One original chunk length τ₀/N: the shift floor. Flooring at the current
    # τ/N is not enough: Stretch shrinks τ, so that floor would shrink with it
    # and the advance per window could fall to one fine step.
    Δτmin::T
    # max(N·h_F, Δτmin): each chunk holds at least one fine step and the window
    # still contains its shift.
    τmin::T
end

# ------------------------------------------------------------ adaptive hooks --

# Unreachable for failures: fixed MoWi never restarts.
adapt_failure!(::Nothing, state, parallelsolver) = parallelsolver
adapt_success!(::Nothing, state, parallelsolver) = parallelsolver

# AMoWi-1 (Stretch): τ, Δτ and the coarse step scale together (N chunks fixed
# ⇒ h_G ∝ τ); the fine step moves only under `fineupdate` (the analysis-model
# variant; otherwise h_F stays fixed).
function adapt_failure!(p::StretchParameters, state, parallelsolver)
    @↓ δτ⁻, fineupdate = p
    hF = nominal_step(parallelsolver.finesolver)
    # The floors are absolute (set when MoWiState is built): under persistent
    # failure τ must not decay to zero and the advance per window must not decay
    # to one fine step.
    τnew = max(state.τ * δτ⁻, state.τmin)
    if fineupdate
        # The model keeps h_F ∝ τ, so the fine step shrinks only as far as τ
        # does: once τ sits on its floor, a restart leaves h_F alone instead of
        # shrinking it, and the work per window, without bound.
        hF *= state.τ * δτ⁻ >= state.τmin ? δτ⁻ : oftype(hF, min(τnew / state.τ, 1))
        parallelsolver = remake(
            parallelsolver; finesolver=withstepsize(parallelsolver.finesolver, hF)
        )
    end
    state.Δτ = clamp(state.Δτ * δτ⁻, state.Δτmin, τnew)
    state.τ = τnew
    parallelsolver = remake(
        parallelsolver;
        coarsesolver=withstepsize(
            parallelsolver.coarsesolver, max(parallelsolver.coarsesolver.stepsize.h * δτ⁻, hF)
        ),
    )
    return parallelsolver
end

function adapt_success!(p::StretchParameters, state, parallelsolver)
    @↓ δτ⁺, fineupdate = p
    state.Δτ = min(state.Δτ * δτ⁺, state.τ)
    state.τ *= δτ⁺
    if fineupdate
        parallelsolver = remake(
            parallelsolver;
            finesolver=withstepsize(
                parallelsolver.finesolver, parallelsolver.finesolver.stepsize.h * δτ⁺
            ),
        )
    end
    parallelsolver = remake(
        parallelsolver;
        coarsesolver=withstepsize(
            parallelsolver.coarsesolver, parallelsolver.coarsesolver.stepsize.h * δτ⁺
        ),
    )
    return parallelsolver
end

# AMoWi-2 (Leap): only the shift moves; growth is capped by the window length
# (the overlap constraint). Shrinkage is floored at one original chunk length
# τ₀/N (`state.Δτmin`): a shift below one chunk re-solves N−1 chunks of known
# data for under a chunk of progress, and without the floor Δτ decays
# geometrically across consecutively failing windows until the loop cannot make
# progress. The floor is a termination guard.
function adapt_failure!(p::LeapParameters, state, parallelsolver)
    return (state.Δτ=max(state.Δτ * p.δΔτ⁻, state.Δτmin); parallelsolver)
end
function adapt_success!(p::LeapParameters, state, parallelsolver)
    return (state.Δτ=min(state.Δτ * p.δΔτ⁺, state.τ); parallelsolver)
end

# AMoWi-3 (Zoom): multiplicative, per restart, on the base weight of the
# weighted defect criterion: after ι restarts the weight is δ^ι·w. Failure
# increases w; success decreases it, with the stored base floored at one.
function zoom_solver(parallelsolver, factor)
    c = parallelsolver.criterion
    return remake(
        parallelsolver;
        criterion=NSDETimeParallel.WeightedDefect(;
            ϵ=c.ϵ, p=c.p, weight=NSDETimeParallel.FixedWeight(max(1, c.weight.w * factor))
        ),
    )
end
adapt_failure!(p::ZoomParameters, state, parallelsolver) = zoom_solver(parallelsolver, p.δw⁻)
adapt_success!(p::ZoomParameters, state, parallelsolver) = zoom_solver(parallelsolver, p.δw⁺)

# Carry the rounding remainder of accepted shifts. Rejected Leap shifts must
# be recomputed from the previous accepted start and its saved remainder.
function window_shift(start, shift, remainder)
    v = shift + remainder
    next = start + v
    error = abs(start) >= abs(v) ? (start - next) + v : (v - next) + start
    return next, error
end

# Stretch/Zoom retry the same start; Leap retries the last accepted shift.
function restart_start(
    ::AbstractAdaptiveMoWiParameters, start, previous, shift, error, previous_error
)
    return (start, error)
end
function restart_start(::LeapParameters, start, previous, shift, error, previous_error)
    return window_shift(previous, shift, previous_error)
end

# Leap leaves τ alone, so if the first window cannot converge no shift policy
# can save it: window 1 is exempt from Leap restarts.
allowsrestart(adaptive, m) = !(adaptive isa Nothing) && !(adaptive isa LeapParameters && m == 1)

# ------------------------------------------------------------- the skeleton --

"""
    mowi!(cache, solution, problem, mowi; verbose=false, kwargs...)

runs the moving-window loop: re-target the shared window problem and parallel
cache, seed the chunk starts from the previous window over the overlap (coarse
elsewhere, through `makeGs`), solve the window with the time-parallel
subroutine, then either accept and shift or restart under the adaptive policy.
MPI runs collect and broadcast each window before the next seed is formed.
The pipeline preserves the seeded window initial value but replaces interior
overlap guesses during its coarse start-up; bulk backends use those guesses.
`nocollect=true` is refused: the driver needs full window trajectories. Extra
`kwargs` (e.g. `history`) go straight to the subroutine.
"""
function mowi!(
    cache::MoWiCache,
    solution::MoWiSolution,
    problem::AbstractInitialValueProblem,
    mowi::MoWi;
    verbose::Bool=false,
    backend::AbstractPararealBackend=SerialBackend(),
    nocollect::Bool=false,
    kwargs...,
)
    nocollect && throw(ArgumentError("MoWi requires collected window trajectories"))
    @↓ windowproblem, windowcache = cache
    typeof(windowproblem.rhs) === typeof(problem.rhs) || throw(
        ArgumentError("MoWi cache RHS type differs; build a new cache for a different RHS type")
    )
    windowproblem.rhs = problem.rhs
    @↓ parallelsolver, adaptive = mowi
    axes(problem.u0) == axes(windowcache.U[1]) && eltype(problem.u0) == eltype(windowcache.U[1]) ||
        throw(
            DimensionMismatch("MoWi cache state shape or element type differs; build a new cache")
        )
    parallelsolver.parameters.N == length(windowcache.U) ||
        throw(DimensionMismatch("MoWi cache chunk count differs; build a new cache"))
    W = windowtime(problem, mowi)
    eltype(windowproblem.tspan) === W ||
        throw(ArgumentError("MoWi cache time type differs; build a new cache"))
    parallelsolver = deepcopy(parallelsolver)
    @↓ windows, restarts, clipped = solution
    fill!(restarts, 0)
    fill!(clipped, false)
    solution.status = :converged
    @↓ (t0, tN) ← tspan = problem
    N = parallelsolver.parameters.N

    Δτmin = W(mowi.τ) / N
    τmin = if adaptive isa StretchParameters
        max(N * nominal_step(parallelsolver.finesolver), Δτmin)
    else
        Δτmin
    end
    # Under `fineupdate` the model variant shrinks h_F too, so this is conservative.
    state = MoWiState{W}(mowi.τ, mowi.Δτ, Δτmin, τmin)
    τ0 = τ0_accepted = W(t0)
    start_error = accepted_start_error = zero(τ0)
    τJ = W(t0) # end of the covered region so far
    time_roundoff = 64eps(max(abs(float(t0)), abs(float(tN))))
    τN = τ0 + state.τ
    isrestart = false

    M = M₀ = length(windows)
    m = 1
    while true
        # Snap tolerance-sized shortfalls as well as clipping genuine overshoots.
        # Use the same compensated addition for the end and the next start:
        # non-overlapping windows then meet without a rounding-sized gap.
        nominal_end, _ = window_shift(τ0, state.τ, start_error)
        endtol = max(time_roundoff, 1e-9 * state.τ)
        τN = nominal_end >= tN - endtol ? W(tN) : nominal_end
        duration = (τN == nominal_end ? state.τ : τN - τ0) / N
        if m > M # grow storage under adaptive shrinking
            append!(windows, Vector{AbstractTimeParallelSolution}(undef, M₀))
            append!(restarts, zeros(Int, M₀))
            append!(clipped, falses(M₀))
            M += M₀
        end
        clipped[m] = nominal_end > tN + endtol
        verbose && @info "MoWi window" m τ0 τJ τN restarts = restarts[m]

        # Re-target the shared window problem and cache:
        tspan = (τ0, τN)
        @↑ windowproblem = tspan
        NSDETimeParallel.shiftwindow!(windowcache, τ0, τN; duration=duration)

        # Seed the chunk starts: previous data over the overlap (the failed
        # attempt on a same-start restart), the coarse solver elsewhere.
        begin
            @↓ makeGs, U, T = windowcache
            usable_retry =
                isrestart &&
                !(adaptive isa LeapParameters) &&
                windows[m].diagnostics.status != :nonfinite
            source = usable_retry ? windows[m] : (m > 1 ? windows[m - 1] : nothing)
            for n in 1:N
                if source === nothing && n == 1
                    copyto!(U[n], problem.u0)
                    makeGs[n] = false
                elseif source !== nothing && T[n] ≤ τJ
                    # copyto!, never rebind: chunk problems reference these buffers.
                    copyto!(U[n], source(T[n]))
                    makeGs[n] = false
                else
                    makeGs[n] = true
                end
            end
        end

        windowproblem.u0 = copy(windowcache.U[1])
        windows[m] = parallelsolver(windowcache, windowproblem; backend, kwargs...)
        if backend isa Union{MPIBackend,PipelinedMPIBackend}
            # Every rank needs the same accepted/retry trajectory to seed the next window.
            mpi = NSDETimeParallel.MPI
            windows[m] = mpi.bcast(
                mpi.Comm_rank(mpi.COMM_WORLD) == 0 ? windows[m] : nothing, 0, mpi.COMM_WORLD
            )
        end

        err = isempty(windows[m].errors) ? Inf : windows[m].errors[end]
        failed = !(windows[m].diagnostics.status in (:converged, :exact))
        if failed && allowsrestart(adaptive, m) && restarts[m] < adaptive.R
            restarts[m] += 1
            parallelsolver = adapt_failure!(adaptive, state, parallelsolver)
            τ0, start_error = restart_start(
                adaptive, τ0, τ0_accepted, state.Δτ, start_error, accepted_start_error
            )
            τN = τ0 + state.τ
            isrestart = true
        else
            if windows[m].diagnostics.status == :nonfinite
                solution.status = :nonfinite
                @warn "MoWi stopped at a nonfinite window" m restarts = restarts[m]
                m += 1 # retain the failed window for diagnostics; do not seed another
                break
            end
            if failed
                solution.status = :budget_exhausted
                # Reaching the restart cap returns an unconverged window.
                # Report this even when verbose progress logging is disabled.
                @warn "MoWi window accepted WITHOUT convergence — convergence is not guaranteed past this window" m error =
                    err restarts = restarts[m]
            end
            # Success hooks fire on success only. Firing them on a forced
            # acceptance would invert the strategies' meaning: Stretch would
            # regrow τ straight after R consecutive shrinks failed to converge
            # the window, and Zoom would tighten w in response to a failure.
            # Leaving the state alone keeps the next window at the conservative
            # post-failure settings.
            if adaptive !== nothing && !failed
                parallelsolver = adapt_success!(adaptive, state, parallelsolver)
            end
            τ0_accepted = τ0
            accepted_start_error = start_error
            τJ = τN
            m += 1
            isrestart = false
            if τJ ≥ tN
                break
            end
            τ0, start_error = window_shift(τ0_accepted, state.Δτ, accepted_start_error)
            τN = τ0 + state.τ
        end
    end

    resize!(windows, m - 1)
    resize!(restarts, m - 1)
    resize!(clipped, m - 1)
    return solution
end
