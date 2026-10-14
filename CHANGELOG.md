# Changelog

## 0.3.0 (unreleased)

Requires NSDEBase 0.4 and NSDETimeParallel 0.3; the tests use NSDERungeKutta
0.3.

### Breaking

- `MoWi` is immutable. Adaptive strategies replace solver settings in run-local
  state, so repeated runs start from the configured settings and never leak
  adapted steps or weights.
- `ZoomParameters` requires `WeightedDefect` with a `FixedWeight`, following
  the typed criteria of NSDETimeParallel 0.3.
- The final time is a hard boundary: no window is integrated past it.
  Genuinely shortened windows are recorded in `solution.clipped` and excluded
  from ensemble statistics unless `include_clipped=true`.
- `ensemblesem` uses sample variance (divisor M−1) and needs at least two
  windows.
- Supported Julia: `1.10` and later.

### Added

- `MoWiSolution` carries a run `status` (`:not_started`, `:converged`,
  `:budget_exhausted` or `:nonfinite`) and the per-window `clipped` flags.
- Ensemble statistics accept vector-valued observables (componentwise) and keep
  complex observables complex.
- Strategy factors and restart limits are validated.

### Changed

- Window starts and ends use compensated summation; Leap retries restart from
  the previous accepted start. Roundoff-sized endpoint gaps snap to the final
  time with a window-scaled tolerance instead of producing a spurious window.
- An unrecoverable nonfinite window stops the run with `:nonfinite` and is kept
  for inspection; ensemble statistics are then refused. A finite window that
  exhausts its budget advances with a warning and `:budget_exhausted`.
- Cache reuse accepts a changed RHS of the same type, and chunk schedules
  follow changes of window length.
- Stored Zoom weights are floored at one.
- Fixed, Leap and Zoom windows accept fine solvers without a `stepsize` field.
- Sources, tests and examples follow the repository JuliaFormatter
  configuration (Blue style), checked in CI.

### Fixed

- Reuse resets seeded starts and restart counts, and the seeded start reaches
  the window problem. MPI runs broadcast each collected window, including with
  `PipelinedMPIBackend`.
- Nonfinite windows with empty residual histories restart cleanly, and retries
  are seeded from valid preceding data.
- Scalar convergence-plot labels and duplicate phase-plot legend entries.

## 0.2.1

### Fixed

- Lower the declared Julia minimum from 1.10 to 1.6, matching the supported
  source and README. The full suite passes on Julia 1.6.7 and 1.10.10.
- Standard-library compatibility bounds admit the unversioned stdlibs used in
  Julia 1.6 test environments.

### Tests and maintenance

- Re-enable Aqua's persistent-task check following package registration.
- Use the same narrow RecipesBase.recipetype exemption as the other packages
  when checking recipe methods for type piracy.
- Refresh workflow comments to describe registered dependencies.

## 0.2.0

Requires NSDEBase 0.3.1, NSDERungeKutta 0.2 and NSDETimeParallel 0.2.

### Changed (breaking)

- One moving-window driver replaces the four numbered ones; the strategy is
  selected by the parameter type: `StretchParameters` (AMoWi-1),
  `LeapParameters` (AMoWi-2), `ZoomParameters` (AMoWi-3). The `isstrategy`
  keyword is gone.
- `AdaptiveMoWiParameters` is removed (export dropped). Migrate to one of the
  three typed parameter sets.
- `MoWi(parallelsolver; adaptive, τ, Δτ=τ)`: `Δτ` defaults to `τ`; `Δτ > τ`
  and `Δτ ≤ 0` throw `ArgumentError`; `ZoomParameters` is refused when the
  tolerance carries no weights or has `updatew = true`.
- Window time spans are moved through `NSDETimeParallel.shiftwindow!`.
- A `NewtonFailure` from an implicit coarse or fine solver reaches the caller;
  it is not counted as a window restart.
- Supported Julia: `1.6` and later (was `1.10`). Verified on 1.6–1.13.

### Added

- `MoWiSolution` interpolation `solution(t)`, `windowspan`, `flatten`;
  `ensemblemean`, `ensemblevariance`, `ensemblesem`.
- Strategies page and API page in the docs; Aqua in the test suite; tests
  with an implicit coarse solver.
