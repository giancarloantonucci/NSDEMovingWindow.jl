# Changelog

## 0.3.0 development revision (unreleased)

- Use compensated window starts and endpoints, with rollback on Leap retries and a window-scaled endpoint snap tolerance. Retain a hard final-time boundary and the clipped-window statistics policy.
- Snap roundoff-sized endpoint gaps to avoid spurious windows. Record genuinely clipped windows and exclude them from ensemble statistics by default (include_clipped=true opts in). Support componentwise vector observables.
- Cap all attempted windows at the final time and accept fine-solver adapters without a step-size field when the strategy does not need one.
- Preserve complex ensemble observables and use sample variance for independence-based SEM estimates, explicitly documenting correlation limits. Pin intended warning paths in tests.
- Allow same-type RHS changes during cache reuse. Fix scalar convergence labels and duplicate phase-plot legend entries. Document that the pipeline preserves only the seeded window initial value.
- Reset seeded starts and restart counts on reuse. Pass the seeded start through the window problem and broadcast collected MPI windows before reuse, including the pipelined backend.
- Stop after unrecoverable nonfinite windows with `solution.status == :nonfinite`; retain the failed window for inspection and refuse ensemble statistics for that run. Finite budget-exhausted windows retain explicit status.
- Floor stored Zoom weights at one and validate all strategy factors and restart limits.
- Handle empty residual histories on nonfinite windows and seed retries from valid preceding data rather than invalid output.
- Window adaptation replaces run-local immutable solver settings. Zoom carries a fixed-weight setting explicitly; repeated runs do not leak adapted steps or weights. Chunk schedules follow changes in window length. Updated dependency versions and documentation.
- All Julia sources, tests, and examples follow the repository JuliaFormatter configuration. Historical implementation notes have been removed from current API descriptions.

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
