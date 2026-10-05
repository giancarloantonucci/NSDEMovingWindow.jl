# Changelog

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


