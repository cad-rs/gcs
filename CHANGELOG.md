# Changelog

All notable changes to `rcad-gcs` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the versioning
policy is [Semantic Versioning](https://semver.org/) (0.x: breaking changes
are allowed with an entry here).

## 0.2.0 — unreleased

Baseline fixes: phase P0 of [docs/CATCHUP_PLAN.md](docs/CATCHUP_PLAN.md).
Several entries under "Changed" alter behaviour or signatures, hence the
minor version.

### Changed

- **A solve that does not converge now puts back the parameters it started
  from** (2D and 3D-geometric families), per connected component: a
  component that converged keeps its solution, a failed one is left
  untouched. 0.1.0 left the last iterate in place.
  `SolveOptions { on_failure: FailurePolicy::KeepLast, .. }` restores the
  old behaviour for a call. `SolveResult::residual` is the residual of the
  state left behind.
- **`Sketch::dof()` and `SpaceSketch::dof()` are rank-based** — the same
  number `diagnose().dof` reports — and so is the `dof` field of
  `DofReport` / `SpaceDofReport`. Redundant equations no longer drive the
  value negative; over-constraint is read from `diagnose()`. The former
  count (free parameters − equations) is `structural_dof_estimate()` and
  the new `structural_dof` report field.
- **A free cylinder reports 5 degrees of freedom, not 7.** The 3D diagnosis
  discounts the two representation freedoms of the cylinder storage (axis
  point sliding along the axis, axis direction rescaling) unless a
  constraint pins them.
- **A step that fails the line search is no longer taken.** 0.1.0 took the
  last halving of a step even when all eight halvings had failed the
  acceptance test. Now the step is recomputed with a stronger
  regularisation (`LAMBDA_ESCALATION`); if that fails too, the iterate is
  put back and the solve stops as `Stalled`. A numerically singular step is
  retried the same way and ends the solve as `Singular` instead of idling
  to the iteration limit. The acceptance test itself is unchanged (the
  residual may grow by a relative 1e-4).
- **Unsatisfied equations with no free parameter are no longer "converged".**
  A component whose parameters are all fixed used to report
  `converged = true` whatever its residual. It now reports
  `converged = false` and `Stalled` when its equations do not hold.
- `solver::solve`, `solver::solve_components`, `space3d::solver::solve_space`,
  `space3d::solver::solve_components` and `numerics::newton_solve` take the
  parameters as `&mut [f64]` instead of `&mut Vec<f64>`.
- `SolveResult` has a new public field `status`, and `DofReport` /
  `SpaceDofReport` have new fields `rank` and `structural_dof` (code that
  builds these structs by literal needs the new fields).
- `numerics::newton_solve_with_jacobian` is crate-private. It took a
  crate-private argument type and could not be called from outside.

### Added

- `SolveStatus` (`Converged`, `MaxIterations`, `Stalled`, `Singular`,
  `Diverged`), reported in `SolveResult::status`.
- `SolveOptions` (tolerance, iteration limit, failure policy) and
  `FailurePolicy` for the 2D and 3D-geometric families, with
  `Sketch::solve_with`, `Sketch::solve_selected_with`,
  `SpaceSketch::solve_with`, `SpaceSketch::solve_selected_with`,
  `SpaceSketch::solve_incremental_with`, `solver::solve_with`,
  `solver::solve_components_with`, `space3d::solver::solve_space_with`,
  `space3d::solver::solve_components_with` and `numerics::newton_solve_with`.
- `Sketch::structural_dof_estimate` and `SpaceSketch::structural_dof_estimate`.
- A criterion benchmark harness (`cargo bench`, `benches/solve.rs`) and the
  recorded numbers in `benches/BASELINE.md`.
- `.gitattributes` pinning LF line endings.
- Tests for the failure policy (`tests/failure_policy.rs`) and for the DOF
  account (`tests/dof_account.rs`).

### Fixed

- The corpus manifest tests failed on checkouts with CRLF line endings
  (Git for Windows' default); they now compare modulo line endings, and the
  repository pins LF.
- `cargo fmt --check` and `cargo clippy --all-targets --all-features
  -- -D warnings`, the first two CI steps, did not pass on the 0.1.0 tree.
- Gaussian elimination panicked on a NaN matrix entry; non-finite values now
  end the solve as `Diverged`.
- The corpus tests never wrote their JSON reports (`CARGO_TARGET_TMPDIR` is
  only set at compile time).
- `cargo doc` built with 59 unresolved-link warnings (bracketed parameter
  offsets in the entity tables read as links); it is clean now.
- The unit test `large_system_makes_progress` built its grid with
  overlapping parameters, which is why it never converged. With the layout
  fixed it converges and is now `grid_system_converges`.

### Performance

- The normal equations are accumulated over the non-zero Jacobian entries
  only (bitwise the same result).
- The rigid-pose solver and diagnosis make unknowns only of the bodies a
  system's constraints reference, and the finite-difference rows no longer
  copy every pose per column.
- The diagnosis builds and factors the Gram matrix of the independent rows
  once instead of once per dependent row. With many redundant constraints
  that is the difference between seconds and minutes (a 780-parameter
  system with 255 redundant rows: about 113 s → 2.4 s).
- Assembly corpus gate on the reference machine: 21.0 s → 1.7 s optimised,
  411 s → 69 s unoptimised, with identical metrics. See
  `benches/BASELINE.md` for per-size numbers.

## 0.1.0 — 2026-11-13

First standalone release, extracted from the `rcad` monorepo
(`libs/rcad-gcs`, itself MIT OR Apache-2.0).

### Added

- **2D sketch family** — `Sketch` with points, lines, circles and arcs, and
  the constraint set: Fixed, Coincident, PointDistance, LineLength,
  Horizontal, Vertical, Parallel, Perpendicular, EqualLength, Angle,
  EqualRadius, PointOnCircle, PointOnLine, Tangent, CircleCircleTangent,
  Radius, ArcArcTangent, Symmetric, Concentric, Midpoint, Diameter.
- **3D geometric family** — `space3d::SpaceSketch` with points, lines,
  planes, spheres and cylinders, and their incidence / distance / length /
  angle / radius / tangency / symmetry constraints.
- **3D rigid-pose family** — `rigid3d` with per-body 6-DoF
  (translation + quaternion manifold), the Fixed / Coincident / Distance /
  Parallel / Angle / Perpendicular constraint set, analytic point-type
  Jacobian columns, localized per-component solves and a warm-start entry
  point.
- **Shared Newton-Raphson least-squares core** — Tikhonov-regularised normal
  equations with a damped line search; under-constrained systems get the
  minimum-norm step.
- **Rank-revealing diagnosis** — redundant / conflicting /
  partially-redundant classification with a true rank-based DOF count.
- **`serde` feature** — `Serialize`/`Deserialize` on the constraint data
  model, off by default (the core depends only on `glam`).
- **Test corpora** — 2D convergence, 3D geometric, and assembly (rigid-pose)
  corpora with metric gates.
