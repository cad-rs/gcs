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
- **The rigid-pose solve runs on the shared Newton core.** It gains the
  line search and the regularisation escalation of the other families.
  `ConstraintSolver::solve` and `solve_from` now solve every connected
  cluster on its own (as `solve_components` always did) instead of the
  whole system at once: `success` means every cluster's own residual norm
  is below the tolerance, and `iterations` is the sum over clusters.
- **A body pinned in full by a `Fixed` constraint is grounded.** When the
  constraint gives both position and orientation, the body is placed at
  exactly that pose and is no longer solved for; a warm start that has it
  elsewhere is corrected, not iterated. Parts mated only to the ground no
  longer couple with one another. Results change slightly for the better:
  on the assembly corpus the worst round-trip drift went from 1.2e-7 to
  1.0e-13, and the nine configurations that flipped branch after a large
  perturbation no longer do. A translation-only `Fixed` is an ordinary
  equation as before.
- `solver::solve`, `solver::solve_components`, `space3d::solver::solve_space`,
  `space3d::solver::solve_components` and `numerics::newton_solve` take the
  parameters as `&mut [f64]` instead of `&mut Vec<f64>`.
- `Sketch::add_constraint` returns the index of the new constraint instead
  of nothing. Code that used the call as the last expression of a `()`
  block or closure needs a semicolon.
- `Sketch::components`, `diagnose`, `mobility` and the DOF accounts ignore
  suppressed and removed constraints.
- `SolveResult` has a new public field `status`, and `DofReport` /
  `SpaceDofReport` have new fields `rank` and `structural_dof` (code that
  builds these structs by literal needs the new fields).
- `numerics::newton_solve_with_jacobian` is crate-private. It took a
  crate-private argument type and could not be called from outside.
- `space3d::diagnosis::SpaceDiagnosis` is a type alias of `Diagnosis` (the
  fields were already the same; serialized key order now follows
  `Diagnosis`), and the 3D threshold constants are re-exports of the 2D
  ones.

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
- **Editing a 2D sketch in place.** Entity ids and constraint indices are
  stable (a removed slot is left behind, nothing shifts):
  - `Sketch::remove_constraint`, `Sketch::set_constraint_enabled` (a
    suppressed dimension is a reference dimension), `Sketch::replace_constraint`,
    `Sketch::constraint_state` and the `ConstraintState` enum;
  - `Sketch::remove_entity`, which removes the constraints on the entity
    with it, and `Sketch::is_entity_removed`;
  - `Sketch::set_constraint_value` / `constraint_value` to drive a dimension
    without rebuilding, and `Sketch::measure` to read one;
  - `Sketch::topology_revision`, which changes only on structural edits;
  - `Sketch::validate` and `SketchError` for dangling references, wrong
    entity kinds and non-finite values;
  - on `Constraint`: `references`, `value`, `set_value`, `measure`; on
    `PointRef`: `entity`, `expected_kinds`.
- `Sketch::structural_dof_estimate` and `SpaceSketch::structural_dof_estimate`.
- `rigid3d::AssemblyDiagnosis::rank` and `dependency_groups` (by mate id):
  the assembly diagnosis now reports what the 2D and 3D ones do.
- A criterion benchmark harness (`cargo bench`, `benches/solve.rs`) and the
  recorded numbers in `benches/BASELINE.md`.
- `.gitattributes` pinning LF line endings.
- Tests for the failure policy (`tests/failure_policy.rs`) and for the DOF
  account (`tests/dof_account.rs`).

### Fixed

- **Two spline constraints of the 3D family had wrong derivatives**, and no
  test covered any of the four.
  - `PointOnSpline` held the nearest point on the curve fixed while
    differentiating, which is right only when the point is already on the
    curve. Debug builds compare the analytic Jacobian with finite
    differences at the start of every solve and panic on a mismatch, so a
    point that started off the curve panicked there; release builds
    iterated with the inexact rows. The rows now follow the nearest point,
    and handle a point that lies beyond an end of the curve.
  - `SplineControlPointSymmetric` had no derivative with respect to its
    mirror line, and incomplete ones for the control points. A free mirror
    line never moved and was not connected to the spline in the
    decomposition. The constraint now shares the rows of `Symmetric`, which
    has the same equations.
  - `PointOnSpline` on a spline without control points underflowed an
    index while building its rows (a panic in debug builds).

  All four spline constraints are now in the per-variant
  Jacobian-versus-finite-difference test and have solve tests.
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

### Internal

- The 2D and 3D-geometric families run through one implementation of the
  decomposition, the per-component solve and the diagnosis
  (`src/system.rs`); each family is now an adapter that supplies residuals
  and Jacobian rows, plus the plane and cylinder gauge for 3D. Results are
  unchanged: the corpus reports are identical to the last digit.
- The rigid-pose family is on that pipeline too: its diagnosis is the
  shared diagnosis run over six columns per body, and its solve is the
  shared per-component solve. The Newton core steps through a `StepSpace`,
  so a pose (translation plus unit quaternion) is stepped by the same loop
  as plain parameters. The family's own Gauss–Newton loop and its copy of
  the QR analysis are gone.

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
- Grounded bodies decouple what hangs off them: a hub with 40 parts pinned
  to it solves about 6× faster than before the change. Small assemblies
  pay a fixed cost for the shared pipeline (a 5-body chain solves in
  109 µs instead of 66 µs), and the assembly diagnosis is 1.4× to 2.3×
  slower while reporting more.
- Assembly corpus gate on the reference machine: 21.0 s → 1.5 s optimised,
  411 s → 38 s unoptimised. See `benches/BASELINE.md` for per-size numbers.

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
