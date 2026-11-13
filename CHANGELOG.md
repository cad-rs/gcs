# Changelog

All notable changes to `rcad-gcs` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the versioning
policy is [Semantic Versioning](https://semver.org/) (0.x: breaking changes
are allowed with an entry here).

## 0.2.0 — unreleased

Phases P0 (baseline fixes) and P1 (one pipeline for the three families, an
editable model) of [docs/CATCHUP_PLAN.md](docs/CATCHUP_PLAN.md). Several
entries under "Changed" alter behaviour or signatures, hence the minor
version; [docs/MIGRATION_0.2.md](docs/MIGRATION_0.2.md) (in Chinese) goes
through them from a host's side.

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
- **Step control is a trust region (Powell's dogleg) instead of a halving
  line search.** 0.1.0 halved the Gauss–Newton step up to eight times,
  accepted a step that let the residual grow by a relative 1e-4, and took
  the last halving even when all eight had failed. Now the Gauss–Newton
  step is taken whole while it fits inside a trust radius that starts out
  unbounded, so an easy solve is unchanged; a step is kept only if it
  lowers the residual, and a rejected one shrinks the radius and bends the
  next trial toward steepest descent, at the cost of a residual evaluation
  and no new factorisation. On the 2D and assembly corpora the work is the
  same to within 0.2% (linear solves and residual evaluations counted) and
  every reported metric is identical, the large-drift branch flips
  included; the 3D corpus passes unchanged. A serial
  3D truss, where errors amplify along the chain, now converges from a
  0.05 perturbation at every length tried; it did not at 160 and 260
  points (`tests/hard_cases.rs`). The constants
  `numerics::LINE_SEARCH_RELAX` and `LINE_SEARCH_STEPS` are gone with the
  line search.
- **A model far from the origin solves.** The tolerance is an absolute
  number (1e-10 by default), and a coordinate of a billion is only stored
  to about 1e-7: no equation on it can be met to 1e-10. Such a solve used
  to end as `Stalled` and put the geometry back. Now each equation is
  weighed by what its own numbers resolve: it is held to the tolerance or
  to four rounding units of the terms it is computed from
  (`numerics::PRECISION_ULPS`), whichever is larger, and the least-squares
  step works in those units, so the rounding of an equation on large
  coordinates is not traded against one that is exact. A solve that ends
  this way reports the new status `SolveStatus::ConvergedToPrecision`,
  counts as converged, and returns a `residual` that can be above the
  tolerance. Near the origin every weight is one and nothing changes: the
  2D corpus is identical. The assembly corpus, which asks for 1e-12 at
  coordinates of some tens, is within reach of the allowance: its worst
  round-trip drift went from 1.0e-13 to 1.6e-13. A contradiction is still
  a failed solve wherever it is larger than that rounding
  (`tests/scale.rs`).
- **The 2D angular relations and equal length are lengths, at any size.**
  `Parallel`, `Perpendicular` and `Angle` were the cross product, the dot
  product and `sin·dot − cos·cross` of the two line vectors, and
  `EqualLength` the difference of the squared lengths: all products of two
  lengths, while every other 2D equation is a length. On small geometry
  they fell under the tolerance before they held (a triangle 1e-5 across
  came back `Converged` with a corner 0.12% of its size out of place; at
  1e-6 with a tolerance of 1e-16 it never converged), and on large geometry
  they dwarfed the rest. Now the three angular relations are stated on the
  unit directions (`sin θ`, `cos θ`, `sin(angle − θ)`) times the geometric
  mean of the two lengths, and equal length is `|d1| − |d2|`. The solution
  sets are the same (an `Angle` is still also met at `angle + π`); the
  residual values, and so a solve's path, are not. The 3D
  `EqualLineLength` gets the same change (at 1e-8 across it used to report
  converged with the line 10% short); the other 3D angular relations were
  already normalised. The same triangle now solves to within 2e-15 of its
  size from 1e-6 to 1e6 across (`tests/scale.rs`). Every metric of the 2D
  and assembly corpora is unchanged; the 3D corpus passes.
- **A 2D component smaller than one unit is held to the tolerance in
  proportion to its size.** The tolerance is an absolute number, and a
  component with freedom left gets no convergence polish, so a small one
  used to stop as soon as its equations were within 1e-10 model units: a
  triangle 1e-4 across, free to turn, was left with its lengths 2.7e-8 of
  its size out, against 7e-16 at unit size. Now the residual of a 2D
  component is measured in units of its size — the larger side of the box
  around the entities its constraints refer to — when that is below one,
  and the same triangle comes back to within 1.4e-15 of its size, in three
  iterations, at every size from 1e-8 to 1e3 (`tests/scale.rs`). A
  component one unit across or larger keeps the absolute tolerance, and
  the four rounding units of its numbers still bound it (so the 2D corpus
  is unchanged, metrics included). The reported `residual` stays in model
  units; with a small component it is accordingly smaller than the
  tolerance. The 3D-geometric family mixes lengths with pure numbers (its
  angular relations are normalised), and the rigid-pose family lengths
  with angles, so both keep the absolute tolerance.
- **The 2D diagnosis tells a conflict from a redundancy on small geometry.**
  Its thresholds were relative to the largest residual of a component but
  never below an absolute floor of one, so on a triangle 1e-6 across a
  distance asking for 0.2% more than the triangle has came back as
  `partially_redundant` rather than `conflicting`. The floor is now the
  size of the component when that is below one (the same unit the solve
  uses). The case is now right at every size from 1e-6 to 1e3 and at up
  to a billion sizes from the origin (`tests/scale.rs`); the corpora are
  unchanged.
- `SolveStatus` has the new variant `ConvergedToPrecision`, between
  `Converged` and `MaxIterations`; a `match` over the statuses needs an arm
  for it. The same status is how a tolerance below what double precision
  holds at all (1e-18 at unit scale, say) now ends, instead of failing.
- **A system that cannot be solved ends as `Stalled`, and soon.** A
  conflict used to wander near its least-squares compromise until the
  iteration limit and report `MaxIterations`. It now stops when no step of
  any length lowers the residual. `MaxIterations` is left for a solve that
  was still making progress when it ran out.
- **Unsatisfied equations with no free parameter are no longer "converged".**
  A component whose parameters are all fixed used to report
  `converged = true` whatever its residual. It now reports
  `converged = false` and `Stalled` when its equations do not hold.
- **The rigid-pose solve runs on the shared Newton core**, with the step
  control of the other families.
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
- **The rigid-pose `Constraint` is `{ id, kind }`**, where `kind` is the new
  `ConstraintKind` enum and carries exactly the data its kind needs. The
  struct of optional fields is gone (`constraint_type`, `component1`,
  `point1`, `component2`, `point2`, `distance`, `angle`,
  `fixed_orientation`, `dir1`, `dir2`, `lock_translation`,
  `lock_rotation`), and with it the directions that were stored in fields
  named `point`. The constructors (`Constraint::fixed`, `coincident`,
  `parallel`, …) keep their signatures, so code that builds constraints
  with them needs no change. Code that read the fields calls
  `constraint_type()`, `component1()` or `component2()`, or matches on
  `kind`. Results are the same to the last bit.
- **A rigid-pose solve that does not converge puts the bodies back too.**
  The bodies of a cluster that fails are returned where the solve started
  (the warm start, or the seed for `solve`), while clusters that converged
  keep their solution; grounded bodies stay at their fixed pose. It used to
  return the last iterate, a least-squares compromise.
  `PoseSolveOptions { on_failure: FailurePolicy::KeepLast, .. }` restores
  that. `PoseSolveResult::error` is the error of the poses returned.
- The rigid-pose options type is `rigid3d::PoseSolveOptions`. The old name
  `rigid3d::SolveOptions` remains as a deprecated alias: it was easy to
  confuse with `rcad_gcs::SolveOptions`, the options of the other two
  families. It has a new field, `on_failure`, and `PoseSolveResult` a new
  field, `status` (code that builds either by naming every field needs
  them).
- **`rigid3d::DOFAnalyzer` and `rigid3d::SolverDiagnostics` are deprecated,
  and now give the answers of the rank diagnosis.** Both used to work from
  the constraint list alone: the first subtracted a fixed number of
  freedoms per constraint, the second matched patterns. They disagreed
  with `diagnose_assembly` whenever constraints overlapped. Each now solves
  the assembly from its seed pose and reads the result there, so their
  results change:
  - `DOFResult::total_dof` is the rank-based count; `free_dofs` lists the
    coordinates of each body that can actually move and `component_dof` is
    how many; `is_over_constrained` means that a constraint is redundant or
    conflicting, not that a sum exceeded six per body.
  - `SolverDiagnostics` reports a conflict or a redundancy with the ids the
    diagnosis names, one issue per dependency group, and remaining freedoms
    with the components that can move. It no longer reports
    `OverConstrained` or `CyclicDependency` (a closed loop of mates is not
    a defect), and `DiagnosticStats::dof_removed` is a rank.
  - Both cost a solve and a factorisation now, where they used to cost a
    pass over the list.
- **A `Parallel` mate is two equations and a `Concentric` mate four**, the
  number of freedoms each removes. They were three and six: the cross
  product of two directions has three components but only two of them
  independent. Solved poses agree with the old ones to within the
  tolerance, not to the last bit, and iteration counts can differ. The
  residual norm of a given pose is the same in both forms.
- `solver::solve`, `solver::solve_components`, `space3d::solver::solve_space`,
  `space3d::solver::solve_components` and `numerics::newton_solve` take the
  parameters as `&mut [f64]` instead of `&mut Vec<f64>`.
- `Sketch::add_constraint`, `SpaceSketch::add_constraint` and
  `rigid3d::ConstraintSystem::add_constraint` return the index of the new
  constraint instead of nothing. Code that used the call as the last
  expression of a `()` block or closure needs a semicolon.
- `rigid3d::ConstraintSystem` has private fields besides `constraints` and
  `num_components`, so it can no longer be built with a struct literal:
  use `ConstraintSystem::with_components`. Its counting helpers
  (`constraints_for_component`, `total_dof_removed`, `max_dof`) leave out
  suppressed and removed constraints and removed components.
- `components`, `diagnose`, `mobility` and the DOF accounts of `Sketch` and
  `SpaceSketch` ignore suppressed and removed constraints.
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
  - on `Constraint`: `references`, `value`, `set_value`, `measure`,
    `has_finite_values`; on `PointRef`: `entity`, `expected_kinds`.
- **Editing a 3D sketch in place**, with the same methods on `SpaceSketch`:
  `remove_constraint`, `set_constraint_enabled`, `replace_constraint`,
  `constraint_state`, `remove_entity`, `is_entity_removed`,
  `set_constraint_value` / `constraint_value`, `measure`,
  `topology_revision` and `validate`. A removed plane or cylinder takes its
  representation freedoms out of the DOF count with it.
  - `SketchError` is generic over the entity-kind type, so both families
    report through one enum; `space3d::SpaceSketchError` is the 3D
    instance. `validate` also reports a spline control point that does not
    exist (`ControlPointOutOfRange`).
  - on `SpaceConstraint`: `references`, `control_point_indices`, `value`,
    `set_value`, `measure`, `has_finite_values`; on `SpacePointRef`:
    `entity`, `expected_kinds`, `point_coords`.
  - Angles are measured as principal values (`[0, π]` between two lines or
    two planes, `[0, π/2]` between a line and a plane), because that is all
    the angle equations themselves distinguish.
- **Editing an assembly in place**, on `rigid3d::ConstraintSystem`, with
  the methods of the two sketches: `remove_constraint`,
  `set_constraint_enabled`, `replace_constraint`, `constraint_state`,
  `set_constraint_value` / `constraint_value`, `measure` (given the poses),
  `topology_revision` and `validate`, plus `add_component`,
  `remove_component`, `is_component_removed`, `live_components`,
  `is_active` and `active_constraints`. A suppressed or removed constraint
  has no equation; a removed component is not solved for and counts no
  freedom, while its slot in the pose arrays stays.
  - `validate` reports `rigid3d::AssemblyError`: a component that does not
    exist or was removed, a non-finite number, a direction or rotation of
    zero length, one component on both sides of a constraint.
  - on `ConstraintKind`: `value`, `set_value`, `measure`,
    `has_finite_values`, `has_usable_directions`.
- `rigid3d::assembly_mobility`: which of each body's six coordinates can
  still move at a pose.
- `rigid3d::DiagnosticSeverity`, `DiagnosticRecommendation` and
  `DiagnosticStats` are exported. They were the types of public fields but
  could not be named from outside the crate.
- `serde` support for the rigid-pose model: `rigid3d::Constraint`,
  `ConstraintKind`, `ConstraintType` and `ConstraintSystem`. They are also
  `PartialEq` now, and `Constraint` is `Copy`.
- `rigid3d::Constraint::new`, and on `ConstraintKind`: `constraint_type`,
  `bodies`, `map_bodies`, `dof_removed`.
- `Sketch::structural_dof_estimate` and `SpaceSketch::structural_dof_estimate`.
- `rigid3d::AssemblyDiagnosis::rank` and `dependency_groups` (by mate id):
  the assembly diagnosis now reports what the 2D and 3D ones do.
- A criterion benchmark harness (`cargo bench`, `benches/solve.rs`) and the
  recorded numbers in `benches/BASELINE.md`.
- `.gitattributes` pinning LF line endings.
- Tests for the failure policy (`tests/failure_policy.rs`) and for the DOF
  account (`tests/dof_account.rs`).

### Fixed

- **`LineLength` on a line of zero length divided zero by zero.** Its
  derivative was `dx / |d|` with no guard, so the solve ended at once as
  `Diverged`, the status for numbers that are not numbers. A line of zero
  length has no direction; its length now has no gradient there, as a
  distance between coincident points already had none, and the solve ends
  as `Stalled` with the geometry put back. `tests/degenerate.rs` surveys
  sixteen degenerate starts (zero-length lines, zero radii, coincident
  points, NaN and infinite coordinates and dimensions) in all three
  families: none panics or leaves a NaN behind, and a failed solve leaves
  the geometry as it was; this was the one wrong status.
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
- **`diagnose_assembly` reported every `Parallel` and `Concentric` mate as
  redundant**, including one that was alone in the assembly, and listed it
  in a dependency group. A lone `Concentric` also put the `Fixed` of the
  body it was mated to on the redundant list. The spare rows of those mates
  were the cause; with as many rows as freedoms the reports are gone, and a
  genuinely repeated mate is still found (`tests/assembly_rank.rs`).
- **Debug builds panicked on large coordinates and on small geometry.**
  The Jacobian cross-check that debug builds run at the start of a solve
  differenced with a fixed step of 1e-7. That is a single float spacing at
  a coordinate of 1e9 and a tenth of a feature 1e-6 across, and in both
  cases the difference disagreed with a correct Jacobian and the check
  panicked. The check now allows for the rounding of the residuals, puts a
  mismatch to a second difference at half the step, and only panics when
  the two differences agree with each other against the analytic value.
  The finite-difference Jacobian of `numerics::newton_solve` widens its
  step where 1e-7 is too few float spacings to difference across.
- **The orientation error of `Fixed` and `Lock` depended on the sign of a
  quaternion.** A quaternion and its negative are the same rotation, but
  the rotation vector was read off as it came: a part 0.001 rad from where
  a lock wanted it reported an error of 0.001 in one sign and 6.282 (a full
  turn less 0.001) in the other. The solves tried still ended at the right
  pose, but the `error` reported for such a pose was wrong, and
  `diagnose_assembly` called two constraints that asked for the same
  orientation in opposite signs conflicting. The rotation is now always
  taken the short way round.
- An orientation error below about 2e-8 rad read as exactly zero, because
  the conversion from a quaternion to a rotation vector gave up there. It
  is now accurate at any size.
- glam's `serde` support was compiled in whether or not this crate's `serde`
  feature was on. It now follows the feature, so a default build pulls in no
  serde at all.
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

- No Jacobian row is finite-differenced any more. The orientation rows of
  `Fixed` and `Lock` were the last ones; they now use the derivative of the
  rotation logarithm, in closed form, where central differences evaluated
  the residual twice per column. The benchmark fixtures hold a single
  constraint of that kind, so their timings do not show the difference.
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
