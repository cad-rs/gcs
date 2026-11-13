# rcad-gcs — roadmap

Possible improvements, not commitments. Roughly ordered by leverage (numerical
core first: everything else inherits its ceiling). Each item lists what exists
today, so the gap is explicit.

The survey of other solvers behind this list, the gap analysis and the phased
plan (P0–P8) are in [docs/CATCHUP_PLAN.md](docs/CATCHUP_PLAN.md) (in Chinese).
Phases P0 (baseline fixes) and P1 (one pipeline for the three families, an
editable model) are done; see the changelog and, for what a host has to
change, [docs/MIGRATION_0.2.md](docs/MIGRATION_0.2.md) (in Chinese).

## 1. Numerical core

- **Sparse linear algebra.** `numerics.rs` accumulates `JᵀJ` over the
  non-zero Jacobian entries, but stores it densely and solves
  `(JᵀJ + λI)Δ = Jᵀ(−F)` by dense Gaussian elimination, i.e. O(n³) in the
  free parameters of one connected component. A sparse Cholesky / QR with a
  fill-reducing ordering would lift the practical size ceiling (see
  `benches/BASELINE.md` for where it sits today). Going through the normal
  equations also squares the condition number; a QR-based minimum-norm step
  would not.
- **Small geometry.** Large coordinates are handled: an equation is held
  to the tolerance or to a few rounding units of its own numbers, whichever
  is larger. Every 2D equation is now a length (the angular relations
  included), and a fully constrained model solves to rounding at any size
  (`the_triangle_solves_exactly_at_any_size`). What is left is the
  tolerance: it is absolute, so a small component with freedom left — which
  gets no convergence polish — stops as soon as its equations are within
  1e-10 model units: a triangle 1e-4 across, free to turn, is left with its
  lengths 2.7e-8 of its size out
  (`cargo test --release --test scale -- --ignored report_small_geometry_free --nocapture`).
  A host should keep features above about 1e-4 units until the tolerance
  is scaled by the size of the component. Arc angles are columns in
  radians, not lengths, and have not been looked at for scale.
- **Hard starts.** The step control is a dogleg trust region, but the
  radius is measured in the raw parameter vector, where a length and an
  angle weigh the same, and there is no homotopy/continuation path for
  starts far from any solution. A long serial chain from a large
  perturbation still runs out of iterations (`tests/hard_cases.rs` prints
  the table); solving such a chain block by block is the better answer
  there.
- **Parallel residual / Jacobian evaluation.** The residual closure is pure;
  `rayon` could parallelise the FD column sweep and the per-constraint
  residuals behind an opt-in feature.

## 2. Diagnostics and degrees of freedom

- **Remove the deprecated assembly analyzers.** Every family reports one
  rank-based DOF account. `rigid3d::DOFAnalyzer` and
  `rigid3d::SolverDiagnostics` used to count and pattern-match; they now
  wrap `diagnose_assembly` and `assembly_mobility`, are deprecated, and can
  go in the release after 0.2.
- **Diagnosis per entity.** `diagnose()` is public and reports
  redundant / conflicting / partially-redundant constraints,
  `dependency_groups` and per-parameter mobility (2D: `Sketch::mobility`;
  assemblies: `assembly_mobility`, per body coordinate; not yet for
  `SpaceSketch`). What a UI needs on top: a status per entity (well /
  under / over-defined), the number of independent motions a body has
  rather than the coordinates they touch, and a geometric description of
  the remaining freedom (translate along, rotate about).
- **Conflict localization.** Dependency groups name the constraints that
  take part in a dependency; a minimal conflicting set and a "remove this
  one" suggestion are not computed yet.
- **Diagnosis away from the current point.** The analysis linearizes at the
  current parameters, so a configuration that is singular only where it
  stands (exact tangency) reads as redundant. Evaluating at a perturbed
  witness configuration would tell the two apart.

## 3. Constraint vocabulary

- **2D:** spline / conic incidence and tangency, arc-length dimensions,
  curvature continuity (G2), equal-tangency, point-on-spline.
- **3D geometric (`space3d`):** torus and B-spline surface incidence, more
  plane/angle families, distance-to-surface.
- **Rigid pose (`rigid3d`):** the host (rcad) currently models
  tangent / symmetric / lock / gear / cam above the solver. Evaluate moving
  the purely positional ones into the solver vocabulary.
- **Expression-driven constraints:** bind a dimension to a named parameter /
  expression so a host can drive the system (rcad does this at the document
  layer; a solver-side hook would let other hosts reuse it).

## 4. API and engineering

- **Reclaim removed slots.** All three containers (`Sketch`, `SpaceSketch`,
  the rigid-pose `ConstraintSystem`) can remove and suppress constraints,
  remove entities or components, change dimension values in place, measure,
  and validate, with stable indices. Removed slots are never reclaimed: a
  `compact()` that returns the index remapping is missing, so a long
  editing session only grows.
- **Dimensions as objects.** A dimension's value lives inside its
  constraint. Sharing one value between several constraints, or binding it
  to a variable or an expression, needs a dimension table on top.
- **Warm start as a first-class API** across all three families (rigid-pose
  already has `solve_from`; 2D/3D geometric solve from the current
  parameters).
- **One options/result type.** All three families share the vocabulary
  (`FailurePolicy`, `SolveStatus`, restore on failure by default), but not
  the types: the 2D and 3D-geometric families use `SolveOptions` /
  `SolveResult`, the rigid-pose family `PoseSolveOptions` /
  `PoseSolveResult`, which measures its tolerance as an L2 norm instead of
  an RMS and returns poses.
- **A builder API** for constraints and failing constraint ids in the
  result (today `SolveStatus` says why a solve failed, not where).
- **Serde coverage audit** under the `serde` feature: every public data type
  should round-trip (the corpus tests cover part of this).
- **Target support:** validate `wasm32-unknown-unknown` (the parent project
  builds it) and document `no_std` feasibility (glam is the only dependency).
- **`f32` variant** behind a feature for memory-constrained hosts.

## 5. Performance

- **A benchmark regression gate.** `cargo bench` runs a criterion harness
  (`benches/solve.rs`) and `benches/BASELINE.md` records the numbers, but
  nothing fails a change that makes them worse yet. Wall-clock gates on
  shared CI runners are noisy; a base-vs-head comparison in one job is the
  usual answer.
- **Reuse across re-solves:** cache the Jacobian sparsity pattern and the
  symbolic factorization between solves of the same system (parametric
  re-solve is the common host pattern).

## 6. Testing and assurance

- **Property-based tests** (`proptest` / `quickcheck`) for invariants:
  convergence implies residual < tol; a solved system re-solved is a fixed
  point; removing a redundant row does not move the solution.
- **Differential testing** against FreeCAD's planegcs (a 2D corpus export
  already exists in the parent project) and, for rigid pose, against
  OndselSolver.
- **Fuzzing** the residual/Jacobian paths with degenerate geometry
  (coincident points, zero-length edges, parallel planes).

## 7. Documentation and community

- **A `mdbook`** with a worked guide: model a sketch, diagnose a conflict,
  drive an assembly; per-family examples (today the README carries one 2D
  example and doc-tests).
- **`CHANGELOG.md` discipline** and an explicit semver policy
  (0.x: breaking changes allowed with a changelog entry).
- **Publish to crates.io** + docs.rs (the manifest already carries the
  metadata; a first `0.1.x` release is the gate).

## 8. Interoperability

- **A neutral constraint-system interchange format** (the serde types are the
  starting point; a versioned JSON schema would let non-Rust hosts exchange
  systems).
- **C ABI / Python bindings** (PyO3) so the solver can be driven from
  Python tooling — the parent project already has Python bindings for its
  other kernels.
