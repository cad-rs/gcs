# rcad-gcs — roadmap

Possible improvements, not commitments. Roughly ordered by leverage (numerical
core first: everything else inherits its ceiling). Each item lists what exists
today, so the gap is explicit.

## 1. Numerical core

- **Analytic Jacobians for every family.** Today only the rigid-pose
  point-type rows are analytic; the 2D and 3D-geometric families (and the
  rigid-pose direction rows) use central finite differences (`FD_H = 1e-7`).
  Analytic rows are faster and free of truncation error — the rigid-pose
  family already proves the pattern (see its analytic-vs-FD gate test).
- **Sparse / iterative linear algebra.** `numerics.rs` solves
  `(JᵀJ + λI)Δ = Jᵀ(−F)` with dense Gaussian elimination on `Vec`, i.e.
  O(n³) in the parameter count. A sparse Cholesky / CG-SMG ladder would lift
  the practical size ceiling (today: "tens of free parameters converge but
  are not guaranteed to").
- **Better step control.** The core uses Tikhonov λ + Armijo-style halving.
  A Levenberg–Marquardt trust region, automatic scaling/equilibration of the
  normal equations, and a homotopy/continuation path for hard starts would
  raise the convergence rate on badly conditioned systems.
- **Parallel residual / Jacobian evaluation.** The residual closure is pure;
  `rayon` could parallelise the FD column sweep and the per-constraint
  residuals behind an opt-in feature.

## 2. Diagnostics and degrees of freedom

- **Unify the two DOF accounts.** `Sketch::dof_analysis` is count-based
  (free parameters − equations) while `diagnosis` already reports the true
  rank-based DOF (free parameters − Jacobian rank). One public account, with
  the rank analysis as the authority.
- **Expose the redundancy analysis.** `diagnosis.rs` classifies
  redundant / conflicting / partially-redundant rows and
  `dependency_groups` names every member of a dependency; this is currently
  internal to the solve path. A public `diagnose(&system)` returning
  typed, name-mapped results is the highest-leverage diagnostic upgrade.
- **Conflict localization.** Report *which* constraints conflict (ids/names),
  not just that the system is inconsistent — the assembly side already maps
  rows back to mate names and can be the model.

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

- **Warm start as a first-class API** across all three families (rigid-pose
  already has `solve_from`; 2D/3D geometric do not).
- **A builder API** for constraints and a typed error enum (today callers get
  `SolveResult` + diagnostics; a structured `SolveError` would carry the
  failing constraint ids).
- **Serde coverage audit** under the `serde` feature: every public data type
  should round-trip (the corpus tests cover part of this).
- **Target support:** validate `wasm32-unknown-unknown` (the parent project
  builds it) and document `no_std` feasibility (glam is the only dependency).
- **`f32` variant** behind a feature for memory-constrained hosts.

## 5. Performance

- **A criterion benchmark harness** with a baseline corpus and a regression
  gate; the existing corpora (2D convergence, 3D geometric, assembly) are the
  natural fixtures.
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
