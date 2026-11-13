# rcad-gcs — geometric constraint solver

A standalone geometric constraint solver (GCS) usable outside of rcad: one
Newton-Raphson least-squares core, three constraint families.

| Family | Module | Unknowns | Typical use |
|--------|--------|----------|-------------|
| 2D sketch | [`sketch`](src/sketch.rs) | flat `(x, y)` parameter vector | sketcher constraints + dimensions |
| 3D geometric | [`space3d`](src/space3d) | flat `(x, y, z)` parameter vector | points / lines / planes / spheres / cylinders / circles / arcs / Bezier splines in space |
| 3D rigid pose | [`rigid3d`](src/rigid3d) | per-body 6-DoF (translation + quaternion manifold) | assembly mate solving, kinematic posing |

## How it solves

Each family reduces to `F(x) = 0` over a parameter vector:

1. the constraint graph is split into connected components, and each
   component is solved on its own (so an edit re-solves only the components
   it touches);
2. residuals and **analytic** Jacobian rows come from the constraint set —
   a CI test keeps every derivative in lockstep with its residual by
   comparing against finite differences; no row is finite-differenced in
   a solve;
3. the step comes from the Tikhonov-regularised normal equations
   `(JᵀJ + λI) Δx = Jᵀ(−F)` — under-constrained systems get the
   minimum-norm step instead of diverging;
4. the step is taken inside a trust radius (Powell's dogleg): whole while
   it fits, otherwise bent toward steepest descent; a step that does not
   lower the residual is discarded and the radius shrunk;
5. iterate until the RMS residual drops below `RESIDUAL_TOL = 1e-10`
   (rigid-pose family: L2 norm < `tolerance`).

A solve that does not converge reports why in `SolveResult::status`
(`MaxIterations`, `Stalled`, `Singular`, `Diverged`) and, by default, puts
the starting geometry back — per component, so components that did converge
keep their solution. `SolveOptions` sets the tolerance, the iteration limit
and that failure policy:

```rust
use rcad_gcs::{FailurePolicy, SolveOptions};

let options = SolveOptions {
    max_iterations: 50,
    on_failure: FailurePolicy::KeepLast, // inspect how far a failed solve got
    ..SolveOptions::default()
};
let result = sk.solve_with(&options);
```

The rigid-pose family follows the same contract: `PoseSolveResult::status`
says why a solve ended, the bodies of a cluster that did not converge are
handed back where they started, and `PoseSolveOptions::on_failure` chooses
otherwise.

Known scaling limit: the linear algebra is dense. The normal equations are
solved by Gaussian elimination, cubic in the number of free parameters of
one connected component: on the reference machine a component of about 200
parameters solves in about 8 ms and one of about 800 in about 0.36 s (see
[benches/BASELINE.md](benches/BASELINE.md)). Independent components do not
compound — 200 separate rectangles (3600 parameters) solve in about 9 ms.
A sparse factorisation is the next planned core upgrade.

## Quick start

```rust
use rcad_gcs::prelude::*;

let mut sk = Sketch::new();
let p1 = sk.add_point(0.0, 0.0);
let p2 = sk.add_point(3.0, 4.0);

// Anchor p1 at the origin, lock the distance p1→p2 to 5.
sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));

let result = sk.solve();
assert!(result.converged);
```

## Editing a sketch

A sketch is edited in place; nothing has to be rebuilt between solves. The
2D `Sketch` and the 3D `SpaceSketch` have the same editing methods. Entity
ids and constraint indices are stable: removing something leaves its slot
behind instead of shifting what comes after it.

```rust
let width = sk.add_constraint(Constraint::LineLength { line, length: 40.0 });

sk.set_constraint_value(width, 55.0);     // drive the dimension
sk.solve();                               // warm start from where it stands

sk.set_constraint_enabled(width, false);  // now a reference dimension
let current = sk.measure(width);          // what the length is right now

sk.remove_constraint(width);
let gone = sk.remove_entity(line);        // and every constraint on it

sk.validate()?;                           // wrong kinds, dangling references
```

`sk.topology_revision()` changes only when the structure does (something
added, removed, enabled, suppressed or fixed), not when a value or the
geometry changes, so structure-derived data can be cached against it.

An assembly (`rigid3d::ConstraintSystem`) is edited the same way, with
`remove_component` in place of `remove_entity`. It holds no poses, so
`measure` takes them, and a re-solve continues from the last result through
`ConstraintSolver::solve_from`.

## Degrees of freedom and diagnosis

`sk.diagnose()` runs a rank-revealing QR per connected component and
reports:

- `dof` — the true degrees of freedom (free parameters − Jacobian rank);
- `redundant` / `conflicting` / `partially_redundant` constraint indices;
- `dependency_groups` — for each dependency, every constraint taking part
  in it (not just one representative).

`sk.dof()` is that same rank-based number, and `sk.mobility()` tells which
individual parameters can still move. `sk.structural_dof_estimate()` is the
cheap count (free parameters − equations); it is exact only when all
equations are independent.

The 3D family discounts storage parameters that are not physical freedoms
(a plane has four for three, a cylinder seven for five).

For assemblies, `rigid3d::diagnose_assembly` is the same diagnosis over six
columns per body, reported by mate id, and `rigid3d::assembly_mobility`
tells which of each body's six coordinates can still move.

## Features

- `serde` — derive `Serialize`/`Deserialize` on the constraint data model
  of all three families. Off by default: the core depends only on `glam`.

## Constraint catalog (summary)

- **2D** (24): Coincident, Fixed, PointDistance, DistanceX, DistanceY,
  LineLength, Horizontal, Vertical, Parallel, Perpendicular, EqualLength,
  Angle, EqualRadius, PointOnCircle, PointOnLine, Tangent (line–circle),
  CircleCircleTangent, Radius, ArcArcTangent, Symmetric, SymmetricPoint,
  Concentric, Midpoint, Diameter.
- **3D geometric** (43): point / line / plane / sphere / cylinder /
  circle / arc / spline incidence, distance / length / angle / radius
  families, tangency, symmetry, midpoint.
- **3D rigid pose** (9): Fixed (± orientation), Coincident, Distance,
  Angle, Parallel, Perpendicular, Concentric, PointOnPlane, Lock.

## Using it from rcad

This crate was extracted from the [`cad-rs/rcad`](https://github.com/cad-rs/rcad)
monorepo, where it is consumed by `rcad-sketcher` (the 2D translation layer
over the element-based sketch model) and by `rcad-assembly` (the rigid-pose
kernel, re-exported under its historical names).

The GCS itself knows nothing about sketches' element ids, mates, BRep, or
documents — hosts translate their domain models into `Constraint`s, exactly
like any third-party consumer would.

## Development

```sh
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features     # about a minute; most of it is the assembly corpus
cargo bench                   # criterion benchmarks, see benches/BASELINE.md
```

## Status and roadmap

The API is usable but still moving; see [ROADMAP.md](ROADMAP.md) for the
planned improvements and [CHANGELOG.md](CHANGELOG.md) for what changed in
each release.
