# rcad-gcs — geometric constraint solver

A standalone geometric constraint solver (GCS) usable outside of rcad: one
Newton-Raphson least-squares core, three constraint families.

| Family | Module | Unknowns | Typical use |
|--------|--------|----------|-------------|
| 2D sketch | [`sketch`](src/sketch.rs) | flat `(x, y)` parameter vector | sketcher constraints + dimensions |
| 3D geometric | [`space3d`](src/space3d) | flat `(x, y, z)` parameter vector | points / lines / planes / spheres / cylinders in space |
| 3D rigid pose | [`rigid3d`](src/rigid3d) | per-body 6-DoF (translation + quaternion manifold) | assembly mate solving, kinematic posing |

## How it solves

Each family reduces to `F(x) = 0` over a parameter vector:

1. residuals from the constraint set,
2. Jacobian by central finite differences (`FD_H = 1e-7`),
3. step from the Tikhonov-regularised normal equations
   `(JᵀJ + λI) Δx = Jᵀ(−F)` — under-constrained systems get the
   minimum-norm step instead of diverging,
4. damped update (step halving), iterate until the RMS residual drops below
   `RESIDUAL_TOL = 1e-10` (rigid-pose family: L2 norm < `tolerance`, with an
   adaptive λ ladder for near-singular systems).

Known scaling limit: dense numerical Jacobians — big systems (tens of free
parameters) converge but are not guaranteed to. Analytic Jacobians and the
sparse solver ladder are the next planned core upgrades.

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

Degrees of freedom: `sk.dof()` / `sk.dof_analysis()` report
free-parameters − equations (count-based; rank-based analysis is planned).

## Features

- `serde` — derive `Serialize`/`Deserialize` on the constraint data model.
  Off by default: the core depends only on `glam`.

## Constraint catalog (summary)

- **2D** (19): Fixed, Coincident, PointDistance, LineLength, Horizontal,
  Vertical, Parallel, Perpendicular, EqualLength, Angle, EqualRadius,
  PointOnCircle, PointOnLine, Tangent, CircleCircleTangent, Radius,
  ArcArcTangent, Symmetric, Concentric, Midpoint, Diameter.
- **3D geometric** (27): point/line/plane/sphere/cylinder incidence,
  distance/length/angle/radius families, tangency, symmetry.
- **3D rigid pose** (6): Fixed (±orientation), Coincident, Distance,
  Parallel, Angle, Perpendicular — with analytic point-type Jacobian columns
  and finite-differenced direction-type rows.

## Using it from rcad

This crate was extracted from the [`cad-rs/rcad`](https://github.com/cad-rs/rcad)
monorepo, where it is consumed by `rcad-sketcher` (the 2D translation layer
over the element-based sketch model) and by `rcad-assembly` (the rigid-pose
kernel, re-exported under its historical names).

The GCS itself knows nothing about sketches' element ids, mates, BRep, or
documents — hosts translate their domain models into `Constraint`s, exactly
like any third-party consumer would.

## Status and roadmap

`0.1.0` is the first standalone release. The API is usable but still moving;
see [ROADMAP.md](ROADMAP.md) for the planned improvements (analytic Jacobians
for every family, a sparse solver ladder, a public rank-based diagnosis API,
more constraint vocabulary, benchmarks, bindings, …) and
[CHANGELOG.md](CHANGELOG.md) for what changed in each release.
