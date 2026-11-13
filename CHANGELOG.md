# Changelog

All notable changes to `rcad-gcs` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the versioning
policy is [Semantic Versioning](https://semver.org/) (0.x: breaking changes
are allowed with an entry here).

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
