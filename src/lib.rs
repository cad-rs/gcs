//! # rcad-gcs — geometric constraint solver
//!
//! A standalone geometric constraint solver (GCS) with three solver families
//! over one shared Newton-Raphson least-squares core:
//!
//! - **2D sketch** ([`Sketch`]): points, lines, circles, arcs — the sketcher
//!   workhorse (coincidence, tangency, dimensions, anchoring, ...).
//! - **3D geometric** ([`space3d::SpaceSketch`]): points, lines, planes,
//!   spheres, cylinders, circles, arcs and Bezier splines in space.
//! - **3D rigid pose** ([`rigid3d`]): per-body 6-DoF position + orientation
//!   (quaternion manifold) — the assembly/mate solving kernel.
//!
//! The core treats a constraint system as `F(x) = 0` over a flat parameter
//! vector. Every constraint supplies its residuals and its analytic Jacobian
//! rows; the step comes from the Tikhonov-regularised normal equations
//! `(JᵀJ + λI) Δx = Jᵀ(−F)` (minimum-norm steps for under-constrained
//! systems), and each update is line-search damped. Systems are split into
//! the connected components of the constraint graph and solved one
//! component at a time. See [`numerics`].
//!
//! A solve that does not converge leaves the geometry as it found it (per
//! component) and says why in [`SolveResult::status`]; [`SolveOptions`]
//! selects the tolerance, the iteration limit and that failure policy.
//!
//! [`Sketch::diagnose`] reports the true (rank-based) degrees of freedom and
//! the redundant / conflicting constraints; see [`diagnosis`].
//!
//! The crate depends only on `glam`; serialization of the constraint data
//! model is available behind the optional `serde` feature.
//!
//! # Quick start (2D)
//!
//! ```
//! use rcad_gcs::Sketch;
//! use rcad_gcs::constraint::Constraint;
//!
//! let mut sk = Sketch::new();
//! let p1 = sk.add_point(0.0, 0.0);
//! let p2 = sk.add_point(3.0, 4.0);
//!
//! // Fix p1 at the origin
//! sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
//! // Constrain the distance p1→p2 to exactly 5
//! sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
//!
//! let result = sk.solve();
//! assert!(result.converged);
//! ```

pub mod constraint;
pub mod diagnosis;
pub mod edit;
pub mod entity;
pub mod graph;
pub mod mobility;
pub mod numerics;
pub mod rigid3d;
pub mod sketch;
pub mod solver;
pub mod space3d;
mod system;

pub use constraint::Constraint;
pub use diagnosis::Diagnosis;
pub use edit::{ConstraintState, SketchError};
pub use entity::{Entity, EntityId, EntityKind, PointRef};
pub use numerics::{FailurePolicy, SolveOptions, SolveResult, SolveStatus};
pub use sketch::{DofReport, Sketch};
pub use space3d::SpaceSketch;

/// The curated public surface in one `use`.
pub mod prelude {
    pub use crate::constraint::Constraint;
    pub use crate::diagnosis::Diagnosis;
    pub use crate::entity::{Entity, EntityId, EntityKind, PointRef};
    pub use crate::numerics::{
        FailurePolicy, SolveOptions, SolveResult, SolveStatus, RESIDUAL_TOL,
    };
    pub use crate::sketch::Sketch;
    pub use crate::solver;
    pub use crate::space3d;
}
