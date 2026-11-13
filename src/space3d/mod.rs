//! 3D geometric constraint solver.
//!
//! Analogous to the 2D [`crate::Sketch`] but operates in 3D space with
//! points, lines, planes, spheres, and cylinders over a flat parameter
//! vector.
//!
//! The 3D family mirrors the 2D solver's capabilities: an analytic Jacobian
//! ([`constraint::SpaceConstraint::jacobian`]), connected-component
//! decomposition ([`graph`]), rank-revealing diagnosis ([`diagnosis`]), and
//! per-component localized solves ([`solver::solve_components`]).
//!
//! # Quick start
//!
//! ```
//! use rcad_gcs::space3d::SpaceSketch;
//! use rcad_gcs::space3d::constraint::SpaceConstraint;
//!
//! let mut sk = SpaceSketch::new();
//! let p1 = sk.add_point(0.0, 0.0, 0.0);
//! let p2 = sk.add_point(3.0, 4.0, 0.0);
//!
//! // Fix p1 at the origin
//! sk.add_constraint(SpaceConstraint::fix_point(p1, 0.0, 0.0, 0.0));
//! // Constrain the distance p1→p2 to exactly 5
//! sk.add_constraint(SpaceConstraint::point_distance(p1, p2, 5.0));
//!
//! let result = sk.solve();
//! assert!(result.converged);
//! ```

pub mod constraint;
pub mod diagnosis;
pub mod entity;
pub mod graph;
pub mod sketch;
pub mod solver;

pub use constraint::SpaceConstraint;
pub use diagnosis::{diagnose as diagnose_space, SpaceDiagnosis};
pub use entity::{SpaceEntity, SpaceEntityId, SpaceEntityKind, SpaceFrame3, SpacePointRef};
pub use graph::constraint_components;
pub use sketch::{SpaceDofReport, SpaceSketch, SpaceSketchError};
pub use solver::{solve_components, solve_space};
