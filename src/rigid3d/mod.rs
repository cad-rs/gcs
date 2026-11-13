//! 3D rigid-pose constraint solver (assembly/mate kernel).
//!
//! Unknowns per body: translation `(tx, ty, tz)` plus incremental orientation
//! via quaternion manifold integration (`q ← normalize(q · dq)`) — 6 columns
//! per component. Jacobian rows are analytic (translation columns ±I,
//! rotation columns via the `ω × r` linearization, direction rows via
//! `δd = ω × d`), except the orientation rows of `Fixed` and `Lock`, which
//! use central finite differences. The solve and the diagnosis are the
//! shared pipeline of the other families, run over these columns: the
//! Newton core of [`crate::numerics`] steps a pose by adding the
//! translation and rotating the quaternion. A body that a `Fixed` constraint
//! pins in full is grounded and is not solved for.
//!
//! Constraint set: Fixed (± orientation), Coincident, Distance, Angle,
//! Parallel, Perpendicular, Concentric, PointOnPlane, Lock.
//!
//! This family is deliberately at the *constraint-system* level: a host
//! (e.g. an assembly layer) resolves geometric references to body-local
//! points/directions first and then translates them into [`Constraint`]s —
//! the solver never sees BRep or document types.

mod constraint;
pub mod diagnosis;
mod diagnostics;
mod dof;
mod equations;
mod solver;

pub use constraint::{ComponentId, Constraint, ConstraintKind, ConstraintSystem, ConstraintType};
pub use diagnosis::{
    diagnose as diagnose_assembly, AssemblyDiagnosis, CONSIST_TOL, PARTIAL_TOL, RANK_REL_TOL,
};
pub use diagnostics::{DiagnosticIssue, DiagnosticType, SolverDiagnostics};
pub use dof::{DOFAnalyzer, DOFResult, DegreeOfFreedom};
#[allow(deprecated)]
pub use solver::SolveOptions;
pub use solver::{ConstraintSolver, PoseSolveOptions, PoseSolveResult};

/// 3D point (alias for `glam::DVec3`).
pub type Point3 = glam::DVec3;
