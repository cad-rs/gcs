//! 3D rigid-pose constraint solver (assembly/mate kernel).
//!
//! Unknowns per body: translation `(tx, ty, tz)` plus incremental orientation
//! via quaternion manifold integration (`q ← normalize(q · dq)`) — 6 columns
//! per component. Residual/Jacobian rows mix analytic forms (point-type
//! constraints: translation columns ±I, rotation columns via the `ω × r`
//! linearization) with central finite differences (direction-type
//! constraints). The step solves `(JᵀJ + λI) Δ = −Jᵀr` with an adaptive λ
//! ladder for near-singular systems (Levenberg–Marquardt-style), on the
//! shared Gaussian-elimination core of [`crate::numerics`].
//!
//! This family is deliberately at the *constraint-system* level: a host
//! (e.g. an assembly layer) resolves geometric references to body-local
//! points/directions first and then translates them into [`Constraint`]s —
//! the solver never sees BRep or document types.

mod constraint;
pub mod diagnosis;
mod diagnostics;
mod dof;
mod solver;

pub use constraint::{ComponentId, Constraint, ConstraintSystem, ConstraintType};
pub use diagnosis::{diagnose as diagnose_assembly, AssemblyDiagnosis, CONSIST_TOL, PARTIAL_TOL, RANK_REL_TOL};
pub use diagnostics::{DiagnosticIssue, DiagnosticType, SolverDiagnostics};
pub use dof::{DegreeOfFreedom, DOFAnalyzer, DOFResult};
pub use solver::{ConstraintSolver, PoseSolveResult, SolveOptions};

/// 3D point (alias for `glam::DVec3`).
pub type Point3 = glam::DVec3;
