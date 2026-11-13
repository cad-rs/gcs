//! Rank-revealing diagnosis for assembly constraint systems (solver
//! roadmap M6.2) — the 2D [`crate::diagnosis`] machinery applied to the
//! rigid-pose Jacobian.
//!
//! The assembly Jacobian (m equations × 6·n pose columns) is analyzed by a
//! column-pivoted QR of its transpose, per BODY-CLUSTER of the constraint
//! graph (constraints joining different bodies are trivially independent,
//! so per-cluster QRs are equivalent to one global QR while scaling with
//! the cluster sizes). Dependent rows are classified by whether their
//! residual agrees with the least-squares step of the independent rows:
//!
//! - **redundant** — dependent and consistent (a duplicate or implied
//!   mate; removable without changing the assembly's solution);
//! - **conflicting** — dependent and inconsistent (no pose satisfies the
//!   whole system at this linearization);
//! - **partially_redundant** — the gray band, informational.
//!
//! `dof` is the true assembly DOF: 6·num_components − Σ rank over
//! clusters (bodies claimed by no constraint count their 6 fully).
//!
//! **Pose discipline**: rank is pose-dependent, so the caller chooses the
//! pose deliberately — redundancy is diagnosed at a TIGHTLY-SOLVED pose
//! (solve with a small tolerance first: dependent-row consistency is only
//! meaningful near the solution), conflicts at the authored/failed pose
//! where the contradiction is visible in the linearization. Same
//! discipline as the 2D diagnosis.

use crate::rigid3d::equations::{pack, RigidEquations, COLUMNS_PER_BODY};
use crate::rigid3d::{ConstraintSystem, Point3};
use crate::system::{self, NoGauge};
use glam::DQuat;

pub use crate::diagnosis::{CONSIST_TOL, PARTIAL_TOL, RANK_REL_TOL};

/// Structured diagnosis of an assembly constraint system. The reported
/// constraint identifiers are each `Constraint`'s `id` field.
#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AssemblyDiagnosis {
    /// True assembly DOF: 6·num_components − Σ cluster rank.
    pub dof: i64,
    /// Rank of the constraint Jacobian over all clusters
    /// (`6·num_components − dof`).
    #[cfg_attr(feature = "serde", serde(default))]
    pub rank: usize,
    /// Mate ids that are dependent but consistent.
    pub redundant: Vec<usize>,
    /// Mate ids whose demands are inconsistent with the rest.
    pub conflicting: Vec<usize>,
    /// Dependent rows in the gray band. Informational.
    pub partially_redundant: Vec<usize>,
    /// Mutually-dependent mate groups, by id: each group is every mate that
    /// takes part in one linear dependency among the rows — the answer to
    /// "which mates are redundant together" that does not depend on which
    /// member the pivoting happened to single out.
    #[cfg_attr(feature = "serde", serde(default))]
    pub dependency_groups: Vec<Vec<usize>>,
}

/// Diagnose an assembly at the given pose (see the pose discipline in the
/// module docs).
///
/// This is the shared diagnosis of the crate-internal `system` module run
/// over the rigid-pose equations: six columns per body, nothing held.
pub fn diagnose(
    system: &ConstraintSystem,
    positions: &[Point3],
    rots: &[DQuat],
) -> AssemblyDiagnosis {
    let params = pack(positions, rots);
    let fixed = vec![false; system.num_components * COLUMNS_PER_BODY];
    let d = system::diagnose(&RigidEquations::new(system), &params, &fixed, &NoGauge);

    // Constraint indices → mate ids, keeping the order and dropping repeats
    // (two constraints may carry the same id).
    let ids = |indices: &[usize]| -> Vec<usize> {
        let mut out: Vec<usize> = Vec::with_capacity(indices.len());
        for &i in indices {
            let id = system.constraints[i].id;
            if !out.contains(&id) {
                out.push(id);
            }
        }
        out
    };
    AssemblyDiagnosis {
        dof: d.dof,
        rank: d.rank,
        redundant: ids(&d.redundant),
        conflicting: ids(&d.conflicting),
        partially_redundant: ids(&d.partially_redundant),
        dependency_groups: d.dependency_groups.iter().map(|g| ids(g)).collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rigid3d::{Constraint, ConstraintSolver, ConstraintSystem};
    use glam::DVec3;

    const TIGHT_TOL: f64 = 1.0e-12;

    fn tight_solver() -> ConstraintSolver {
        ConstraintSolver::with_options(crate::rigid3d::PoseSolveOptions {
            tolerance: TIGHT_TOL,
            ..crate::rigid3d::PoseSolveOptions::default()
        })
    }

    // ── seed-fault corpus (M6.2 acceptance): 10/10/10 ────────────────────

    fn build(kind: usize, v: f64) -> ConstraintSystem {
        let mut s = ConstraintSystem::with_components(3);
        match kind {
            // ── conflicts 0..=9: linearly visible (duplicate structures
            //    disagreeing on values) ──
            // 0-3: double fixed with disagreeing positions (all offsets
            // non-zero: a zero offset would be CONSISTENT-redundant).
            0..=3 => {
                let (dx, dy, dz) = [
                    (0.5, 0.5, 0.5),
                    (1.0, 1.0, 1.0),
                    (2.0, 0.0, 0.0),
                    (0.0, 0.0, 3.0),
                ][kind % 4];
                s.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
                s.add_constraint(Constraint::fixed(1, 0, DVec3::new(dx, dy, dz)));
            }
            // 4-6: duplicate distance with different values.
            4..=6 => {
                s.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
                s.add_constraint(Constraint::distance(
                    1,
                    0,
                    DVec3::ZERO,
                    1,
                    DVec3::new(2.0, 0.0, 0.0),
                    2.0,
                ));
                s.add_constraint(Constraint::distance(
                    2,
                    0,
                    DVec3::ZERO,
                    1,
                    DVec3::new(2.0, 0.0, 0.0),
                    2.0 + v * 2.0,
                ));
            }
            // 7-9: duplicate angle with different values.
            7..=9 => {
                s.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
                s.add_constraint(Constraint::angle(
                    1,
                    0,
                    DVec3::new(1.0, 0.0, 0.0),
                    1,
                    DVec3::new(1.0, 0.0, 0.0),
                    0.5,
                ));
                s.add_constraint(Constraint::angle(
                    2,
                    0,
                    DVec3::new(1.0, 0.0, 0.0),
                    1,
                    DVec3::new(1.0, 0.0, 0.0),
                    0.5 + v,
                ));
            }
            _ => unreachable!(),
        }
        s
    }

    #[test]
    fn conflict_seed_set_all_flagged() {
        for kind in 0..10usize {
            let v = match kind {
                4..=6 => (kind - 3) as f64 * 0.5,
                7..=9 => (kind - 6) as f64 * 0.3,
                _ => 0.0,
            };
            let system = build(kind, v);
            let (pos, rot) = ConstraintSolver::seed_fixed(&system);
            let d = diagnose(&system, &pos, &rot);
            assert!(
                !d.conflicting.is_empty(),
                "kind {kind} v {v}: expected a conflict, got {d:?}"
            );
        }
    }

    #[test]
    fn redundant_seed_set_all_flagged_no_false_conflicts() {
        // Ten consistent-redundant assemblies. Every case is solved TIGHT
        // first (pose discipline), then diagnosed at the solved pose.
        let build = |kind: usize| -> ConstraintSystem {
            let mut s = ConstraintSystem::with_components(2);
            s.add_constraint(Constraint::fixed_with_orientation(
                0,
                0,
                DVec3::ZERO,
                DQuat::IDENTITY,
            ));
            match kind {
                // 0-1: duplicate coincident.
                0 | 1 => {
                    s.add_constraint(Constraint::coincident(
                        1,
                        0,
                        DVec3::ZERO,
                        1,
                        DVec3::new(1.0, 0.5, 0.0),
                    ));
                    s.add_constraint(Constraint::coincident(
                        2,
                        0,
                        DVec3::ZERO,
                        1,
                        DVec3::new(1.0, 0.5, 0.0),
                    ));
                }
                // 2-3: duplicate distance, same value.
                2 | 3 => {
                    s.add_constraint(Constraint::distance(
                        1,
                        0,
                        DVec3::ZERO,
                        1,
                        DVec3::new(2.0, 0.0, 0.0),
                        2.0,
                    ));
                    s.add_constraint(Constraint::distance(
                        2,
                        0,
                        DVec3::ZERO,
                        1,
                        DVec3::new(2.0, 0.0, 0.0),
                        2.0,
                    ));
                }
                // 4-5: duplicate angle, same value. Body 1's local dir is
                // NOT parallel to body 0's — an angle between currently
                // parallel dirs has a zero first-order gradient (rotating a
                // vector does not change its dot with itself) and the tight
                // solve stalls.
                4 | 5 => {
                    s.add_constraint(Constraint::angle(
                        1,
                        0,
                        DVec3::new(1.0, 0.0, 0.0),
                        1,
                        DVec3::new(0.0, 1.0, 0.0),
                        0.6,
                    ));
                    s.add_constraint(Constraint::angle(
                        2,
                        0,
                        DVec3::new(1.0, 0.0, 0.0),
                        1,
                        DVec3::new(0.0, 1.0, 0.0),
                        0.6,
                    ));
                }
                // 6-7: angle 90° + perpendicular (identical demands).
                6 | 7 => {
                    s.add_constraint(Constraint::angle(
                        1,
                        0,
                        DVec3::new(1.0, 0.0, 0.0),
                        1,
                        DVec3::new(0.0, 1.0, 0.0),
                        std::f64::consts::FRAC_PI_2,
                    ));
                    s.add_constraint(Constraint::perpendicular(
                        2,
                        0,
                        DVec3::new(1.0, 0.0, 0.0),
                        1,
                        DVec3::new(0.0, 1.0, 0.0),
                    ));
                }
                // 8-9: duplicate parallel.
                _ => {
                    s.add_constraint(Constraint::parallel(
                        1,
                        0,
                        DVec3::new(0.0, 1.0, 0.0),
                        1,
                        DVec3::new(0.0, 1.0, 0.0),
                    ));
                    s.add_constraint(Constraint::parallel(
                        2,
                        0,
                        DVec3::new(0.0, 1.0, 0.0),
                        1,
                        DVec3::new(0.0, 1.0, 0.0),
                    ));
                }
            }
            s
        };

        let solver = tight_solver();
        for kind in 0..10usize {
            let system = build(kind);
            let r = solver.solve(&system);
            assert!(r.success, "kind {kind}: redundant case must converge");
            let d = diagnose(&system, &r.positions, &r.rotations);
            assert!(
                !d.redundant.is_empty(),
                "kind {kind}: expected redundant rows, got {d:?}"
            );
            assert!(
                d.conflicting.is_empty(),
                "kind {kind}: redundant misclassified as conflicting: {d:?}"
            );
        }
    }

    #[test]
    fn under_seed_set_all_positive_dof() {
        let build = |kind: usize| -> ConstraintSystem {
            let mut s = ConstraintSystem::with_components(2);
            match kind {
                // 0-1: a single free body, no constraints.
                0 | 1 => {}
                // 2-3: one coincident between two free bodies.
                2 | 3 => {
                    s.add_constraint(Constraint::coincident(
                        0,
                        0,
                        DVec3::ZERO,
                        1,
                        DVec3::new(1.0, 0.0, 0.0),
                    ));
                }
                // 4-5: one distance between two free bodies.
                4 | 5 => {
                    s.add_constraint(Constraint::distance(
                        0,
                        0,
                        DVec3::ZERO,
                        1,
                        DVec3::new(2.0, 0.0, 0.0),
                        2.0,
                    ));
                }
                // 6-7: one concentric between two free bodies.
                6 | 7 => {
                    s.add_constraint(Constraint::concentric(
                        0,
                        0,
                        DVec3::ZERO,
                        DVec3::new(0.0, 0.0, 1.0),
                        1,
                        DVec3::new(2.0, 0.0, 0.0),
                        DVec3::new(0.0, 0.0, 1.0),
                    ));
                }
                // 8-9: coincident + parallel.
                _ => {
                    s.add_constraint(Constraint::coincident(
                        0,
                        0,
                        DVec3::ZERO,
                        1,
                        DVec3::new(1.0, 0.0, 0.0),
                    ));
                    s.add_constraint(Constraint::parallel(
                        1,
                        0,
                        DVec3::new(0.0, 1.0, 0.0),
                        1,
                        DVec3::new(0.0, 1.0, 0.0),
                    ));
                }
            }
            s
        };

        for kind in 0..10usize {
            let system = build(kind);
            let (pos, rot) = ConstraintSolver::seed_fixed(&system);
            let d = diagnose(&system, &pos, &rot);
            assert!(
                d.dof > 0 && d.conflicting.is_empty(),
                "kind {kind}: expected under-constrained, got {d:?}"
            );
        }
    }

    /// A shaft in a bore: coincident + concentric is under-constrained by
    /// exactly the spin about the axis (dof 1 on the two-body assembly
    /// with a grounded base) — the mate-redundancy-adjacent DOF report.
    #[test]
    fn shaft_in_bore_reports_spin_dof() {
        let mut s = ConstraintSystem::with_components(2);
        s.add_constraint(Constraint::fixed_with_orientation(
            0,
            0,
            DVec3::ZERO,
            DQuat::IDENTITY,
        ));
        s.add_constraint(Constraint::coincident(
            1,
            0,
            DVec3::new(1.0, 0.0, 1.0),
            1,
            DVec3::new(0.0, 0.0, 0.0),
        ));
        s.add_constraint(Constraint::concentric(
            2,
            0,
            DVec3::new(1.0, 0.0, 1.0),
            DVec3::new(0.0, 0.0, 1.0),
            1,
            DVec3::ZERO,
            DVec3::new(0.0, 0.0, 1.0),
        ));
        let (pos, rot) = ConstraintSolver::seed_fixed(&s);
        let d = diagnose(&s, &pos, &rot);
        assert_eq!(d.dof, 1, "the spin about the shaft axis: {d:?}");
        assert!(d.conflicting.is_empty(), "{d:?}");
    }
}
