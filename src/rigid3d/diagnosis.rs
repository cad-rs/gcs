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

use crate::numerics::{householder_qr_column_pivoting, solve_normal_equations};
use crate::rigid3d::solver::ConstraintSolver;
use crate::rigid3d::{ConstraintSystem, Point3};
use glam::DQuat;

/// Rank threshold, relative to the largest initial column norm of the
/// row-equilibrated transpose Jacobian. All mate rows are analytic (or
/// self-consistent FD for ground-orientation rows), so this sits well
/// above roundoff.
pub const RANK_REL_TOL: f64 = 1.0e-9;

/// Dependent-row consistency threshold (relative to the residual scale).
pub const CONSIST_TOL: f64 = 1.0e-9;

/// Upper end of the gray band between redundant and conflicting.
pub const PARTIAL_TOL: f64 = 1.0e-6;

/// Structured diagnosis of an assembly constraint system. The reported
/// constraint identifiers are each `Constraint`'s `id` field.
#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AssemblyDiagnosis {
    /// True assembly DOF: 6·num_components − Σ cluster rank.
    pub dof: i64,
    /// Mate ids that are dependent but consistent.
    pub redundant: Vec<usize>,
    /// Mate ids whose demands are inconsistent with the rest.
    pub conflicting: Vec<usize>,
    /// Dependent rows in the gray band. Informational.
    pub partially_redundant: Vec<usize>,
}

/// Diagnose an assembly at the given pose (see the pose discipline in the
/// module docs).
pub fn diagnose(
    system: &ConstraintSystem,
    positions: &[Point3],
    rots: &[DQuat],
) -> AssemblyDiagnosis {
    let free_total = (system.num_components * 6) as i64;
    let n_eq = ConstraintSolver::equation_count(system);
    if n_eq == 0 {
        return AssemblyDiagnosis {
            dof: free_total,
            ..Default::default()
        };
    }

    let clusters = ConstraintSolver::components(system);
    let mut d = AssemblyDiagnosis::default();
    let mut total_rank = 0i64;

    for group in &clusters {
        // The cluster as its own system over only the bodies it references
        // (sorted, so the column mapping is stable): residuals and Jacobian
        // rows depend on no other body, and an untouched body would only
        // add zero columns.
        let sub_system = ConstraintSystem {
            constraints: group
                .iter()
                .map(|&i| system.constraints[i].clone())
                .collect(),
            num_components: system.num_components,
        };
        let bodies = ConstraintSolver::referenced_bodies(&sub_system);
        let compact = ConstraintSolver::compact_system(&sub_system, &bodies);
        let cluster_positions: Vec<Point3> = bodies.iter().map(|&b| positions[b]).collect();
        let cluster_rots: Vec<DQuat> = bodies.iter().map(|&b| rots[b]).collect();
        let rj =
            ConstraintSolver::build_residual_jacobian(&compact, &cluster_positions, &cluster_rots);

        let n_eq_c = rj.r.len();
        let n_free_c = bodies.len() * 6;
        // Residuals and Jacobian in cluster row order.
        let mut f_c = rj.r;
        let mut j = rj.j;

        // Row-equilibrate (J and f together) so a distance in mm and an
        // angle in radians contribute comparable column norms.
        for (row_j, fi) in j.iter_mut().zip(f_c.iter_mut()) {
            let s = row_j.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            if s > 0.0 {
                for v in row_j.iter_mut() {
                    *v /= s;
                }
                *fi /= s;
            }
        }

        // A = Jᵀ: column i of A is cluster row i (equilibrated).
        let mut a = vec![vec![0.0_f64; n_eq_c]; n_free_c];
        for (i, row) in j.iter().enumerate() {
            for (p, v) in row.iter().enumerate() {
                a[p][i] = *v;
            }
        }
        let (pivots, rank) = householder_qr_column_pivoting(&mut a, RANK_REL_TOL);
        total_rank += rank as i64;

        let independent = &pivots[..rank];
        let dependent: Vec<usize> = (0..n_eq_c).filter(|&i| !independent.contains(&i)).collect();
        if dependent.is_empty() || rank == 0 {
            continue;
        }

        // Least-squares step of the independent rows.
        let j_sub: Vec<Vec<f64>> = independent.iter().map(|&i| j[i].clone()).collect();
        let f_sub: Vec<f64> = independent.iter().map(|&i| f_c[i]).collect();
        let z = solve_normal_equations(&j_sub, &f_sub, n_free_c);

        let scale = f_c.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        // Row → constraint id within this cluster.
        let mut row_id: Vec<usize> = Vec::with_capacity(n_eq_c);
        for &gi in group {
            let eqs = ConstraintSolver::equation_count_of(&system.constraints[gi]);
            for _ in 0..eqs {
                row_id.push(system.constraints[gi].id);
            }
        }
        for &r in &dependent {
            let drive: f64 = j[r].iter().zip(z.iter()).map(|(ji, zi)| ji * zi).sum();
            let e = (f_c[r] + drive).abs();
            let global = row_id[r];
            let bucket: &mut Vec<usize> = if e <= CONSIST_TOL * scale {
                &mut d.redundant
            } else if e <= PARTIAL_TOL * scale {
                &mut d.partially_redundant
            } else {
                &mut d.conflicting
            };
            if !bucket.contains(&global) {
                bucket.push(global);
            }
        }
    }
    d.dof = free_total - total_rank;
    d
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rigid3d::{Constraint, ConstraintSolver, ConstraintSystem};
    use glam::DVec3;

    const TIGHT_TOL: f64 = 1.0e-12;

    fn tight_solver() -> ConstraintSolver {
        ConstraintSolver::with_options(crate::rigid3d::SolveOptions {
            tolerance: TIGHT_TOL,
            ..crate::rigid3d::SolveOptions::default()
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
