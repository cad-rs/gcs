//! Rank-revealing diagnosis of a constraint system (solver roadmap M1/M2).
//!
//! Before (or instead of) iterating, the constraint Jacobian is analyzed by
//! a column-pivoted Householder QR of its transpose: the pivot columns are
//! the linearly independent constraint rows, so the complement is the set
//! of dependent constraints. Each dependent row is then classified by
//! whether its residual agrees with the least-squares step of the
//! independent rows:
//!
//! - **redundant** — dependent and consistent (adds no information; the
//!   solve can proceed and the row can be removed without changing the
//!   solution);
//! - **conflicting** — dependent and inconsistent (no solution satisfies
//!   the whole system at this linearization; the solve would thrash, so
//!   callers fail fast);
//! - **partially_redundant** — dependent with a weakly inconsistent
//!   residual (the gray band: reported but not fatal).
//!
//! `dof` is the true degrees of freedom: free parameters minus the
//! Jacobian rank (redundant rows consume none).
//!
//! Since M2 the Jacobian is the ANALYTIC one (exact up to floating-point
//! roundoff ~1e-14), so the thresholds are tight: true dependencies drop
//! to roundoff, independent directions sit orders of magnitude above.

use crate::constraint::Constraint;
use crate::entity::Entity;
use crate::numerics::{
    gaussian_elimination, householder_qr_column_pivoting, solve_normal_equations, LAMBDA,
};
use crate::solver::{analytic_jacobian, eval_residuals};
use std::collections::HashMap;

/// Column-pivot rank threshold, relative to the largest initial column
/// norm of the (row-equilibrated) transpose Jacobian. With analytic
/// derivatives an exact dependency projects to floating-point roundoff
/// (~1e-14 relative), so 1e-9 sits comfortably above it.
pub const RANK_REL_TOL: f64 = 1.0e-9;

/// Dependent-row consistency threshold (relative to the residual scale):
/// at or below this the row agrees with the independent rows' least-squares
/// step — redundant.
pub const CONSIST_TOL: f64 = 1.0e-9;

/// Upper end of the gray band between redundant and conflicting
/// (relative). Rows in between are reported as partially redundant.
pub const PARTIAL_TOL: f64 = 1.0e-6;

/// Structured diagnosis of a constraint system.
#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Diagnosis {
    /// True degrees of freedom: free parameters − Jacobian rank.
    pub dof: i64,
    /// Rank of the constraint Jacobian over the free parameters, from the
    /// same rank-revealing QR that decides `dof` (`rank = free − dof`).
    /// Deterministic: column pivoting ties break by position.
    #[cfg_attr(feature = "serde", serde(default))]
    pub rank: usize,
    /// Constraints that are dependent but consistent — removable without
    /// changing the solution.
    pub redundant: Vec<usize>,
    /// Constraints whose demands are inconsistent with the rest — no
    /// solution exists at this linearization.
    pub conflicting: Vec<usize>,
    /// Dependent rows in the gray band (weakly inconsistent; possibly FD
    /// noise, possibly genuine nonlinearity). Informational.
    pub partially_redundant: Vec<usize>,
    /// Mutually-dependent constraint groups (sketch_singular v2): each
    /// group is the full support of one linear dependency among rows,
    /// extracted from the pivoted QR's R factor — the flip-free answer
    /// to "which constraints are redundant together" (a single
    /// representative would depend on pivot tie-breaking).
    #[cfg_attr(feature = "serde", serde(default))]
    pub dependency_groups: Vec<Vec<usize>>,
}

/// Diagnose a constraint system at the current parameters.
///
/// Constraint indices are the positions in `constraints` (equations stack
/// in order, 1 or 2 rows per constraint). Since solver roadmap M3 the
/// analysis runs PER COMPONENT of the constraint graph — dependencies can
/// only occur inside a component (cross-component rows share no
/// parameters, hence are trivially independent), so the per-component QRs
/// are mathematically equivalent to one global QR while scaling with the
/// component sizes. The reported indices are global.
pub fn diagnose(
    params: &[f64],
    fixed: &[bool],
    entities: &[Entity],
    constraints: &[Constraint],
) -> Diagnosis {
    let free_total = params.iter().zip(fixed).filter(|(_, &f)| !f).count();
    let n_eq: usize = constraints.iter().map(|c| c.equation_count()).sum();

    if n_eq == 0 {
        return Diagnosis {
            dof: free_total as i64,
            ..Default::default()
        };
    }

    // Global residuals once; per-component slices below.
    let f = eval_residuals(params, entities, constraints, n_eq);
    // Row offset of each constraint in the stacked residual vector.
    let mut row_of: Vec<usize> = Vec::with_capacity(constraints.len() + 1);
    let mut acc = 0usize;
    for c in constraints {
        row_of.push(acc);
        acc += c.equation_count();
    }

    let components = crate::graph::constraint_components(params, entities, constraints);
    let mut d = Diagnosis::default();
    // True DOF = free_total − Σ rank_c: parameters claimed by no component
    // are unconstrained and count fully (e.g. an entity whose only
    // constraint pins one endpoint).
    let mut total_rank = 0i64;

    for group in &components {
        let sub: Vec<&Constraint> = group.iter().map(|&i| &constraints[i]).collect();
        let n_eq_c: usize = sub.iter().map(|c| c.equation_count()).sum();
        // Local residual slice, in component row order.
        let mut f_c = Vec::with_capacity(n_eq_c);
        for &gi in group {
            let (s, e) = (row_of[gi], row_of[gi] + constraints[gi].equation_count());
            f_c.extend_from_slice(&f[s..e]);
        }

        // Touched free parameters of this component.
        let mut touched: Vec<usize> = Vec::new();
        {
            let mut scratch: Vec<Vec<(usize, f64)>> = Vec::new();
            for c in &sub {
                scratch.clear();
                scratch.resize_with(c.equation_count(), Vec::new);
                c.jacobian(params, entities, &mut scratch);
                for row in &scratch {
                    for &(p, _) in row {
                        if !fixed[p] && !touched.contains(&p) {
                            touched.push(p);
                        }
                    }
                }
            }
            touched.sort_unstable();
        }
        let n_free_c = touched.len();

        // Map component-local rows back to global constraint indices.
        let mut row_constraint: Vec<usize> = Vec::with_capacity(n_eq_c);
        for &gi in group {
            for _ in 0..constraints[gi].equation_count() {
                row_constraint.push(gi);
            }
        }

        if n_free_c == 0 {
            // Nothing in this component can move: any unsatisfied residual
            // is a genuine conflict.
            let scale = f_c.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            for (r, fi) in f_c.iter().enumerate() {
                if fi.abs() > CONSIST_TOL * scale {
                    d.conflicting.push(row_constraint[r]);
                }
            }
            continue;
        }

        let col_of: HashMap<usize, usize> =
            touched.iter().enumerate().map(|(c, &p)| (p, c)).collect();
        let mut j = vec![vec![0.0_f64; n_free_c]; n_eq_c];
        analytic_jacobian(params, &col_of, entities, &sub, &mut j);

        // Row-equilibrate the Jacobian and the residuals together so a
        // distance in mm and an angle in radians contribute comparable
        // column norms to the rank decision.
        for (row, fi) in j.iter_mut().zip(f_c.iter_mut()) {
            let s = row.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            if s > 0.0 {
                for v in row.iter_mut() {
                    *v /= s;
                }
                *fi /= s;
            }
        }

        // A = Jᵀ: column i of A is component row i (equilibrated).
        let mut a = vec![vec![0.0_f64; n_eq_c]; n_free_c];
        for (i, row) in j.iter().enumerate() {
            for (p, v) in row.iter().enumerate() {
                a[p][i] = *v;
            }
        }
        // The QR mutates `a` in place; the group extraction needs the
        // ORIGINAL columns, so keep a copy for it.
        let a_original = a.clone();
        let (pivots, rank) = householder_qr_column_pivoting(&mut a, RANK_REL_TOL);
        total_rank += rank as i64;

        let independent = &pivots[..rank];
        let dependent: Vec<usize> = (0..n_eq_c).filter(|&i| !independent.contains(&i)).collect();
        if dependent.is_empty() || rank == 0 {
            continue;
        }

        // Null-space dependency groups (sketch_singular v2): the
        // dependency lives in ROW space — find coefficients y with
        // SUM y_j*row_j = dep_row via the row-space Gram system
        // (G + lambda*I) y = -(rows*dep), G[i][j] = row_i.row_j. The
        // coefficient support is the whole mutually-dependent row group,
        // invariant to pivot tie-breaking. The Gram solve reuses the
        // shared gaussian_elimination (the r x r system is tiny).
        let rows: Vec<Vec<f64>> = independent
            .iter()
            .map(|&col| (0..n_free_c).map(|rr| a_original[rr][col]).collect())
            .collect();
        let r = rows.len();
        for &dep in &dependent {
            let dep_row: Vec<f64> = (0..n_free_c).map(|rr| a_original[rr][dep]).collect();
            let mut gram = vec![vec![0.0_f64; r]; r];
            let mut rhs = vec![0.0_f64; r];
            for (i, ri) in rows.iter().enumerate() {
                rhs[i] = ri
                    .iter()
                    .zip(dep_row.iter())
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
                for (j, rj) in rows.iter().enumerate() {
                    gram[i][j] = ri.iter().zip(rj.iter()).map(|(a, b)| a * b).sum::<f64>();
                }
                gram[i][i] += LAMBDA;
            }
            let y = gaussian_elimination(&mut gram, &mut rhs).unwrap_or_else(|| vec![0.0_f64; r]);
            let max_c = y.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            let support: Vec<usize> = (0..r).filter(|&jj| y[jj].abs() > 1.0e-6 * max_c).collect();
            // Self-verification: the reported group must be a REAL
            // relation — reconstruct the dependent row from the support
            // and reject the group when the reconstruction residual is
            // not at roundoff (an ill-conditioned basis must never
            // surface a fake dependency group).
            let residual: f64 = (0..n_free_c)
                .map(|rr| {
                    let recon: f64 = support.iter().map(|&jj| y[jj] * rows[jj][rr]).sum();
                    (recon - dep_row[rr]).abs()
                })
                .fold(0.0_f64, f64::max);
            let scale = dep_row.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            if !(residual <= 1.0e-9 * scale) {
                continue;
            }
            let mut group: Vec<usize> = support
                .iter()
                .map(|&jj| row_constraint[independent[jj]])
                .chain(std::iter::once(row_constraint[dep]))
                .collect();
            group.sort_unstable();
            group.dedup();
            if group.len() > 1 && !d.dependency_groups.contains(&group) {
                d.dependency_groups.push(group);
            }
        }

        // Least-squares step of the independent rows (minimizes |Jₛ z + fₛ|).
        let j_sub: Vec<Vec<f64>> = independent.iter().map(|&i| j[i].clone()).collect();
        let f_sub: Vec<f64> = independent.iter().map(|&i| f_c[i]).collect();
        let z = solve_normal_equations(&j_sub, &f_sub, n_free_c);

        let scale = f_c.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        for &r in &dependent {
            let drive: f64 = j[r].iter().zip(z.iter()).map(|(ji, zi)| ji * zi).sum();
            let e = (f_c[r] + drive).abs();
            let global = row_constraint[r];
            if e <= CONSIST_TOL * scale {
                if !d.redundant.contains(&global) {
                    d.redundant.push(global);
                }
            } else if e <= PARTIAL_TOL * scale {
                if !d.partially_redundant.contains(&global) {
                    d.partially_redundant.push(global);
                }
            } else if !d.conflicting.contains(&global) {
                d.conflicting.push(global);
            }
        }
    }
    d.dof = free_total as i64 - total_rank;
    d.rank = total_rank.max(0) as usize;
    d
}

/// Rank of the constraint Jacobian over the free parameters at `params`.
///
/// Thin convenience wrapper over [`diagnose`]: same per-component
/// rank-revealing QR on the ANALYTIC Jacobian, same deterministic
/// column-pivoting tie-breaks, same linearization point (the current
/// parameters). Callers that also need the redundant/conflicting
/// classification should call [`diagnose`] directly.
pub fn rank(
    params: &[f64],
    fixed: &[bool],
    entities: &[Entity],
    constraints: &[Constraint],
) -> usize {
    diagnose(params, fixed, entities, constraints).rank
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::{EntityId, PointRef};
    use crate::sketch::Sketch;

    // ── direct QR kernel test ────────────────────────────────────────────

    #[test]
    fn qr_reveals_rank_and_pivots() {
        // Columns: col1 has the largest norm and pivots first; then cols 0
        // and 2 tie at 1/√2 (the earlier position wins), col 3 is zero.
        // col2 = col1 − col0, so the rank-revealing set is {1, 0}.
        let mut a = vec![
            vec![1.0, 1.0, 0.0, 0.0],
            vec![0.0, 1.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0],
        ];
        let (pivots, rank) = householder_qr_column_pivoting(&mut a, 1.0e-10);
        assert_eq!(rank, 2);
        assert_eq!(pivots[..2], [1, 0]);
        assert_eq!(pivots[2..], [2, 3]);
    }

    /// v2 dependency groups: the 3-4-5 directional family's whole
    /// mutually-dependent set — DistanceX + DistanceY + Length — lands in
    /// ONE group, not a pivot-chosen single member (the member the flat
    /// `redundant` set names flips with pivot tie-breaking; the group
    /// does not).
    #[test]
    fn dependency_groups_name_the_whole_relation() {
        let mut sk = Sketch::new();
        let line = sk.add_line(0.0, 0.0, 3.0, 4.0);
        sk.add_constraint(crate::constraint::Constraint::DistanceX {
            p1: PointRef::LineEnd(line),
            p2: PointRef::LineStart(line),
            distance: 3.0,
        });
        sk.add_constraint(crate::constraint::Constraint::DistanceY {
            p1: PointRef::LineEnd(line),
            p2: PointRef::LineStart(line),
            distance: 4.0,
        });
        sk.add_constraint(crate::constraint::Constraint::LineLength { line, length: 5.0 });
        // Pin the start point through the solver mask.
        sk.fix_param(0);
        sk.fix_param(1);

        let d = crate::diagnosis::diagnose(
            &sk.params_snapshot(),
            &sk.fixed,
            &sk.entities,
            sk.constraints(),
        );
        assert_eq!(d.dof, 0, "fully constrained: {:?}", d);
        assert_eq!(
            d.dependency_groups,
            [vec![0, 1, 2]],
            "the whole DX/DY/Length relation is one group: {:?}",
            d.dependency_groups
        );
        assert!(!d.redundant.is_empty(), "a dependent row is classified");
    }

    /// Rank API (M3.2): the rank is free parameters − dof, and counting a
    /// candidate row at the current pose shows the increment directly. A
    /// row whose support is already determined (masked entity) adds none.
    #[test]
    fn rank_reports_free_minus_dof_and_candidate_increments() {
        let mut sk = Sketch::new();
        let line = sk.add_line(0.0, 0.0, 3.0, 4.0);

        let d0 = sk.diagnose();
        assert_eq!((d0.rank, d0.dof), (0, 4), "{d0:?}");
        assert_eq!(d0.rank as i64, 4 - d0.dof, "rank = free − dof");

        sk.add_constraint(Constraint::LineLength { line, length: 5.0 });
        let d1 = sk.diagnose();
        assert_eq!((d1.rank, d1.dof), (1, 3), "{d1:?}");

        sk.add_constraint(Constraint::DistanceX {
            p1: PointRef::LineEnd(line),
            p2: PointRef::LineStart(line),
            distance: 3.0,
        });
        let d2 = sk.diagnose();
        assert_eq!((d2.rank, d2.dof), (2, 2), "{d2:?}");

        // The y separation is linearly dependent on Length + DistanceX at
        // this pose (the 3-4-5 family spans two directions only).
        sk.add_constraint(Constraint::DistanceY {
            p1: PointRef::LineEnd(line),
            p2: PointRef::LineStart(line),
            distance: 4.0,
        });
        let d3 = sk.diagnose();
        assert_eq!(d3.rank, 2, "dependent row adds no rank: {d3:?}");
        assert!(!d3.redundant.is_empty(), "{d3:?}");

        // A fully masked entity gives a zero free-parameter row space.
        let locked = sk.add_line(10.0, 0.0, 14.0, 0.0);
        sk.fix_entity(locked);
        let before = sk.diagnose().rank;
        sk.add_constraint(Constraint::LineLength {
            line: locked,
            length: 4.0,
        });
        let after = sk.diagnose().rank;
        assert_eq!(after, before, "a masked row adds no rank");
    }

    // 10 conflicting / 10 redundant / 10 under-constrained, each visible in
    // the initial linearization (the honest M1 scope — later-manifesting
    // nonlinear conflicts are the M5 corpus's looser targets).

    type Build = Box<dyn FnOnce(&mut Sketch)>;

    fn point(sk: &mut Sketch, x: f64, y: f64) -> EntityId {
        sk.add_point(x, y)
    }

    fn fixed(sk: &mut Sketch, id: EntityId, x: f64, y: f64) {
        sk.add_constraint(Constraint::fix_point(PointRef::Point(id), x, y));
    }

    #[test]
    fn seed_fault_set_classification() {
        let conflicting: Vec<Build> = vec![
            // C1: one point demanded at two places.
            Box::new(|sk: &mut Sketch| {
                let p = point(sk, 0.0, 0.0);
                fixed(sk, p, 0.0, 0.0);
                fixed(sk, p, 1.0, 1.0);
            }),
            // C2: radius and diameter disagreeing on the same parameter.
            Box::new(|sk: &mut Sketch| {
                let c = sk.add_circle(0.0, 0.0, 10.0);
                sk.add_constraint(Constraint::Radius {
                    circle: c,
                    radius: 3.0,
                });
                sk.add_constraint(Constraint::Diameter {
                    circle: c,
                    diameter: 50.0,
                });
            }),
            // C3: fixed points further apart than the demanded distance.
            Box::new(|sk: &mut Sketch| {
                let p1 = point(sk, 0.0, 0.0);
                let p2 = point(sk, 1.0, 1.0);
                fixed(sk, p1, 0.0, 0.0);
                fixed(sk, p2, 1.0, 1.0);
                sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
            }),
            // C4: directed x-distance against two fixed points.
            Box::new(|sk: &mut Sketch| {
                let p1 = point(sk, 0.0, 0.0);
                let p2 = point(sk, 1.0, 1.0);
                fixed(sk, p1, 0.0, 0.0);
                fixed(sk, p2, 1.0, 1.0);
                sk.add_constraint(Constraint::distance_x(p2, p1, 5.0));
            }),
            // C5: a line whose endpoints are pinned at a different length.
            Box::new(|sk: &mut Sketch| {
                let l = sk.add_line(0.0, 0.0, 3.0, 0.0);
                sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
                sk.add_constraint(Constraint::fix_point(PointRef::LineEnd(l), 3.0, 0.0));
                sk.add_constraint(Constraint::LineLength {
                    line: l,
                    length: 7.0,
                });
            }),
            // C6: concentric circles whose centers are pinned apart.
            Box::new(|sk: &mut Sketch| {
                let c1 = sk.add_circle(0.0, 0.0, 1.0);
                let c2 = sk.add_circle(4.0, 0.0, 1.0);
                sk.add_constraint(Constraint::fix_point(PointRef::Center(c1), 0.0, 0.0));
                sk.add_constraint(Constraint::fix_point(PointRef::Center(c2), 4.0, 0.0));
                sk.add_constraint(Constraint::Concentric(c1, c2));
            }),
            // C7: midpoint of a pinned line demanded at a different point.
            Box::new(|sk: &mut Sketch| {
                let l = sk.add_line(0.0, 0.0, 4.0, 0.0);
                sk.fix_entity(l);
                let p = point(sk, 9.0, 9.0);
                fixed(sk, p, 9.0, 9.0);
                sk.add_constraint(Constraint::Midpoint {
                    point: PointRef::Point(p),
                    line: l,
                });
            }),
            // C8: equal lengths on two pinned lines of different lengths.
            Box::new(|sk: &mut Sketch| {
                let l1 = sk.add_line(0.0, 0.0, 1.0, 0.0);
                let l2 = sk.add_line(5.0, 5.0, 9.0, 5.0);
                sk.fix_entity(l1);
                sk.fix_entity(l2);
                sk.add_constraint(Constraint::EqualLength(l1, l2));
            }),
            // C9: parallel and perpendicular demanded of the same pair —
            // generic (non-axis-aligned) position so both demands load the
            // single free parameter (l2's end x).
            Box::new(|sk: &mut Sketch| {
                let l1 = sk.add_line(0.0, 0.0, 3.0, 4.0);
                let l2 = sk.add_line(0.0, 5.0, 1.0, 8.0);
                sk.fix_entity(l1);
                sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 0.0, 5.0));
                let y_param = sk.entities[l2].param_start + 3; // l2's end y
                sk.fix_param(y_param);
                sk.add_constraint(Constraint::Parallel(l1, l2));
                sk.add_constraint(Constraint::Perpendicular(l1, l2));
            }),
            // C10: distance whose components disagree (a²+b² ≠ d²).
            Box::new(|sk: &mut Sketch| {
                let p1 = point(sk, 0.0, 0.0);
                let p2 = point(sk, 3.0, 4.0);
                fixed(sk, p1, 0.0, 0.0);
                sk.add_constraint(Constraint::point_distance(p1, p2, 10.0));
                sk.add_constraint(Constraint::distance_x(p2, p1, 3.0));
                sk.add_constraint(Constraint::distance_y(p2, p1, 4.0));
            }),
        ];

        let redundant: Vec<Build> = vec![
            // R1: duplicate fixed.
            Box::new(|sk: &mut Sketch| {
                let p = point(sk, 0.0, 0.0);
                fixed(sk, p, 3.0, 4.0);
                fixed(sk, p, 3.0, 4.0);
            }),
            // R2: duplicate coincident.
            Box::new(|sk: &mut Sketch| {
                let p1 = point(sk, 0.0, 0.0);
                let p2 = point(sk, 5.0, 5.0);
                fixed(sk, p1, 1.0, 1.0);
                sk.add_constraint(Constraint::coincident(p1, p2));
                sk.add_constraint(Constraint::coincident(p1, p2));
            }),
            // R3: distance fully implied by consistent components.
            Box::new(|sk: &mut Sketch| {
                let p1 = point(sk, 0.0, 0.0);
                let p2 = point(sk, 3.0, 4.0);
                fixed(sk, p1, 0.0, 0.0);
                sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
                sk.add_constraint(Constraint::distance_x(p2, p1, 3.0));
                sk.add_constraint(Constraint::distance_y(p2, p1, 4.0));
            }),
            // R4: consistent radius + diameter (same demand).
            Box::new(|sk: &mut Sketch| {
                let c = sk.add_circle(1.0, 1.0, 5.0);
                sk.add_constraint(Constraint::Radius {
                    circle: c,
                    radius: 5.0,
                });
                sk.add_constraint(Constraint::Diameter {
                    circle: c,
                    diameter: 10.0,
                });
            }),
            // R5: point already at the pinned line's midpoint.
            Box::new(|sk: &mut Sketch| {
                let l = sk.add_line(0.0, 0.0, 4.0, 0.0);
                sk.fix_entity(l);
                let p = point(sk, 2.0, 0.0);
                fixed(sk, p, 2.0, 0.0);
                sk.add_constraint(Constraint::Midpoint {
                    point: PointRef::Point(p),
                    line: l,
                });
            }),
            // R6: duplicate horizontal.
            Box::new(|sk: &mut Sketch| {
                let l = sk.add_line(0.0, 0.0, 5.0, 0.0);
                sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
                sk.add_constraint(Constraint::Horizontal(l));
                sk.add_constraint(Constraint::Horizontal(l));
            }),
            // R7: concentric circles whose centers are both already pinned
            // at the same spot (c1 by the fixed mask, c2 by constraints).
            Box::new(|sk: &mut Sketch| {
                let c1 = sk.add_circle(0.0, 0.0, 1.0);
                let c2 = sk.add_circle(0.0, 0.0, 2.0);
                sk.fix_entity(c1);
                sk.add_constraint(Constraint::fix_point(PointRef::Center(c2), 0.0, 0.0));
                sk.add_constraint(Constraint::Concentric(c1, c2));
            }),
            // R8: both points already pinned symmetric about a pinned line —
            // the Symmetric row is fully implied.
            Box::new(|sk: &mut Sketch| {
                let l = sk.add_line(0.0, 0.0, 0.0, 8.0);
                sk.fix_entity(l);
                let p1 = point(sk, 3.0, 2.0);
                let p2 = point(sk, -3.0, 2.0);
                fixed(sk, p1, 3.0, 2.0);
                fixed(sk, p2, -3.0, 2.0);
                sk.add_constraint(Constraint::Symmetric {
                    p1: PointRef::Point(p1),
                    p2: PointRef::Point(p2),
                    line: l,
                });
            }),
            // R9: directed x-distance duplicating a fixed displacement.
            Box::new(|sk: &mut Sketch| {
                let p1 = point(sk, 0.0, 0.0);
                let p2 = point(sk, 4.0, 7.0);
                fixed(sk, p1, 0.0, 0.0);
                fixed(sk, p2, 4.0, 7.0);
                sk.add_constraint(Constraint::distance_x(p2, p1, 4.0));
            }),
            // R10: duplicate point-on-line.
            Box::new(|sk: &mut Sketch| {
                let l = sk.add_line(0.0, 0.0, 6.0, 0.0);
                sk.fix_entity(l);
                let p = point(sk, 2.0, 0.0);
                sk.add_constraint(Constraint::PointOnLine {
                    point: PointRef::Point(p),
                    line: l,
                });
                sk.add_constraint(Constraint::PointOnLine {
                    point: PointRef::Point(p),
                    line: l,
                });
            }),
        ];

        let under: Vec<(&str, Build)> = vec![
            (
                "no constraints",
                Box::new(|sk: &mut Sketch| {
                    point(sk, 1.0, 2.0);
                }),
            ),
            (
                "line with fixed start",
                Box::new(|sk: &mut Sketch| {
                    let l = sk.add_line(0.0, 0.0, 3.0, 4.0);
                    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
                }),
            ),
            (
                "circle with radius only",
                Box::new(|sk: &mut Sketch| {
                    let c = sk.add_circle(1.0, 1.0, 2.0);
                    sk.add_constraint(Constraint::Radius {
                        circle: c,
                        radius: 2.0,
                    });
                }),
            ),
            (
                "coincident pair",
                Box::new(|sk: &mut Sketch| {
                    let p1 = point(sk, 0.0, 0.0);
                    let p2 = point(sk, 5.0, 5.0);
                    sk.add_constraint(Constraint::coincident(p1, p2));
                }),
            ),
            (
                "horizontal line",
                Box::new(|sk: &mut Sketch| {
                    let l = sk.add_line(0.0, 0.0, 5.0, 1.0);
                    sk.add_constraint(Constraint::Horizontal(l));
                }),
            ),
            (
                "single distance",
                Box::new(|sk: &mut Sketch| {
                    let p1 = point(sk, 0.0, 0.0);
                    let p2 = point(sk, 3.0, 4.0);
                    sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
                }),
            ),
            (
                "midpoint of free line",
                Box::new(|sk: &mut Sketch| {
                    let l = sk.add_line(0.0, 0.0, 4.0, 0.0);
                    let p = point(sk, 2.0, 0.0);
                    sk.add_constraint(Constraint::Midpoint {
                        point: PointRef::Point(p),
                        line: l,
                    });
                }),
            ),
            (
                "equal lengths of free lines",
                Box::new(|sk: &mut Sketch| {
                    let l1 = sk.add_line(0.0, 0.0, 1.0, 0.0);
                    let l2 = sk.add_line(5.0, 5.0, 8.0, 5.0);
                    sk.add_constraint(Constraint::EqualLength(l1, l2));
                }),
            ),
            (
                "point on free circle",
                Box::new(|sk: &mut Sketch| {
                    let c = sk.add_circle(0.0, 0.0, 3.0);
                    let p = point(sk, 3.0, 0.0);
                    sk.add_constraint(Constraint::PointOnCircle {
                        point: PointRef::Point(p),
                        circle: c,
                    });
                }),
            ),
            (
                "free symmetric pair",
                Box::new(|sk: &mut Sketch| {
                    let l = sk.add_line(0.0, 0.0, 0.0, 8.0);
                    let p1 = point(sk, 3.0, 2.0);
                    let p2 = point(sk, -3.0, 2.0);
                    sk.add_constraint(Constraint::Symmetric {
                        p1: PointRef::Point(p1),
                        p2: PointRef::Point(p2),
                        line: l,
                    });
                }),
            ),
        ];

        let conflict_names = [
            "C1 point-demanded-twice",
            "C2 radius-vs-diameter",
            "C3 distance-beyond-fixed-span",
            "C4 distance_x-vs-fixed",
            "C5 line-length-vs-pinned-ends",
            "C6 concentric-vs-pinned-centers",
            "C7 midpoint-vs-pinned",
            "C8 equal-length-vs-pinned",
            "C9 parallel+perpendicular",
            "C10 distance-vs-components",
        ];
        let mut conflicts_flagged = 0usize;
        for (name, build) in conflict_names.into_iter().zip(conflicting) {
            let mut sk = Sketch::new();
            build(&mut sk);
            let d = sk.diagnose();
            assert!(
                !d.conflicting.is_empty(),
                "{name}: expected a conflict, got {d:?}"
            );
            conflicts_flagged += 1;
        }
        assert_eq!(conflicts_flagged, 10, "conflict recall target is 1.0");

        let redundant_names = [
            "R1 duplicate-fixed",
            "R2 duplicate-coincident",
            "R3 distance-implied",
            "R4 consistent-radius-diameter",
            "R5 point-at-pinned-midpoint",
            "R6 duplicate-horizontal",
            "R7 concentric-on-pinned",
            "R8 symmetric-satisfied",
            "R9 distance_x-duplicating",
            "R10 duplicate-point-on-line",
        ];
        for (name, build) in redundant_names.into_iter().zip(redundant) {
            let mut sk = Sketch::new();
            build(&mut sk);
            let d = sk.diagnose();
            assert!(
                !d.redundant.is_empty(),
                "{name}: expected a redundant row, got {d:?}"
            );
            assert!(
                d.conflicting.is_empty(),
                "redundant case misclassified as conflicting: {d:?}"
            );
        }

        for (name, build) in under {
            let mut sk = Sketch::new();
            build(&mut sk);
            let d = sk.diagnose();
            assert!(
                d.dof > 0 && d.conflicting.is_empty(),
                "{name}: expected under-constrained, got {d:?}"
            );
        }
    }
}

#[cfg(test)]
mod degenerate_and_bench {
    use crate::constraint::Constraint;
    use crate::entity::PointRef;
    use crate::sketch::Sketch;
    use std::time::Instant;

    /// Degenerate family (M2.3 acceptance): rank judgment must be jitter-
    /// free at exact degeneracies — tangency, collinearity, zero margins —
    /// and stay stable under small perturbations of the same configuration.
    #[test]
    fn degenerate_families_rank_stable() {
        // Exactly tangent circle-line, point exactly on the line, duplicate
        // coincident rows: the dependent structure must be reported
        // identically at the exact point and at nearby perturbations.
        let build = |jitter: f64| {
            let mut sk = Sketch::new();
            let l = sk.add_line(0.0, 0.0, 6.0, 0.0);
            sk.fix_entity(l);
            let c = sk.add_circle(3.0, 1.0 + jitter, 1.0);
            sk.add_constraint(Constraint::Tangent { circle: c, line: l });
            let p = sk.add_point(2.0 + jitter, 0.0);
            sk.add_constraint(Constraint::PointOnLine {
                point: PointRef::Point(p),
                line: l,
            });
            let q = sk.add_point(4.0, jitter);
            sk.add_constraint(Constraint::PointOnLine {
                point: PointRef::Point(q),
                line: l,
            });
            sk.add_constraint(Constraint::PointOnLine {
                point: PointRef::Point(q),
                line: l,
            });
            sk.diagnose()
        };
        let exact = build(0.0);
        for eps in [1e-10, 1e-8, 1e-6] {
            let perturbed = build(eps);
            assert_eq!(
                exact.dof, perturbed.dof,
                "dof jitters at eps={eps}: {exact:?} vs {perturbed:?}"
            );
            assert_eq!(
                exact.conflicting, perturbed.conflicting,
                "conflicting set jitters at eps={eps}"
            );
            assert_eq!(
                exact.redundant.len(),
                perturbed.redundant.len(),
                "redundant count jitters at eps={eps}"
            );
        }
        // The duplicate PointOnLine must be seen as redundant, never as a
        // conflict.
        assert!(exact.conflicting.is_empty(), "{exact:?}");
        assert!(!exact.redundant.is_empty(), "{exact:?}");
    }

    /// M2.2 acceptance measurement (ignored by default; run explicitly):
    /// `cargo test -p rcad-gcs --lib -- --ignored --nocapture jacobian_bench`
    /// A ~thousand-parameter sketch: analytic assembly vs finite-difference
    /// build, single pass each.
    #[test]
    #[ignore]
    fn jacobian_bench_thousand_params() {
        // 500 points chained by distance constraints: 1000 params, 499
        // constraints (+1 fixed point).
        let n = 500usize;
        let mut sk = Sketch::new();
        for i in 0..n {
            sk.add_point(i as f64, 0.0);
        }
        sk.add_constraint(Constraint::fix_point(PointRef::Point(0), 0.0, 0.0));
        for i in 0..n - 1 {
            sk.add_constraint(Constraint::point_distance(
                PointRef::Point(i),
                PointRef::Point(i + 1),
                1.0,
            ));
        }

        let t0 = Instant::now();
        let d = sk.diagnose(); // one analytic Jacobian + QR
        let analytic_us = t0.elapsed().as_secs_f64() * 1e6;
        println!(
            "analytic (diagnose, incl. QR): {analytic_us:.0} µs, dof {}",
            d.dof
        );

        // Same-size finite-difference build (the pre-M2 path): one pass =
        // 2 × n_free full residual evaluations. Both builders are timed
        // ALONE — the QR / normal-equations cost is algorithmic (dense at
        // this size; M3's decomposition targets it), not Jacobian-build
        // cost.
        use crate::entity::{Entity, EntityKind};
        use crate::numerics::numerical_jacobian_into;
        use crate::solver::analytic_jacobian;
        let entities: Vec<Entity> = (0..n)
            .map(|i| Entity::new(EntityKind::Point, i * 2))
            .collect();
        let mut constraints = vec![Constraint::fix_point(PointRef::Point(0), 0.0, 0.0)];
        for i in 0..n - 1 {
            constraints.push(Constraint::point_distance(
                PointRef::Point(i),
                PointRef::Point(i + 1),
                1.0,
            ));
        }
        let n_eq: usize = constraints.iter().map(|c| c.equation_count()).sum();
        let free: Vec<usize> = (0..n * 2).collect();
        let col_of: std::collections::HashMap<usize, usize> =
            free.iter().enumerate().map(|(c, &p)| (p, c)).collect();
        let mut ja = vec![vec![0.0_f64; free.len()]; n_eq];
        let params: Vec<f64> = (0..n).flat_map(|i| [i as f64, 0.0]).collect();
        let t0 = Instant::now();
        analytic_jacobian(&params, &col_of, &entities, &constraints, &mut ja);
        let analytic_us = t0.elapsed().as_secs_f64() * 1e6;
        println!("analytic assembly:                {analytic_us:.0} µs");

        let mut eval = |p: &[f64], out: &mut [f64]| {
            let v = crate::solver::eval_residuals(p, &entities, &constraints, n_eq);
            out.copy_from_slice(&v);
        };
        let mut jfd = vec![vec![0.0_f64; free.len()]; n_eq];
        let t1 = Instant::now();
        numerical_jacobian_into(&params, &free, n_eq, &mut eval, &mut jfd);
        let fd_us = t1.elapsed().as_secs_f64() * 1e6;
        println!("finite differences (same size):   {fd_us:.0} µs");
        println!("build speedup: {:.1}×", fd_us / analytic_us.max(1e-9));
    }
}

#[cfg(test)]
mod decomposition_acceptance {
    use crate::constraint::Constraint;
    use crate::entity::PointRef;
    use crate::sketch::Sketch;

    /// Two self-contained groups (anchored driven line + anchored driven
    /// circle), built twice: solved ALL at once vs solved per component.
    /// The per-component path must land on the same parameters — the
    /// two-pass construction is what makes incremental re-solve
    /// result-identical by construction.
    fn build() -> Sketch {
        let mut sk = Sketch::new();
        // Group A: fixed point, coincident line start, driven length.
        // Entities: point id 0 (params 0-1), line id 1 (params 2-5).
        let pa = sk.add_point(0.0, 0.0);
        let la = sk.add_line(0.0, 0.0, 5.0, 0.0);
        sk.add_constraint(Constraint::fix_point(PointRef::Point(pa), 0.0, 0.0));
        sk.add_constraint(Constraint::coincident(PointRef::LineStart(la), pa));
        sk.add_constraint(Constraint::Horizontal(la));
        sk.add_constraint(Constraint::LineLength {
            line: la,
            length: 5.0,
        });
        // Group B: pinned circle center + driven radius.
        // Entity: circle id 2 (params 6-8). Constraints 4-5.
        let cb = sk.add_circle(9.0, 9.0, 2.0);
        sk.add_constraint(Constraint::fix_point(PointRef::Center(cb), 9.0, 9.0));
        sk.add_constraint(Constraint::Radius {
            circle: cb,
            radius: 2.0,
        });
        sk
    }

    #[test]
    fn per_component_solve_matches_all_at_once() {
        let mut all = build();
        let r_all = all.solve();
        assert!(r_all.converged, "residual {}", r_all.residual);

        let mut piecewise = build();
        let comps = piecewise.components();
        // Three components: group A (constraints 0-3), the circle's center pin
        // (4) and its radius (5) share no parameters — genuinely independent.
        assert_eq!(comps.len(), 3, "{comps:?}");
        for c in 0..comps.len() {
            let r = crate::solver::solve_components(
                &mut piecewise.params,
                &piecewise.fixed,
                &piecewise.entities,
                &piecewise.constraints,
                &comps,
                Some(&[c]),
            );
            assert!(r.converged, "component {c}: residual {}", r.residual);
        }
        for (a, b) in all.params.iter().zip(piecewise.params.iter()) {
            assert!(
                (a - b).abs() < 1e-12,
                "per-component params diverge: {a} vs {b}"
            );
        }
    }

    /// Localized solving (M3.2): re-solving ONLY group A's component
    /// leaves group B's parameters bitwise untouched, even after group A
    /// was drifted away from its solution and its driving length edited.
    #[test]
    fn untouched_component_params_bitwise_stable() {
        let mut sk = build();
        assert_eq!(sk.components().len(), 3);
        assert!(sk.solve().converged);
        let before = sk.params.clone();

        // Drift group A's line and shrink its driving length (constraint
        // index 3, line entity id 1).
        sk.params[2] += 0.8;
        sk.params[4] += 1.5;
        sk.constraints[3] = Constraint::LineLength {
            line: 1,
            length: 3.0,
        };

        let comps = sk.components();
        let a_component = comps
            .iter()
            .position(|g| g.contains(&3))
            .expect("edited constraint's component");
        let r = crate::solver::solve_components(
            &mut sk.params,
            &sk.fixed,
            &sk.entities,
            &sk.constraints,
            &comps,
            Some(&[a_component]),
        );
        assert!(r.converged, "residual {}", r.residual);

        // Group B's params (circle id 2 → params 6-8) are bitwise unchanged.
        let circle_start = sk.entities[2].param_start;
        for (i, p) in sk.params.iter().enumerate().skip(circle_start) {
            assert_eq!(*p, before[i], "group B param {i} moved");
        }
        // And group A actually re-solved to its new length.
        let dx = sk.params[4] - sk.params[2];
        let dy = sk.params[5] - sk.params[3];
        assert!(
            ((dx * dx + dy * dy).sqrt() - 3.0).abs() < 1e-9,
            "group A re-solved to length 3, got {}",
            (dx * dx + dy * dy).sqrt()
        );
    }

    /// A conflict inside one component is reported with that component's
    /// global indices only; the healthy component contributes no rows and
    /// the DOF is the sum over components (never lost, never negative).
    #[test]
    fn conflict_isolated_to_its_component() {
        let mut sk = build();
        // Sabotage group B: a second, disagreeing radius (constraint 6).
        sk.add_constraint(Constraint::Radius {
            circle: 2,
            radius: 7.0,
        });
        let d = sk.diagnose();
        assert_eq!(d.conflicting, vec![6], "{d:?}");
        assert!(d.redundant.is_empty(), "{d:?}");
        assert!(d.dof >= 0, "{d:?}");
    }
}
