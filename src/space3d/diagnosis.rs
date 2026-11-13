//! Rank-revealing diagnosis of a 3D constraint system.
//!
//! Mirrors the 2D [`crate::diagnosis`] module over the analytic 3D
//! Jacobian: the row space is analyzed by a column-pivoted Householder QR of
//! the transposed Jacobian, the pivot columns are the linearly independent
//! constraint rows, and the complement is classified by whether its residual
//! agrees with the least-squares step of the independent rows:
//!
//! - **redundant** — dependent and consistent (removable without changing
//!   the solution);
//! - **conflicting** — dependent and inconsistent (no solution satisfies the
//!   whole system at this linearization);
//! - **partially_redundant** — the gray band between the two.
//!
//! The analysis runs PER COMPONENT of the 3D constraint graph (dependencies
//! can only occur inside a component), and reported indices are positions in
//! the `constraints` slice.
//!
//! # Plane gauge
//!
//! A plane is stored as `(nx, ny, nz, d)`; uniform scaling of the pair
//! leaves every residual invariant, so only three of the four coordinates
//! are physical. The plane's free columns are therefore replaced by the two
//! tangent directions of its unit normal plus the offset direction (the
//! radial scaling direction is dropped), and the DOF count subtracts one
//! gauge dimension per fully free plane. Without this, every free plane
//! would contribute a phantom DOF and a fully constrained plane could show a
//! spurious dependency.

use crate::numerics::{
    gaussian_elimination, householder_qr_column_pivoting, solve_normal_equations, LAMBDA,
};
use crate::space3d::constraint::SpaceConstraint;
use crate::space3d::entity::{SpaceEntity, SpaceEntityKind};
use crate::space3d::solver::{analytic_jacobian, eval_residuals, SOLVER_NORM_TOL};
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

/// Structured diagnosis of a 3D constraint system.
#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SpaceDiagnosis {
    /// Jacobian rank: number of linearly independent constraint rows over
    /// all components (redundant rows consume none).
    pub rank: usize,
    /// True degrees of freedom: effective free parameters − Jacobian rank.
    pub dof: i64,
    /// Constraints that are dependent but consistent — removable without
    /// changing the solution.
    pub redundant: Vec<usize>,
    /// Constraints whose demands are inconsistent with the rest — no
    /// solution exists at this linearization.
    pub conflicting: Vec<usize>,
    /// Dependent rows in the gray band (weakly inconsistent; possibly
    /// genuine nonlinearity). Informational.
    pub partially_redundant: Vec<usize>,
    /// Mutually-dependent constraint groups: each group is the full support
    /// of one linear dependency among rows, extracted from the pivoted QR's
    /// R factor (flip-free: the whole relation, not a pivot-chosen member).
    pub dependency_groups: Vec<Vec<usize>>,
}

/// Diagnose a 3D constraint system at the current parameters.
///
/// Constraint indices are the positions in `constraints` (equations stack in
/// order, 1–3 rows per constraint). Reported indices are global.
pub fn diagnose(
    params: &[f64],
    fixed: &[bool],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
) -> SpaceDiagnosis {
    let free_total = params.iter().zip(fixed).filter(|(_, &f)| !f).count();
    let n_eq: usize = constraints.iter().map(|c| c.equation_count()).sum();

    // One gauge dimension per fully free plane with a usable normal.
    let gauge_reduction: usize = entities
        .iter()
        .filter(|e| e.kind == SpaceEntityKind::Plane)
        .filter(|e| {
            let idx = [e.param(0), e.param(1), e.param(2), e.param(3)];
            idx.iter().all(|&i| !fixed[i])
                && plane_tangents(
                    params[idx[0]],
                    params[idx[1]],
                    params[idx[2]],
                )
                .is_some()
        })
        .count();

    if n_eq == 0 {
        return SpaceDiagnosis {
            dof: free_total as i64 - gauge_reduction as i64,
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

    let components =
        crate::space3d::graph::constraint_components(params, entities, constraints);
    let mut d = SpaceDiagnosis::default();
    let mut total_rank = 0i64;

    for group in &components {
        let sub: Vec<&SpaceConstraint> = group.iter().map(|&i| &constraints[i]).collect();
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

        // Map component-local rows back to global constraint indices.
        let mut row_constraint: Vec<usize> = Vec::with_capacity(n_eq_c);
        for &gi in group {
            for _ in 0..constraints[gi].equation_count() {
                row_constraint.push(gi);
            }
        }

        if touched.is_empty() {
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

        // Raw Jacobian over the touched local columns.
        let col_of: HashMap<usize, usize> =
            touched.iter().enumerate().map(|(c, &p)| (p, c)).collect();
        let mut j_raw = vec![vec![0.0_f64; touched.len()]; n_eq_c];
        analytic_jacobian(params, &col_of, entities, &sub, &mut j_raw);

        // Effective free directions: plane gauge is removed here.
        let columns = effective_columns(params, fixed, entities, &touched);
        let n_eff = columns.len();
        if n_eff == 0 {
            let scale = f_c.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            for (r, fi) in f_c.iter().enumerate() {
                if fi.abs() > CONSIST_TOL * scale {
                    d.conflicting.push(row_constraint[r]);
                }
            }
            continue;
        }
        let mut j = vec![vec![0.0_f64; n_eff]; n_eq_c];
        for (row, jr) in j_raw.iter().enumerate() {
            for (c, col) in columns.iter().enumerate() {
                let mut v = 0.0_f64;
                for &(lc, coeff) in col {
                    v += jr[lc] * coeff;
                }
                j[row][c] = v;
            }
        }

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

        // A = J^T: column i of A is component row i (equilibrated).
        let mut a = vec![vec![0.0_f64; n_eq_c]; n_eff];
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
        let dependent: Vec<usize> =
            (0..n_eq_c).filter(|&i| !independent.contains(&i)).collect();
        if dependent.is_empty() || rank == 0 {
            continue;
        }

        // Null-space dependency groups: the dependency lives in ROW space —
        // find coefficients y with SUM y_j*row_j = dep_row via the row-space
        // Gram system (G + lambda*I) y = -(rows*dep), G[i][j] = row_i.row_j.
        // The coefficient support is the whole mutually-dependent row group,
        // invariant to pivot tie-breaking.
        let rows: Vec<Vec<f64>> = independent
            .iter()
            .map(|&col| (0..n_eff).map(|rr| a_original[rr][col]).collect())
            .collect();
        let r = rows.len();
        for &dep in &dependent {
            let dep_row: Vec<f64> = (0..n_eff).map(|rr| a_original[rr][dep]).collect();
            let mut gram = vec![vec![0.0_f64; r]; r];
            let mut rhs = vec![0.0_f64; r];
            for (i, ri) in rows.iter().enumerate() {
                rhs[i] = ri.iter().zip(dep_row.iter()).map(|(a, b)| a * b).sum::<f64>();
                for (j, rj) in rows.iter().enumerate() {
                    gram[i][j] = ri.iter().zip(rj.iter()).map(|(a, b)| a * b).sum::<f64>();
                }
                gram[i][i] += LAMBDA;
            }
            let y = gaussian_elimination(&mut gram, &mut rhs)
                .unwrap_or_else(|| vec![0.0_f64; r]);
            let max_c = y.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            let support: Vec<usize> = (0..r)
                .filter(|&jj| y[jj].abs() > 1.0e-6 * max_c)
                .collect();
            // Self-verification: the reported group must be a REAL relation —
            // reconstruct the dependent row from the support and reject the
            // group when the reconstruction residual is not at roundoff.
            let residual: f64 = (0..n_eff)
                .map(|rr| {
                    let recon: f64 =
                        support.iter().map(|&jj| y[jj] * rows[jj][rr]).sum();
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

        // Least-squares step of the independent rows (minimizes |J_s z + f_s|).
        let j_sub: Vec<Vec<f64>> = independent.iter().map(|&i| j[i].clone()).collect();
        let f_sub: Vec<f64> = independent.iter().map(|&i| f_c[i]).collect();
        let z = solve_normal_equations(&j_sub, &f_sub, n_eff);

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
    d.rank = total_rank as usize;
    d.dof = free_total as i64 - gauge_reduction as i64 - total_rank;
    d
}

/// Effective free directions over the component-local `touched` columns.
///
/// Fully free planes contribute two normal-tangent directions and the offset
/// direction instead of four raw coordinate columns (the radial scaling
/// direction is gauge). Everything else contributes its unit column. The
/// returned vectors are sparse combinations of local column indices.
fn effective_columns(
    params: &[f64],
    fixed: &[bool],
    entities: &[SpaceEntity],
    touched: &[usize],
) -> Vec<Vec<(usize, f64)>> {
    let mut consumed = vec![false; touched.len()];
    let mut columns: Vec<Vec<(usize, f64)>> = Vec::new();
    for e in entities {
        if e.kind != SpaceEntityKind::Plane {
            continue;
        }
        let idx = [e.param(0), e.param(1), e.param(2), e.param(3)];
        let local: Vec<(usize, usize)> = touched
            .iter()
            .enumerate()
            .filter_map(|(lc, &p)| idx.iter().position(|&i| i == p).map(|off| (lc, off)))
            .collect();
        if local.is_empty() {
            continue;
        }
        if idx.iter().any(|&i| fixed[i]) {
            // Partially fixed plane: keep the raw free coordinates.
            continue;
        }
        let Some((t1, t2)) = plane_tangents(params[idx[0]], params[idx[1]], params[idx[2]])
        else {
            continue;
        };
        let offset = (0.0, 0.0, 0.0, 1.0);
        let t1 = (t1.0, t1.1, t1.2, 0.0);
        let t2 = (t2.0, t2.1, t2.2, 0.0);
        for basis in [t1, t2, offset] {
            let mut col: Vec<(usize, f64)> = Vec::new();
            for &(lc, off) in &local {
                let coeff = match off {
                    0 => basis.0,
                    1 => basis.1,
                    2 => basis.2,
                    _ => basis.3,
                };
                if coeff != 0.0 {
                    col.push((lc, coeff));
                }
            }
            columns.push(col);
        }
        for &(lc, _) in &local {
            consumed[lc] = true;
        }
    }
    for (lc, used) in consumed.iter().enumerate() {
        if !used {
            columns.push(vec![(lc, 1.0)]);
        }
    }
    columns
}

/// Two unit vectors spanning the tangent plane of `n = (nx, ny, nz)`.
/// Deterministic: the helper axis is the coordinate axis least aligned with
/// the normal, so the basis does not jitter between equivalent normals.
fn plane_tangents(
    nx: f64,
    ny: f64,
    nz: f64,
) -> Option<((f64, f64, f64), (f64, f64, f64))> {
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if !(len > SOLVER_NORM_TOL) {
        return None;
    }
    let u = (nx / len, ny / len, nz / len);
    let ax = u.0.abs();
    let ay = u.1.abs();
    let az = u.2.abs();
    let helper = if ax <= ay && ax <= az {
        (1.0, 0.0, 0.0)
    } else if ay <= az {
        (0.0, 1.0, 0.0)
    } else {
        (0.0, 0.0, 1.0)
    };
    let c1 = cross(u, helper);
    let l1 = (c1.0 * c1.0 + c1.1 * c1.1 + c1.2 * c1.2).sqrt();
    if !(l1 > 0.0) {
        return None;
    }
    let t1 = (c1.0 / l1, c1.1 / l1, c1.2 / l1);
    let t2 = cross(u, t1);
    Some((t1, t2))
}

#[inline]
fn cross(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (
        a.1 * b.2 - a.2 * b.1,
        a.2 * b.0 - a.0 * b.2,
        a.0 * b.1 - a.1 * b.0,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space3d::constraint::SpaceConstraint;
    use crate::space3d::entity::SpacePointRef;
    use crate::space3d::sketch::SpaceSketch;

    #[test]
    fn qr_reveals_rank_and_pivots() {
        let mut a = vec![
            vec![1.0, 1.0, 0.0, 0.0],
            vec![0.0, 1.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0],
        ];
        let (pivots, rank) = householder_qr_column_pivoting(&mut a, 1.0e-10);
        assert_eq!(rank, 2);
        assert_eq!(pivots[..2], [1, 0]);
    }

    #[test]
    fn fully_constrained_plane_has_zero_dof() {
        // PlaneNormal (2 rows) + PointOnPlane (1 row) determine the plane's
        // three physical DOF; the fourth storage coordinate is gauge.
        let mut sk = SpaceSketch::new();
        let plane = sk.add_plane(0.3, -0.4, 0.5, 1.2);
        let p = sk.add_point(0.0, 0.0, 0.0);
        sk.fix_entity(p);
        sk.add_constraint(SpaceConstraint::PlaneNormal {
            plane,
            nx: 0.0,
            ny: 0.0,
            nz: 1.0,
        });
        sk.add_constraint(SpaceConstraint::PointOnPlane {
            point: SpacePointRef::Point(p),
            plane,
        });
        let d = sk.diagnose();
        assert_eq!(d.dof, 0, "plane must be fully constrained: {d:?}");
        assert!(d.conflicting.is_empty(), "{d:?}");
        assert!(d.redundant.is_empty(), "gauge must not surface as a dependency: {d:?}");
    }

    #[test]
    fn duplicate_sphere_radius_is_redundant() {
        let mut sk = SpaceSketch::new();
        let s = sk.add_sphere(0.0, 0.0, 0.0, 2.0);
        sk.add_constraint(SpaceConstraint::SphereRadius { sphere: s, radius: 2.0 });
        sk.add_constraint(SpaceConstraint::SphereRadius { sphere: s, radius: 2.0 });
        let d = sk.diagnose();
        assert!(!d.redundant.is_empty(), "{d:?}");
        assert!(d.conflicting.is_empty(), "{d:?}");
        assert!(
            d.dependency_groups.iter().any(|g| g == &vec![0, 1]),
            "duplicate rows form one dependency group: {d:?}"
        );
    }

    #[test]
    fn contradiction_is_conflicting() {
        let mut sk = SpaceSketch::new();
        let p = sk.add_point(0.0, 0.0, 0.0);
        sk.add_constraint(SpaceConstraint::fix_point(p, 0.0, 0.0, 0.0));
        sk.add_constraint(SpaceConstraint::fix_point(p, 1.0, 1.0, 1.0));
        let d = sk.diagnose();
        assert!(!d.conflicting.is_empty(), "{d:?}");
    }

}
