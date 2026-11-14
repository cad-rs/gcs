//! Per-parameter mobility from the linearized constraint system.
//!
//! A parameter is **movable** when it has a non-negligible component in at
//! least one null-space vector of the Jacobian at the current configuration:
//! there is an infinitesimal motion satisfying every constraint that changes
//! it. Masked parameters (`fixed[i] == true`) and parameters pinned by their
//! own independent rows are not movable.
//!
//! The analysis runs per connected component of the constraint graph, so the
//! cost scales with component size. It is a commit/solve-time signal (the
//! sketcher caches it after a converged solve), not a per-frame full-system
//! QR.
//!
//! This is deliberately stronger than "the component has DOF": in
//! `x = 5, y + z = 7` the system has one DOF, but `x` is pinned — only `y`
//! and `z` are movable. The null-space support distinguishes them.
//!
//! Performance: the Jacobian is sparse (each constraint touches only a few
//! parameters), so the RREF uses a sparse row representation that skips zero
//! entries. This reduces the per-component cost from O(n^3) to O(n * nnz)
//! where nnz is the number of non-zero entries.

use crate::constraint::Constraint;
use crate::diagnosis::RANK_REL_TOL;
use crate::entity::Entity;
use crate::solver::SketchEquations;
use crate::system::{self, Equations};
use std::collections::HashMap;

/// Marks every parameter that can move at the current configuration.
///
/// Returns one flag per entry of `params`. Masked parameters are never
/// movable; free parameters that no constraint touches are movable (nothing
/// holds them).
///
/// Performance: uses a sparse RREF that skips zero entries in the Jacobian.
/// The sparse representation stores each row as a Vec of (column, value)
/// pairs, so the elimination cost is proportional to the number of non-zero
/// entries rather than the full matrix size.
pub fn mobility(
    params: &[f64],
    fixed: &[bool],
    entities: &[Entity],
    constraints: &[Constraint],
) -> Vec<bool> {
    mobility_of(&SketchEquations::new(entities, constraints), params, fixed)
}

/// [`mobility`] of any set of equations over plain (additively stepped)
/// parameters.
pub(crate) fn mobility_of<E: Equations + ?Sized>(
    eqs: &E,
    params: &[f64],
    fixed: &[bool],
) -> Vec<bool> {
    let mut movable: Vec<bool> = fixed.iter().map(|fixed| !fixed).collect();
    if params.is_empty() {
        return movable;
    }

    for group in system::components(eqs, params) {
        let n_eq = system::equation_total(eqs, &group);
        if n_eq == 0 {
            continue;
        }

        // Touched free parameters of this component, in ascending index
        // order (the Jacobian touch set includes structural zeros).
        let mut touched = system::touched_params(eqs, &group, params);
        touched.retain(|&param| !fixed[param]);
        if touched.is_empty() {
            continue;
        }

        let col_of: HashMap<usize, usize> = touched
            .iter()
            .enumerate()
            .map(|(col, &param)| (param, col))
            .collect();
        let mut jacobian = vec![vec![0.0_f64; touched.len()]; n_eq];
        system::analytic_jacobian(eqs, &group, params, &col_of, &mut jacobian);

        // Convert to sparse representation: each row stores only non-zero
        // entries as (column, value) pairs.
        let sparse: Vec<Vec<(usize, f64)>> = jacobian
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(_, &v)| v != 0.0)
                    .map(|(col, &v)| (col, v))
                    .collect()
            })
            .collect();

        let support = null_space_support_sparse(&sparse, touched.len(), RANK_REL_TOL);
        for (col, &param) in touched.iter().enumerate() {
            movable[param] = support[col];
        }
    }
    movable
}

/// Sparse null-space support computation.
///
/// Takes a sparse matrix (each row is a Vec of (column, value) pairs),
/// the number of columns, and returns a Vec<bool> indicating which columns
/// have a non-negligible component in at least one null-space vector.
///
/// The algorithm mirrors `numerics::null_space_support` but operates on the
/// sparse representation, skipping zero entries during elimination. This
/// reduces the cost from O(n^3) to O(n * nnz) where nnz is the number of
/// non-zero entries.
fn null_space_support_sparse(j: &[Vec<(usize, f64)>], n_cols: usize, rel_tol: f64) -> Vec<bool> {
    let mut movable = vec![false; n_cols];
    if n_cols == 0 || j.is_empty() {
        return movable;
    }

    // Row equilibration: scale each row so its largest entry is 1.0.
    let mut r: Vec<Vec<(usize, f64)>> = j.to_vec();
    for row in &mut r {
        let scale = row.iter().fold(0.0_f64, |m, (_, v)| m.max(v.abs()));
        if scale > 0.0 {
            for (_, v) in row.iter_mut() {
                *v /= scale;
            }
        }
    }

    let n_rows = r.len();
    let mut remaining: Vec<usize> = (0..n_cols).collect();
    let mut pivot_cols: Vec<usize> = Vec::new();
    let mut pivot_row = 0usize;

    while pivot_row < n_rows && !remaining.is_empty() {
        // Pivot on the remaining column with the largest residual norm.
        let mut best_pos = 0usize;
        let mut best = -1.0_f64;
        for (pos, &col) in remaining.iter().enumerate() {
            let norm = r[pivot_row..]
                .iter()
                .map(|row| {
                    row.iter()
                        .filter(|(c, _)| *c == col)
                        .map(|(_, v)| v * v)
                        .sum::<f64>()
                })
                .sum::<f64>()
                .sqrt();
            if norm > best {
                best = norm;
                best_pos = pos;
            }
        }
        if best <= rel_tol {
            break;
        }
        let pcol = remaining.remove(best_pos);

        // Bring the row with the largest entry in the pivot column up to
        // the pivot position: the column was chosen over all the rows left,
        // and the row that happens to stand there may not contain it at all
        // (eliminating with such a row would subtract rows that have no
        // pivot to cancel).
        let entry = |row: &[(usize, f64)]| {
            row.iter()
                .find(|(c, _)| *c == pcol)
                .map_or(0.0, |(_, v)| v.abs())
        };
        let best_row = (pivot_row..n_rows)
            .max_by(|&a, &b| entry(&r[a]).total_cmp(&entry(&r[b])))
            .expect("rows remain");
        r.swap(pivot_row, best_row);

        // Normalize the pivot row.
        let pivot = r[pivot_row]
            .iter()
            .find(|(c, _)| *c == pcol)
            .map(|(_, v)| *v)
            .unwrap_or(1.0);
        if pivot != 0.0 {
            for (_, v) in r[pivot_row].iter_mut() {
                *v /= pivot;
            }
        }

        // Eliminate the pivot column from every other row (RREF).
        for i in 0..n_rows {
            if i == pivot_row {
                continue;
            }
            let factor = r[i]
                .iter()
                .find(|(c, _)| *c == pcol)
                .map(|(_, v)| *v)
                .unwrap_or(0.0);
            if factor != 0.0 {
                // Sparse elimination: only touch columns that appear in
                // either row.
                let mut new_row: Vec<(usize, f64)> = Vec::new();
                let mut j_idx = 0;
                let mut p_idx = 0;
                while j_idx < r[i].len() || p_idx < r[pivot_row].len() {
                    let j_col = r[i].get(j_idx).map(|(c, _)| *c);
                    let p_col = r[pivot_row].get(p_idx).map(|(c, _)| *c);
                    match (j_col, p_col) {
                        (Some(jc), Some(pc)) if jc == pc => {
                            let v = r[i][j_idx].1 - factor * r[pivot_row][p_idx].1;
                            if v != 0.0 {
                                new_row.push((jc, v));
                            }
                            j_idx += 1;
                            p_idx += 1;
                        }
                        (Some(jc), Some(pc)) if jc < pc => {
                            new_row.push(r[i][j_idx]);
                            j_idx += 1;
                        }
                        (Some(_), Some(pc)) => {
                            let v = -factor * r[pivot_row][p_idx].1;
                            if v != 0.0 {
                                new_row.push((pc, v));
                            }
                            p_idx += 1;
                        }
                        (Some(_), None) => {
                            new_row.push(r[i][j_idx]);
                            j_idx += 1;
                        }
                        (None, Some(_)) => {
                            let v = -factor * r[pivot_row][p_idx].1;
                            if v != 0.0 {
                                new_row.push((p_col.unwrap(), v));
                            }
                            p_idx += 1;
                        }
                        (None, None) => break,
                    }
                }
                r[i] = new_row;
            }
        }
        pivot_cols.push(pcol);
        pivot_row += 1;
    }

    // A column is movable if it is free (not a pivot column) or if it has a
    // non-negligible component in at least one null-space basis vector.
    for &free in &remaining {
        movable[free] = true;
        for (row, &pivot) in pivot_cols.iter().enumerate() {
            let val = r[row]
                .iter()
                .find(|(c, _)| *c == free)
                .map(|(_, v)| v.abs())
                .unwrap_or(0.0);
            if val > rel_tol {
                movable[pivot] = true;
            }
        }
    }
    movable
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::PointRef;
    use crate::sketch::Sketch;

    #[test]
    fn fully_constrained_point_is_locked_and_free_point_is_movable() {
        let mut sk = Sketch::new();
        let anchor = sk.add_point(0.0, 0.0);
        let _free = sk.add_point(5.0, 7.0);
        sk.add_constraint(Constraint::fix_point(PointRef::Point(anchor), 0.0, 0.0));

        let m = mobility(&sk.params, &sk.fixed, &sk.entities, &sk.constraints);
        assert_eq!(&m[0..2], &[false, false], "the anchored point is locked");
        assert_eq!(
            &m[2..4],
            &[true, true],
            "the untouched point is fully movable"
        );
    }

    #[test]
    fn an_x_driven_point_is_locked_only_along_x() {
        let mut sk = Sketch::new();
        let anchor = sk.add_point(0.0, 0.0);
        let driven = sk.add_point(5.0, 7.0);
        sk.add_constraint(Constraint::fix_point(PointRef::Point(anchor), 0.0, 0.0));
        sk.add_constraint(Constraint::distance_x(
            PointRef::Point(driven),
            PointRef::Point(anchor),
            5.0,
        ));

        let m = mobility(&sk.params, &sk.fixed, &sk.entities, &sk.constraints);
        assert_eq!(
            &m[0..2],
            &[false, false],
            "the anchor is pinned in both axes"
        );
        assert!(!m[2], "x is fixed by the DistanceX row");
        assert!(m[3], "y remains free");
    }

    #[test]
    fn an_anchor_stays_locked_under_an_open_chain() {
        // The anchor's rows hold one column each; the pivot column chosen
        // first is a chain point's, which the anchor's row does not
        // contain. Eliminating with it used to mark the anchor's y movable.
        let mut sk = Sketch::new();
        let anchor = sk.add_point(0.0, 0.0);
        let p1 = sk.add_point(0.96, 0.29);
        let p2 = sk.add_point(1.92, 0.57);
        sk.add_constraint(Constraint::fix_point(PointRef::Point(anchor), 0.0, 0.0));
        sk.add_constraint(Constraint::point_distance(anchor, p1, 1.0));
        sk.add_constraint(Constraint::point_distance(p1, p2, 1.0));

        let m = mobility(&sk.params, &sk.fixed, &sk.entities, &sk.constraints);
        assert_eq!(&m[0..2], &[false, false], "the anchor is pinned: {m:?}");
        assert!(
            m[2..6].iter().all(|&movable| movable),
            "the chain swings: {m:?}"
        );
    }

    #[test]
    fn masked_parameters_are_never_movable() {
        let mut sk = Sketch::new();
        let line = sk.add_line(0.0, 0.0, 10.0, 0.0);
        let other = sk.add_line(0.0, 5.0, 10.0, 5.0);
        sk.fix_entity(line);

        let m = mobility(&sk.params, &sk.fixed, &sk.entities, &sk.constraints);
        assert_eq!(&m[0..4], &[false, false, false, false], "masked line");
        assert_eq!(&m[4..8], &[true, true, true, true], "free line");
        let _ = other;
    }

    #[test]
    fn a_free_horizontal_line_moves_along_its_free_directions() {
        let mut sk = Sketch::new();
        let line = sk.add_line(0.0, 0.0, 10.0, 0.0);
        sk.add_constraint(Constraint::Horizontal(line));

        let m = mobility(&sk.params, &sk.fixed, &sk.entities, &sk.constraints);
        assert_eq!(m.len(), 4);
        assert!(m.iter().all(|&movable| movable), "H leaves 3 DOF: {m:?}");
    }

    #[test]
    fn no_constraints_means_everything_is_movable() {
        let mut sk = Sketch::new();
        sk.add_point(1.0, 2.0);
        sk.add_line(0.0, 0.0, 1.0, 1.0);
        let m = mobility(&sk.params, &sk.fixed, &sk.entities, &sk.constraints);
        assert!(m.iter().all(|&movable| movable));
    }
}
