//! How a driven (reference) dimension follows the driving ones.
//!
//! A driven dimension is a dimensional constraint that is measured, not
//! enforced: `m(x)`, the value its residual would need. Where the system
//! stands (solved), the driving dimensions `v` hold `r(x, v) = 0`; moving
//! one of them moves the solution by `dx` with `J dx = −∂r/∂v · dv`, and
//! the driven one by `∇m · dx`. When `∇m` lies in the row space of `J`
//! (`∇m = Jᵀ y`) that change is the same whichever solution the motion
//! picks: `dm = −yᵀ (∂r/∂v) dv`, so one least-squares solve gives the
//! derivative with respect to every driving dimension at once. When it
//! does not, the driven dimension can change with no driving one changing
//! — it is not determined — and there is no derivative to give. Over the
//! effective columns of the gauge, like the diagnosis. Where driving
//! dimensions depend on each other (redundant), they cannot change alone;
//! the derivative is split between them by the least-squares solution of
//! smallest norm.

use std::collections::HashMap;

use crate::system::{self, Equations, Gauge, SparseRow};

/// `∇m` must be reproduced by the rows to this, relative to its size, to
/// count as determined.
const DETERMINED_TOL: f64 = 1.0e-6;

/// The derivative of a driven dimension with respect to every active
/// driving dimension of `eqs` at `params`: `driven_row` is the gradient of
/// the driven dimension's residual over the stored parameters,
/// `driven_slope` its derivative with respect to the driven value, and
/// `slope_of(j)` the derivative of constraint `j`'s residual with respect
/// to its value (`None` for a constraint without one). `None` when the
/// driven dimension is not determined; otherwise `(constraint, ∂m/∂v)`
/// for every driving dimension it depends on, ascending.
pub(crate) fn sensitivities<E: Equations + ?Sized, G: Gauge>(
    eqs: &E,
    params: &[f64],
    fixed: &[bool],
    gauge: &G,
    driven_row: &SparseRow,
    driven_slope: f64,
    slope_of: &dyn Fn(usize) -> Option<f64>,
) -> Option<Vec<(usize, f64)>> {
    if driven_slope == 0.0 || !driven_slope.is_finite() {
        return None;
    }
    // ∇m over the free parameters.
    let mut gradient: HashMap<usize, f64> = HashMap::new();
    for &(p, v) in driven_row {
        if !fixed[p] {
            *gradient.entry(p).or_insert(0.0) -= v / driven_slope;
        }
    }
    gradient.retain(|_, v| *v != 0.0);
    if gradient.is_empty() {
        return Some(Vec::new());
    }
    // The components the driven dimension reaches.
    let components = system::components(eqs, params);
    let mut component_of = vec![usize::MAX; eqs.column_count(params)];
    for (k, group) in components.iter().enumerate() {
        for p in system::touched_params(eqs, group, params) {
            component_of[p] = k;
        }
    }
    let mut merged: Vec<usize> = Vec::new();
    for &p in gradient.keys() {
        let c = component_of[p];
        if c == usize::MAX {
            // A free parameter no constraint touches: it moves freely.
            return None;
        }
        merged.push(c);
    }
    merged.sort_unstable();
    merged.dedup();
    let group: Vec<usize> = merged
        .iter()
        .flat_map(|&c| components[c].iter().copied())
        .collect();
    let mut touched = system::touched_params(eqs, &group, params);
    touched.retain(|&p| !fixed[p]);
    let col_of: HashMap<usize, usize> = touched.iter().enumerate().map(|(c, &p)| (p, c)).collect();
    let n_eq = system::equation_total(eqs, &group);
    let mut j_raw = vec![vec![0.0_f64; touched.len()]; n_eq];
    system::analytic_jacobian(eqs, &group, params, &col_of, &mut j_raw);
    let columns: Vec<SparseRow> = gauge
        .effective_columns(params, fixed, &touched)
        .unwrap_or_else(|| (0..touched.len()).map(|c| vec![(c, 1.0)]).collect());
    let n_cols = columns.len();
    // Rows and ∇m over the effective columns; rows equilibrated.
    let mut scale = vec![0.0_f64; n_eq];
    let rows: Vec<Vec<f64>> = j_raw
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r: Vec<f64> = columns
                .iter()
                .map(|col| col.iter().map(|&(lc, w)| row[lc] * w).sum())
                .collect();
            let s = r.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            if s > 0.0 {
                for v in &mut r {
                    *v /= s;
                }
            }
            scale[i] = s;
            r
        })
        .collect();
    let g: Vec<f64> = columns
        .iter()
        .map(|col| {
            col.iter()
                .map(|&(lc, w)| gradient.get(&touched[lc]).copied().unwrap_or(0.0) * w)
                .sum()
        })
        .collect();
    let g_norm = g.iter().map(|v| v * v).sum::<f64>().sqrt();
    if g_norm == 0.0 {
        return Some(Vec::new());
    }
    // Jᵀ y = ∇m without forming normal equations (a near-singular row
    // squares badly): an orthonormal basis Q of the rows (pivoted
    // Gram–Schmidt, as the diagnosis's rank) says whether ∇m is in their
    // span; then with L = J Qᵀ the smallest y with Lᵀ y = Q ∇m comes from
    // a thin QR of L.
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let max0 = rows
        .iter()
        .map(|r| dot(r, r).sqrt())
        .fold(0.0_f64, f64::max);
    let tol = crate::diagnosis::RANK_REL_TOL * max0;
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut rest: Vec<Vec<f64>> = rows.clone();
    loop {
        let (best, norm) = rest
            .iter()
            .enumerate()
            .map(|(i, r)| (i, dot(r, r).sqrt()))
            .fold((0, -1.0), |acc, x| if x.1 > acc.1 { x } else { acc });
        if norm <= tol {
            break;
        }
        let u: Vec<f64> = rest[best].iter().map(|v| v / norm).collect();
        for r in &mut rest {
            for _ in 0..2 {
                let d = dot(r, &u);
                for (a, b) in r.iter_mut().zip(&u) {
                    *a -= d * b;
                }
            }
        }
        q.push(u);
    }
    let c: Vec<f64> = q.iter().map(|u| dot(u, &g)).collect();
    let miss = {
        let mut left = g.clone();
        for (u, &ci) in q.iter().zip(&c) {
            for (a, b) in left.iter_mut().zip(u) {
                *a -= ci * b;
            }
        }
        dot(&left, &left).sqrt()
    };
    if miss > DETERMINED_TOL * g_norm {
        return None;
    }
    // L = J Qᵀ (rows by basis vectors), its columns orthonormalized:
    // L = Q_L R_L; then y = Q_L R_L⁻ᵀ c.
    let r = q.len();
    let mut l_cols: Vec<Vec<f64>> = (0..r)
        .map(|k| rows.iter().map(|row| dot(row, &q[k])).collect())
        .collect();
    let mut r_l = vec![vec![0.0_f64; r]; r];
    for k in 0..r {
        for i in 0..k {
            let d = dot(&l_cols[k], &l_cols[i]);
            r_l[i][k] += d;
            let (head, tail) = l_cols.split_at_mut(k);
            for (a, b) in tail[0].iter_mut().zip(&head[i]) {
                *a -= d * b;
            }
        }
        for i in 0..k {
            let d = dot(&l_cols[k], &l_cols[i]);
            r_l[i][k] += d;
            let (head, tail) = l_cols.split_at_mut(k);
            for (a, b) in tail[0].iter_mut().zip(&head[i]) {
                *a -= d * b;
            }
        }
        let n = dot(&l_cols[k], &l_cols[k]).sqrt();
        r_l[k][k] = n;
        if n > 0.0 {
            for v in &mut l_cols[k] {
                *v /= n;
            }
        }
    }
    // R_Lᵀ t = c, forward.
    let mut t = vec![0.0_f64; r];
    for k in 0..r {
        let s: f64 = (0..k).map(|i| r_l[i][k] * t[i]).sum();
        t[k] = if r_l[k][k] != 0.0 {
            (c[k] - s) / r_l[k][k]
        } else {
            0.0
        };
    }
    let mut y = vec![0.0_f64; n_eq];
    for (k, &tk) in t.iter().enumerate() {
        for (yi, li) in y.iter_mut().zip(&l_cols[k]) {
            *yi += tk * li;
        }
    }
    let _ = n_cols;

    // dm/dv_j = −y_j · ∂r_j/∂v_j for every single-equation dimension.
    let mut out = Vec::new();
    let mut row = 0;
    for &c in &group {
        let count = eqs.equation_count(c);
        if count == 1 {
            if let Some(slope) = slope_of(c) {
                if scale[row] > 0.0 {
                    let d = -(y[row] / scale[row]) * slope;
                    if d.abs() > 1.0e-12 {
                        out.push((c, d));
                    }
                }
            }
        }
        row += count;
    }
    out.sort_unstable_by_key(|&(c, _)| c);
    Some(out)
}

/// The derivative of a residual with respect to the value it is driven
/// to, by central differences of `residual(value)`.
pub(crate) fn value_slope(value: f64, residual: &dyn Fn(f64) -> f64) -> f64 {
    let h = 1.0e-6 * value.abs().max(1.0);
    (residual(value + h) - residual(value - h)) / (2.0 * h)
}
