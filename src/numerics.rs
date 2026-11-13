//! The shared Newton-Raphson least-squares core.
//!
//! Every solver family (2D sketch, 3D geometric, 3D rigid pose) reduces to
//! the same loop over a flat parameter vector: evaluate residuals, build the
//! Jacobian by central finite differences, solve the regularised normal
//! equations `(JᵀJ + λI) Δx = Jᵀ(−F)`, and take a damped step. Only the
//! residual evaluation differs per family, so it is injected as a closure.
//!
//! The Moore-Penrose pseudo-inverse J⁺ is realised through the normal
//! equations with Tikhonov regularisation λ: under-constrained systems get
//! the minimum-norm step instead of diverging.

/// Result returned by every [`solve`][`crate::solver::solve`] entry point.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub struct SolveResult {
    /// `true` if all constraint residuals are below [`RESIDUAL_TOL`].
    pub converged: bool,
    /// RMS residual at termination.
    pub residual: f64,
    /// Number of Newton iterations performed.
    pub iterations: usize,
}

/// Convergence tolerance on the RMS constraint residual.
pub const RESIDUAL_TOL: f64 = 1e-10;
/// Maximum Newton iterations.
pub const MAX_ITER: usize = 100;
/// Step size for numerical Jacobian (central differences).
pub const FD_H: f64 = 1e-7;
/// Tikhonov regularisation coefficient.
pub const LAMBDA: f64 = 1e-10;
/// Pivot threshold for Gaussian elimination.
pub const PIVOT_TOL: f64 = 1e-14;
/// Normalization guard threshold — prevents division by near-zero vectors.
pub const SOLVER_NORM_TOL: f64 = 1e-10;
/// Relative non-worsening tolerance for the damped line search (Armijo-type:
/// accept a step whose residual is not larger than the old one by this factor).
pub const LINE_SEARCH_RELAX: f64 = 1e-4;
/// Halving attempts of the damped line search.
pub const LINE_SEARCH_STEPS: usize = 8;
/// Bounded strictly-improving convergence-polish steps after convergence
/// (unregularised, accepted only on strict residual improvement). Each
/// unregularised refinement lands within linear-solve rounding (~1e-16
/// relative) of the fixed point, collapsing the RESIDUAL_TOL-scale
/// terminal path-dependence by ~64 orders. 0.0-targets have a zero-width
/// absorption window and keep descending a sub-noise power-of-two ladder
/// (measured: ~1e-6.4 per rung) — the representation last mile for those
/// belongs to the persistence zero-band (rcad-sketcher write-back), not
/// to more polish rungs. The strict-improvement guard terminates early
/// on stalls (finite-difference Jacobian families stall near their
/// ~1e-8 noise floor after a couple of rungs).
pub const POLISH_STEPS: usize = 4;

/// Run the damped Newton-Raphson loop over a flat parameter vector.
///
/// - `params`: full parameter vector (free entries updated in place).
/// - `fixed`: mask — `fixed[i] == true` means `params[i]` is held constant.
/// - `n_eq`: total number of constraint equations (= residual components).
/// - `eval_residuals`: writes the `n_eq` residual components for a given
///   parameter vector into the output buffer (called with the buffer sized
///   `n_eq`).
pub fn newton_solve<F>(
    params: &mut Vec<f64>,
    fixed: &[bool],
    n_eq: usize,
    mut eval_residuals: F,
) -> SolveResult
where
    F: FnMut(&[f64], &mut [f64]),
{
    newton_solve_with_jacobian(
        params,
        fixed,
        n_eq,
        &mut eval_residuals,
        Jacobian::FiniteDifference,
    )
}

/// Where the Newton loop sources its Jacobian each iteration.
pub(crate) enum Jacobian<'a> {
    /// Central finite differences over the residual closure (the fallback
    /// for solver families without analytic derivatives yet).
    FiniteDifference,
    /// Analytic assembly (solver roadmap M2): the closure writes the
    /// m × n_free dense matrix, pre-zeroed by the loop, at the given
    /// parameters.
    Analytic(&'a mut dyn FnMut(&[f64], &mut [Vec<f64>])),
}

/// [`newton_solve`] with an explicit Jacobian source.
pub fn newton_solve_with_jacobian<F>(
    params: &mut Vec<f64>,
    fixed: &[bool],
    n_eq: usize,
    mut eval_residuals: F,
    jacobian: Jacobian<'_>,
) -> SolveResult
where
    F: FnMut(&[f64], &mut [f64]),
{
    let free_params: Vec<usize> = (0..params.len()).filter(|&i| !fixed[i]).collect();
    let n_free = free_params.len();

    if n_eq == 0 || n_free == 0 {
        return SolveResult {
            converged: true,
            residual: 0.0,
            iterations: 0,
        };
    }

    let mut residual;
    let mut iters = 0;
    // Hoisted once; re-filled every iteration (FD overwrites fully, the
    // analytic builder receives it zeroed).
    let mut j = vec![vec![0.0_f64; n_free]; n_eq];
    // The Jacobian source is matched once (it is consumed by value); the
    // per-iteration branch checks the option.
    let mut analytic: Option<&mut dyn FnMut(&[f64], &mut [Vec<f64>])> = match jacobian {
        Jacobian::FiniteDifference => None,
        Jacobian::Analytic(build) => Some(build),
    };

    for _ in 0..MAX_ITER {
        let f = eval_vec(&mut eval_residuals, params, n_eq);
        residual = rms(&f);
        if residual < RESIDUAL_TOL {
            break;
        }

        // Jacobian J (n_eq × n_free).
        if let Some(build) = analytic.as_mut() {
            for row in j.iter_mut() {
                row.fill(0.0);
            }
            build(params, &mut j);
        } else {
            numerical_jacobian_into(params, &free_params, n_eq, &mut eval_residuals, &mut j);
        }

        // Solve (J^T J + λI) Δx = J^T (−F)  →  Δx
        let delta = solve_normal_equations(&j, &f, n_free);

        // Damped update: halve the step while the residual grows.
        let mut alpha = 1.0_f64;
        let old_params: Vec<f64> = free_params.iter().map(|&gi| params[gi]).collect();
        for _ in 0..LINE_SEARCH_STEPS {
            for (k, &gi) in free_params.iter().enumerate() {
                params[gi] = old_params[k] + alpha * delta[k];
            }
            let f_new = eval_vec(&mut eval_residuals, params, n_eq);
            if rms(&f_new) < residual * (1.0 + LINE_SEARCH_RELAX) {
                break;
            }
            alpha *= 0.5;
        }

        iters += 1;
    }

    let mut f_final = eval_vec(&mut eval_residuals, params, n_eq);
    residual = rms(&f_final);

    // Convergence polish (solver roadmap follow-up): the loop stops at
    // the first iterate under RESIDUAL_TOL, which is path-dependent at
    // that scale. Accepted-if-not-worse steps on the UNREGULARISED normal
    // equations nail nonzero exactly-representable fixed points to
    // machine zero (half-ulp absorption), collapsing the terminal
    // path-dependence that used to change serialized floats. The main
    // loop's Tikhonov λ would shorten every step by ~λ/σ², and even
    // unregularised single steps carry linear-solve rounding (~1e-16
    // relative), which strands 0.0-targets — zero-width absorption
    // window — on a sub-noise power-of-two ladder (measured: 2^-110)
    // instead of exact zeros; the representation last mile for those is
    // the persistence zero band in the rcad-sketcher write-back. Hence
    // the bounded strictly-improving run. Singular systems
    // (rank-deficient or under-constrained components) yield a zero step
    // and the polish declines to act.
    let mut polish_steps = 0usize;
    while residual > 0.0 && residual < RESIDUAL_TOL && polish_steps < POLISH_STEPS {
        if let Some(build) = analytic.as_mut() {
            for row in j.iter_mut() {
                row.fill(0.0);
            }
            build(params, &mut j);
        } else {
            numerical_jacobian_into(params, &free_params, n_eq, &mut eval_residuals, &mut j);
        }
        // Unregularised step: see the block comment above.
        let delta = solve_normal_equations_with_lambda(&j, &f_final, n_free, 0.0);
        let old_params: Vec<f64> = free_params.iter().map(|&gi| params[gi]).collect();
        for (k, &gi) in free_params.iter().enumerate() {
            params[gi] = old_params[k] + delta[k];
        }
        let f_polish = eval_vec(&mut eval_residuals, params, n_eq);
        let residual_polish = rms(&f_polish);
        if residual_polish < residual {
            residual = residual_polish;
            f_final = f_polish;
            polish_steps += 1;
        } else {
            // Worsened (or NaN, or no strict gain) — reject the polish
            // step, restore, and stop.
            for (k, &gi) in free_params.iter().enumerate() {
                params[gi] = old_params[k];
            }
            break;
        }
    }

    SolveResult {
        converged: residual < RESIDUAL_TOL,
        residual,
        iterations: iters,
    }
}

fn eval_vec<F>(eval: &mut F, params: &[f64], n_eq: usize) -> Vec<f64>
where
    F: FnMut(&[f64], &mut [f64]),
{
    let mut f = vec![0.0_f64; n_eq];
    eval(params, &mut f);
    f
}

/// Jacobian `j[row][col] = ∂F_row / ∂x_col` by central finite differences,
/// one parameter column at a time, written into the caller's buffer.
pub(crate) fn numerical_jacobian_into<F>(
    params: &[f64],
    free_params: &[usize],
    n_eq: usize,
    eval_residuals: &mut F,
    out: &mut [Vec<f64>],
) where
    F: FnMut(&[f64], &mut [f64]),
{
    let mut p = params.to_vec();

    for (col, &gi) in free_params.iter().enumerate() {
        let orig = p[gi];

        p[gi] = orig + FD_H;
        let f_plus = eval_vec(eval_residuals, &p, n_eq);

        p[gi] = orig - FD_H;
        let f_minus = eval_vec(eval_residuals, &p, n_eq);

        p[gi] = orig;

        for (row, o) in out.iter_mut().enumerate() {
            o[col] = (f_plus[row] - f_minus[row]) / (2.0 * FD_H);
        }
    }
}

/// Householder QR with column pivoting (Golub & Van Loan, Alg. 5.4.1),
/// rank-revealing: the matrix `a` (n_rows × n_cols) is reduced in place,
/// and columns are permuted so that the leading `rank` columns are the
/// linearly independent ones (each step pivots the column with the largest
/// remaining norm). Returns `(pivots, rank)` where `pivots[j]` is the
/// original column index now at position j.
///
/// `rel_tol` is relative to the largest initial column norm: a column whose
/// remaining norm falls at or below `rel_tol × max_initial` ends the
/// factorization. With finite-difference Jacobians the threshold must sit
/// above the differentiation noise floor (~1e-8 relative) — see
/// `diagnosis::RANK_REL_TOL`.
pub(crate) fn householder_qr_column_pivoting(
    a: &mut [Vec<f64>],
    rel_tol: f64,
) -> (Vec<usize>, usize) {
    let n_rows = a.len();
    let n_cols = if n_rows == 0 { 0 } else { a[0].len() };
    let mut pivots: Vec<usize> = (0..n_cols).collect();

    let col_norm_from = |a: &[Vec<f64>], j: usize, from: usize| -> f64 {
        a[from..]
            .iter()
            .map(|row| row[j] * row[j])
            .sum::<f64>()
            .sqrt()
    };

    let max0 = (0..n_cols)
        .map(|j| col_norm_from(a, j, 0))
        .fold(0.0_f64, f64::max);
    if !(max0 > 0.0) {
        return (pivots, 0);
    }

    let mut rank = 0usize;
    for k in 0..n_rows.min(n_cols) {
        // Pivot: the column with the largest remaining norm.
        let mut p = k;
        let mut best = col_norm_from(a, k, k);
        for j in (k + 1)..n_cols {
            let nj = col_norm_from(a, j, k);
            if nj > best {
                best = nj;
                p = j;
            }
        }
        if best <= rel_tol * max0 {
            break;
        }

        a[k..].iter_mut().for_each(|row| row.swap(k, p));
        pivots.swap(k, p);

        // Householder reflector H = I − 2 v vᵀ / vᵀv over rows k.., column k.
        let sigma: f64 = a[k..].iter().map(|row| row[k] * row[k]).sum::<f64>().sqrt();
        let alpha = if a[k][k] > 0.0 { -sigma } else { sigma };
        let mut v: Vec<f64> = a[k..].iter().map(|row| row[k]).collect();
        v[0] -= alpha;
        let v_norm2: f64 = v.iter().map(|x| x * x).sum();
        if v_norm2 > 0.0 {
            for j in (k + 1)..n_cols {
                let dot: f64 = a[k..]
                    .iter()
                    .zip(v.iter())
                    .map(|(row, vi)| row[j] * vi)
                    .sum();
                let beta = 2.0 * dot / v_norm2;
                for (row, vi) in a[k..].iter_mut().zip(v.iter()) {
                    row[j] -= beta * vi;
                }
            }
        }
        a[k][k] = alpha;
        for row in a[(k + 1)..].iter_mut() {
            row[k] = 0.0;
        }
        rank += 1;
    }
    (pivots, rank)
}

/// Per-column mobility from the null space of `j` (rows × cols).
///
/// The Jacobian is row-equilibrated, then reduced to RREF with column
/// pivoting: the pivot columns span the row space and every remaining column
/// is free. Each free column yields one null-space basis vector (1 at the
/// free column, the negated RREF entry at every pivot column), so a parameter
/// is movable iff it has a non-negligible component in at least one basis
/// vector.
///
/// The distinction matters: for `x = 5, y + z = 7` the columns of `y` and `z`
/// are coupled by one free direction, while `x` is pinned by its own pivot
/// column even though the system still has one DOF. `rel_tol` is relative to
/// the equilibrated RREF entries (exact dependencies land at roundoff).
pub(crate) fn null_space_support(j: &[Vec<f64>], rel_tol: f64) -> Vec<bool> {
    let n_cols = j.first().map_or(0, |row| row.len());
    let mut movable = vec![false; n_cols];
    if n_cols == 0 || j.is_empty() {
        return movable;
    }
    // Row equilibration: a distance row in mm and an angle row in radians
    // then contribute comparable pivot norms.
    let mut r = j.to_vec();
    for row in &mut r {
        let scale = row.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        if scale > 0.0 {
            for v in row.iter_mut() {
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
                .map(|row| row[col] * row[col])
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
        let pivot = r[pivot_row][pcol];
        for v in r[pivot_row].iter_mut() {
            *v /= pivot;
        }
        // Eliminate the pivot column from every other row (RREF).
        for i in 0..n_rows {
            if i == pivot_row {
                continue;
            }
            let factor = r[i][pcol];
            if factor != 0.0 {
                for col in 0..n_cols {
                    r[i][col] -= factor * r[pivot_row][col];
                }
            }
        }
        pivot_cols.push(pcol);
        pivot_row += 1;
    }
    for &free in &remaining {
        movable[free] = true;
        for (row, &pivot) in pivot_cols.iter().enumerate() {
            if r[row][free].abs() > rel_tol {
                movable[pivot] = true;
            }
        }
    }
    movable
}

#[cfg(test)]
mod null_space_tests {
    use super::null_space_support;

    #[test]
    fn pinned_column_is_not_movable() {
        let j = vec![vec![1.0, 0.0, 0.0]];
        assert_eq!(null_space_support(&j, 1e-9), vec![false, true, true]);
    }

    #[test]
    fn coupled_columns_are_both_movable() {
        let j = vec![vec![1.0, 1.0]];
        assert_eq!(null_space_support(&j, 1e-9), vec![true, true]);
    }

    #[test]
    fn pinned_pair_leaves_only_the_untouched_column() {
        let j = vec![vec![1.0, 1.0, 0.0], vec![1.0, -1.0, 0.0]];
        assert_eq!(null_space_support(&j, 1e-9), vec![false, false, true]);
    }

    #[test]
    fn zero_and_empty_rows_mark_every_column_movable() {
        let j = vec![vec![0.0, 0.0]];
        assert_eq!(null_space_support(&j, 1e-9), vec![true, true]);
        let empty: Vec<Vec<f64>> = Vec::new();
        assert!(null_space_support(&empty, 1e-9).is_empty());
    }
}

/// Regularised normal equations `(JᵀJ + λI) Δx = Jᵀ(−F)`; a singular system
/// (Gaussian elimination pivot below [`PIVOT_TOL`]) yields a zero step.
pub(crate) fn solve_normal_equations(j: &[Vec<f64>], f: &[f64], n: usize) -> Vec<f64> {
    solve_normal_equations_with_lambda(j, f, n, LAMBDA)
}

/// [`solve_normal_equations`] with an explicit regularisation strength;
/// `lambda == 0.0` solves the unregularised normal equations (the
/// convergence-polish step: any positive λ scales the step short of the
/// fixed point by ~λ/σ², so the terminal rounding never lands exactly on
/// it). A singular system yields a zero step.
pub(crate) fn solve_normal_equations_with_lambda(
    j: &[Vec<f64>],
    f: &[f64],
    n: usize,
    lambda: f64,
) -> Vec<f64> {
    let m = j.len();

    // A = J^T J + λI  (n×n)
    let mut a = vec![vec![0.0_f64; n]; n];
    // b = J^T (−F)    (n-vector)
    let mut b = vec![0.0_f64; n];

    for i in 0..m {
        for col in 0..n {
            b[col] -= j[i][col] * f[i];
            for k in 0..n {
                a[col][k] += j[i][col] * j[i][k];
            }
        }
    }
    // Tikhonov regularisation
    for i in 0..n {
        a[i][i] += lambda;
    }

    gaussian_elimination(&mut a, &mut b).unwrap_or_else(|| vec![0.0; n])
}

/// Gaussian elimination with partial pivoting; `None` when a pivot falls
/// below [`PIVOT_TOL`] (singular system).
pub(crate) fn gaussian_elimination(a: &mut Vec<Vec<f64>>, b: &mut Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    debug_assert_eq!(a.len(), n);

    for col in 0..n {
        // Find pivot row (max absolute value in column)
        let pivot_row =
            (col..n).max_by(|&i, &j| a[i][col].abs().partial_cmp(&a[j][col].abs()).unwrap())?;
        a.swap(col, pivot_row);
        b.swap(col, pivot_row);

        let pivot = a[col][col];
        if pivot.abs() < PIVOT_TOL {
            return None;
        }

        // Eliminate below
        for row in (col + 1)..n {
            let factor = a[row][col] / pivot;
            for k in col..n {
                let v = a[col][k] * factor;
                a[row][k] -= v;
            }
            let bv = b[col] * factor;
            b[row] -= bv;
        }
    }

    // Back substitution
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        x[i] = b[i];
        for j in (i + 1)..n {
            let v = a[i][j] * x[j];
            x[i] -= v;
        }
        x[i] /= a[i][i];
    }
    Some(x)
}

/// Root-mean-square of a residual vector.
pub(crate) fn rms(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gaussian_2x2() {
        // 2x + y = 5, x + 3y = 10  →  x=1, y=3
        let mut a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let mut b = vec![5.0, 10.0];
        let x = gaussian_elimination(&mut a, &mut b).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-10, "x={}", x[0]);
        assert!((x[1] - 3.0).abs() < 1e-10, "y={}", x[1]);
    }

    #[test]
    fn gaussian_3x3() {
        // x + y + z = 6, 2x + y = 5, y + 3z = 10  →  x=2, y=1, z=3
        let mut a = vec![
            vec![1.0, 1.0, 1.0],
            vec![2.0, 1.0, 0.0],
            vec![0.0, 1.0, 3.0],
        ];
        let mut b = vec![6.0, 5.0, 10.0];
        let x = gaussian_elimination(&mut a, &mut b).unwrap();
        assert!((x[0] - 2.0).abs() < 1e-9, "x={}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-9, "y={}", x[1]);
        assert!((x[2] - 3.0).abs() < 1e-9, "z={}", x[2]);
    }

    #[test]
    fn newton_core_solves_scalar_system() {
        // One parameter, one equation: x² - 4 = 0 from x = 5 → x = 2.
        let mut params = vec![5.0_f64];
        let result = newton_solve(&mut params, &[false], 1, |p, f| {
            f[0] = p[0] * p[0] - 4.0;
        });
        assert!(result.converged, "residual {}", result.residual);
        assert!((params[0] - 2.0).abs() < 1e-7, "x={}", params[0]);
    }

    #[test]
    fn newton_core_respects_fixed_mask() {
        let mut params = vec![5.0_f64, 9.0];
        let result = newton_solve(&mut params, &[true, false], 1, |p, f| {
            f[0] = p[1] - 1.0;
        });
        assert!(result.converged);
        assert_eq!(params[0], 5.0, "fixed param untouched");
        assert!((params[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn polish_nails_exactly_representable_fixed_point_bitwise() {
        // x² - 4 = 0 from x = 5. The loop terminates on the RMS threshold,
        // so the raw terminal iterate may sit a rounding step off the root;
        // the convergence polish must pin the exactly-representable fixed
        // point to machine zero — BITWISE 2.0, not merely within tolerance.
        let mut params = vec![5.0_f64];
        let result = newton_solve(&mut params, &[false], 1, |p, f| {
            f[0] = p[0] * p[0] - 4.0;
        });
        assert!(result.converged, "residual {}", result.residual);
        assert_eq!(
            result.residual, 0.0,
            "polished residual must be machine zero"
        );
        assert_eq!(params[0], 2.0, "polish must land on the exact root");
    }
}
