//! The family-independent pipeline.
//!
//! The 2D-sketch and 3D-geometric families differ in their entity and
//! constraint vocabularies and in nothing else. Everything downstream of
//! "constraint `i` has these residuals and these Jacobian rows" is shared
//! and lives here, written once against the [`Equations`] trait:
//!
//! - [`components`] — connected components of the constraint graph;
//! - [`solve_components`] — the per-component Newton solve;
//! - [`diagnose`] — the rank-revealing diagnosis;
//! - the helpers they share ([`eval_residuals`], [`analytic_jacobian`],
//!   [`touched_params`]).
//!
//! A family plugs in with an [`Equations`] adapter over its entity and
//! constraint lists, plus a [`Gauge`] if some of its storage parameters are
//! not physical freedoms (the 3D planes and cylinders).

use std::collections::HashMap;

use crate::diagnosis::{Diagnosis, CONSIST_TOL, PARTIAL_TOL, RANK_REL_TOL};
use crate::numerics::{
    self, gt, householder_qr_column_pivoting, le, solve_normal_equations, AdditiveSpace,
    CoreOptions, RowSpaceGram, SolveResult, SolveStatus, StepSpace, LAMBDA,
};

/// A sparse row or direction: `(parameter or column index, value)` pairs.
pub(crate) type SparseRow = Vec<(usize, f64)>;

/// A list of constraints, each contributing scalar equations over a flat
/// parameter vector.
pub(crate) trait Equations {
    /// Number of constraints.
    fn constraint_count(&self) -> usize;

    /// Number of scalar equations constraint `index` contributes.
    fn equation_count(&self, index: usize) -> usize;

    /// Residuals of constraint `index` at `params`, written to `out`
    /// (`equation_count(index)` long).
    fn residuals(&self, index: usize, params: &[f64], out: &mut [f64]);

    /// Jacobian rows of constraint `index` at `params`: `out[r]` receives
    /// the `(column, derivative)` pairs of equation `r`. Every column the
    /// constraint touches is listed, whatever the derivative value — the
    /// decomposition reads the touched set from here.
    fn jacobian(&self, index: usize, params: &[f64], out: &mut [SparseRow]);

    /// Number of columns: the coordinates that Jacobian entries, the fixed
    /// mask and the degrees of freedom are counted in. They are the stored
    /// parameters themselves unless the family steps on a manifold, where a
    /// rigid body has seven stored numbers and six columns.
    fn column_count(&self, params: &[f64]) -> usize {
        params.len()
    }

    /// The step space of the free columns `free` (sorted ascending): how a
    /// step in those columns moves the stored parameters. Additive by
    /// default — column `i` is parameter `i`.
    fn step_space(&self, free: &[usize]) -> Box<dyn StepSpace + '_> {
        Box::new(AdditiveSpace::new(free.to_vec()))
    }

    /// Whether a solve cross-checks the Jacobian against finite differences
    /// in debug builds (two residual evaluations per free column).
    fn debug_check_jacobian(&self) -> bool {
        true
    }
}

/// The indices of every constraint, in order.
pub(crate) fn all_indices<E: Equations + ?Sized>(eqs: &E) -> Vec<usize> {
    (0..eqs.constraint_count()).collect()
}

/// Number of scalar equations of the constraints `indices`.
pub(crate) fn equation_total<E: Equations + ?Sized>(eqs: &E, indices: &[usize]) -> usize {
    indices.iter().map(|&i| eqs.equation_count(i)).sum()
}

/// The stacked residual vector of the constraints `indices` at `params`.
pub(crate) fn eval_residuals<E: Equations + ?Sized>(
    eqs: &E,
    indices: &[usize],
    params: &[f64],
) -> Vec<f64> {
    let mut f = vec![0.0_f64; equation_total(eqs, indices)];
    let mut row = 0;
    for &i in indices {
        let count = eqs.equation_count(i);
        eqs.residuals(i, params, &mut f[row..row + count]);
        row += count;
    }
    f
}

/// Assemble the dense Jacobian (equations of `indices` × the columns of
/// `col_of`) from the per-constraint sparse rows. Parameters without a
/// column (those held by the fixed mask) are dropped. `out` arrives zeroed.
pub(crate) fn analytic_jacobian<E: Equations + ?Sized>(
    eqs: &E,
    indices: &[usize],
    params: &[f64],
    col_of: &HashMap<usize, usize>,
    out: &mut [Vec<f64>],
) {
    let mut rows: Vec<SparseRow> = Vec::new();
    let mut row = 0;
    for &i in indices {
        let count = eqs.equation_count(i);
        rows.clear();
        rows.resize_with(count, Vec::new);
        eqs.jacobian(i, params, &mut rows);
        for (k, entries) in rows.iter().enumerate() {
            for &(param, v) in entries {
                if let Some(&col) = col_of.get(&param) {
                    out[row + k][col] += v;
                }
            }
        }
        row += count;
    }
}

/// Sorted, de-duplicated set of the parameters the constraints `indices`
/// touch (fixed or not).
pub(crate) fn touched_params<E: Equations + ?Sized>(
    eqs: &E,
    indices: &[usize],
    params: &[f64],
) -> Vec<usize> {
    let mut touched: Vec<usize> = Vec::new();
    let mut scratch: Vec<SparseRow> = Vec::new();
    for &i in indices {
        scratch.clear();
        scratch.resize_with(eqs.equation_count(i), Vec::new);
        eqs.jacobian(i, params, &mut scratch);
        for row in &scratch {
            touched.extend(row.iter().map(|&(p, _)| p));
        }
    }
    touched.sort_unstable();
    touched.dedup();
    touched
}

/// Connected components of the constraint graph.
///
/// Constraints and parameters form a bipartite graph (a constraint is
/// adjacent to every parameter it touches); its components partition the
/// constraints into subsystems that share no parameter and can be solved,
/// diagnosed and skipped independently. Returns groups of constraint
/// indices, ordered by first appearance, so the partition is deterministic.
/// A constraint that touches no parameter forms a group of its own.
pub(crate) fn components<E: Equations + ?Sized>(eqs: &E, params: &[f64]) -> Vec<Vec<usize>> {
    components_of(eqs, &all_indices(eqs), params)
}

/// [`components`] of the constraints `indices` alone: the connected
/// components of the graph those constraints span, as groups of the given
/// indices in order of first appearance.
pub(crate) fn components_of<E: Equations + ?Sized>(
    eqs: &E,
    indices: &[usize],
    params: &[f64],
) -> Vec<Vec<usize>> {
    // Union-find over columns.
    let mut parent: Vec<usize> = (0..eqs.column_count(params)).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        let mut root = x;
        while parent[root] != root {
            root = parent[root];
        }
        let mut cur = x;
        while parent[cur] != root {
            let next = parent[cur];
            parent[cur] = root;
            cur = next;
        }
        root
    }
    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[ra] = rb;
        }
    }

    // Pass 1: every constraint unions the parameters it touches. Its first
    // touched parameter is remembered for the grouping pass.
    let mut first_param: Vec<Option<usize>> = Vec::with_capacity(indices.len());
    let mut scratch: Vec<SparseRow> = Vec::new();
    for &i in indices {
        scratch.clear();
        scratch.resize_with(eqs.equation_count(i), Vec::new);
        eqs.jacobian(i, params, &mut scratch);
        let mut first: Option<usize> = None;
        for row in &scratch {
            for &(p, _) in row {
                if let Some(f) = first {
                    union(&mut parent, f, p);
                } else {
                    first = Some(p);
                }
            }
        }
        first_param.push(first);
    }

    // Pass 2: group constraints by the FINAL root of their first touched
    // parameter. Grouping during pass 1 would key earlier constraints by
    // stale roots when a later constraint merges two established groups
    // (regression-tested).
    let mut order: Vec<usize> = Vec::new();
    let mut by_root: HashMap<usize, Vec<usize>> = HashMap::new();
    for (&index, first) in indices.iter().zip(first_param) {
        // A constraint with no equations (suppressed or removed) is not part
        // of the system: it gets no group.
        if eqs.equation_count(index) == 0 {
            continue;
        }
        match first {
            Some(param) => {
                let root = find(&mut parent, param);
                by_root
                    .entry(root)
                    .or_insert_with(|| {
                        order.push(root);
                        Vec::new()
                    })
                    .push(index);
            }
            // A constraint touching no parameters has no equations to
            // satisfy through the solver; it forms its own (trivial) group
            // so indices are never dropped.
            None => {
                order.push(usize::MAX - index);
                by_root.insert(usize::MAX - index, vec![index]);
            }
        }
    }
    order
        .into_iter()
        .filter_map(|root| by_root.remove(&root))
        .collect()
}

/// Solve the selected connected components (`only` holds ordinals into
/// `components`; `None` selects all), one Newton solve each. Unselected
/// components are not iterated and their parameters are not touched.
///
/// `after_solve` runs once the selected components are done, before the
/// final residual is taken; it receives the parameters and the sorted set
/// of parameters the solved components touch (the 3D family re-pins its
/// plane representation there).
///
/// The reported status is the most severe one among the solved components;
/// the failure policy of `options` applies to each component separately.
pub(crate) fn solve_components<E, F>(
    eqs: &E,
    params: &mut [f64],
    fixed: &[bool],
    components: &[Vec<usize>],
    only: Option<&[usize]>,
    options: &CoreOptions,
    after_solve: F,
) -> SolveResult
where
    E: Equations + ?Sized,
    F: FnOnce(&mut [f64], &[usize]),
{
    let mut iterations = 0usize;
    let mut converged = true;
    let mut status = SolveStatus::Converged;
    let mut solved_touched: Vec<usize> = Vec::new();

    for (ordinal, group) in components.iter().enumerate() {
        if let Some(selected) = only {
            if !selected.contains(&ordinal) {
                continue;
            }
        }
        let (result, touched) = solve_one(eqs, group, params, fixed, options);
        iterations += result.iterations;
        converged &= result.converged;
        status = status.max(result.status);
        // Components share no parameters, so the union stays duplicate-free.
        solved_touched.extend(touched);
    }
    solved_touched.sort_unstable();
    after_solve(params, &solved_touched);

    let f = eval_residuals(eqs, &all_indices(eqs), params);
    SolveResult {
        converged,
        residual: numerics::rms(&f),
        iterations,
        status,
    }
}

/// Solve one connected subsystem: only the free parameters its constraints
/// touch may move. Returns the result and the touched parameter set.
fn solve_one<E: Equations + ?Sized>(
    eqs: &E,
    group: &[usize],
    params: &mut [f64],
    fixed: &[bool],
    options: &CoreOptions,
) -> (SolveResult, Vec<usize>) {
    let touched = touched_params(eqs, group, params);
    // Only this component's own free columns move: the Newton loop's free
    // list is exactly its movable coordinates.
    let free: Vec<usize> = touched.iter().copied().filter(|&p| !fixed[p]).collect();
    let col_of: HashMap<usize, usize> = free.iter().enumerate().map(|(c, &p)| (p, c)).collect();
    let space = eqs.step_space(&free);

    let n_eq = equation_total(eqs, group);
    let mut build = |p: &[f64], out: &mut [Vec<f64>]| {
        analytic_jacobian(eqs, group, p, &col_of, out);
    };

    // One-shot integration check (debug builds): the analytic assembly
    // must agree with central differences at the starting point. The
    // exhaustive per-variant gates live with the constraint definitions;
    // this catches per-call assembly bugs (column mapping, mask filtering).
    #[cfg(debug_assertions)]
    if n_eq > 0 && !free.is_empty() && eqs.debug_check_jacobian() {
        let mut eval = |p: &[f64], f: &mut [f64]| {
            f.copy_from_slice(&eval_residuals(eqs, group, p));
        };
        let mut j_analytic = vec![vec![0.0_f64; free.len()]; n_eq];
        build(params, &mut j_analytic);
        let mut j_fd = vec![vec![0.0_f64; free.len()]; n_eq];
        numerics::numerical_jacobian_for(params, &*space, n_eq, &mut eval, &mut j_fd);
        for (ri, (ra, rf)) in j_analytic.iter().zip(j_fd.iter()).enumerate() {
            for (ci, (va, vf)) in ra.iter().zip(rf.iter()).enumerate() {
                assert!(
                    (va - vf).abs() <= 1.0e-4 * (1.0 + vf.abs()),
                    "analytic Jacobian drifts from finite differences at row {ri} col {ci}: {va} vs {vf}"
                );
            }
        }
    }

    let result = numerics::newton_solve_in(
        params,
        &*space,
        n_eq,
        |p, f| {
            f.copy_from_slice(&eval_residuals(eqs, group, p));
        },
        numerics::Jacobian::Analytic(&mut build),
        options,
    );
    (result, touched)
}

/// Storage parameters that are not physical freedoms.
///
/// A family whose entities are stored with more parameters than they have
/// freedoms describes the surplus here, so the diagnosis does not report it
/// as degrees of freedom. The default is "none" (see [`NoGauge`]).
pub(crate) trait Gauge {
    /// Number of gauge dimensions that [`Gauge::effective_columns`] removes
    /// outright across the whole system.
    fn removed_dimensions(&self, _params: &[f64], _fixed: &[bool]) -> usize {
        0
    }

    /// The effective free directions of one component, as sparse
    /// combinations of its local columns (`touched[k]` is local column
    /// `k`). `None` means the columns are used as they are.
    fn effective_columns(
        &self,
        _params: &[f64],
        _fixed: &[bool],
        _touched: &[usize],
    ) -> Option<Vec<SparseRow>> {
        None
    }

    /// Directions (unit length, in global parameter indices) that are
    /// representation freedom unless a constraint row pins them. Each one
    /// that stays unpinned is subtracted from the degrees of freedom.
    fn soft_directions(&self, _params: &[f64], _fixed: &[bool]) -> Vec<SparseRow> {
        Vec::new()
    }
}

/// Every storage parameter is a physical freedom.
pub(crate) struct NoGauge;

impl Gauge for NoGauge {}

/// Rank-revealing diagnosis of the constraints at `params`.
///
/// Per connected component, the transposed Jacobian over the free
/// parameters is factored by a column-pivoted Householder QR: the pivot
/// columns are the independent constraint rows, and each remaining
/// (dependent) row is classified by whether its residual agrees with the
/// least-squares step of the independent rows — redundant, partially
/// redundant, or conflicting. Dependencies can only occur inside a
/// component, so the per-component analysis is equivalent to a global one.
/// See [`crate::diagnosis`] for the vocabulary.
pub(crate) fn diagnose<E, G>(eqs: &E, params: &[f64], fixed: &[bool], gauge: &G) -> Diagnosis
where
    E: Equations + ?Sized,
    G: Gauge + ?Sized,
{
    // Free columns: the fixed mask is indexed by column.
    let free_total = fixed.iter().filter(|held| !**held).count();
    let all = all_indices(eqs);
    let n_eq = equation_total(eqs, &all);

    let removed = gauge.removed_dimensions(params, fixed);
    // Soft gauge directions; how many of them stay unpinned is settled per
    // component below.
    let soft = gauge.soft_directions(params, fixed);
    let mut soft_seen = vec![false; soft.len()];
    let mut soft_free = 0usize;

    if n_eq == 0 {
        return Diagnosis {
            dof: free_total as i64 - removed as i64 - soft.len() as i64,
            ..Default::default()
        };
    }

    // Global residuals once; per-component slices below.
    let f = eval_residuals(eqs, &all, params);
    // Row offset of each constraint in the stacked residual vector.
    let mut row_of: Vec<usize> = Vec::with_capacity(all.len());
    let mut acc = 0usize;
    for &i in &all {
        row_of.push(acc);
        acc += eqs.equation_count(i);
    }

    let mut d = Diagnosis::default();
    // True DOF = free_total − Σ rank_c: parameters claimed by no component
    // are unconstrained and count fully (e.g. an entity whose only
    // constraint pins one endpoint).
    let mut total_rank = 0i64;

    for group in &components(eqs, params) {
        let n_eq_c = equation_total(eqs, group);
        // Local residual slice, in component row order.
        let mut f_c = Vec::with_capacity(n_eq_c);
        for &gi in group {
            f_c.extend_from_slice(&f[row_of[gi]..row_of[gi] + eqs.equation_count(gi)]);
        }

        // Touched free parameters of this component.
        let mut touched = touched_params(eqs, group, params);
        touched.retain(|&p| !fixed[p]);

        // Map component-local rows back to global constraint indices.
        let mut row_constraint: Vec<usize> = Vec::with_capacity(n_eq_c);
        for &gi in group {
            for _ in 0..eqs.equation_count(gi) {
                row_constraint.push(gi);
            }
        }

        // Nothing in this component can move: any unsatisfied residual is a
        // genuine conflict.
        let flag_unsatisfied = |d: &mut Diagnosis, f_c: &[f64]| {
            let scale = f_c.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            for (r, fi) in f_c.iter().enumerate() {
                if fi.abs() > CONSIST_TOL * scale {
                    d.conflicting.push(row_constraint[r]);
                }
            }
        };
        if touched.is_empty() {
            flag_unsatisfied(&mut d, &f_c);
            continue;
        }

        // Raw Jacobian over the touched local columns.
        let col_of: HashMap<usize, usize> =
            touched.iter().enumerate().map(|(c, &p)| (p, c)).collect();
        let mut j_raw = vec![vec![0.0_f64; touched.len()]; n_eq_c];
        analytic_jacobian(eqs, group, params, &col_of, &mut j_raw);

        // Soft gauge directions this component touches but does not pin.
        soft_free += unpinned_gauge_dims(&j_raw, &col_of, &soft, &mut soft_seen);

        // Effective free directions: hard gauge is removed here.
        let mut j = match gauge.effective_columns(params, fixed, &touched) {
            None => j_raw,
            Some(columns) => {
                let mut j = vec![vec![0.0_f64; columns.len()]; n_eq_c];
                for (row, jr) in j_raw.iter().enumerate() {
                    for (c, col) in columns.iter().enumerate() {
                        let mut v = 0.0_f64;
                        for &(lc, coeff) in col {
                            v += jr[lc] * coeff;
                        }
                        j[row][c] = v;
                    }
                }
                j
            }
        };
        let n_cols = j.first().map_or(0, |row| row.len());
        if n_cols == 0 {
            flag_unsatisfied(&mut d, &f_c);
            continue;
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

        // A = Jᵀ: column i of A is component row i (equilibrated).
        let mut a = vec![vec![0.0_f64; n_eq_c]; n_cols];
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

        // Dependency groups: the dependency lives in ROW space — find
        // coefficients y with Σ y_j·row_j = dep_row via the row-space Gram
        // system (G + λI) y = rows·dep, G[i][j] = row_i·row_j. The
        // coefficient support is the whole mutually-dependent row group,
        // invariant to pivot tie-breaking. One Gram matrix serves every
        // dependent row: it is built and factored once.
        let rows: Vec<Vec<f64>> = independent
            .iter()
            .map(|&col| (0..n_cols).map(|rr| a_original[rr][col]).collect())
            .collect();
        let r = rows.len();
        let gram = RowSpaceGram::new(&rows, n_cols, LAMBDA);
        for &dep in &dependent {
            let dep_row: Vec<f64> = (0..n_cols).map(|rr| a_original[rr][dep]).collect();
            let y = gram.coefficients(&dep_row);
            let max_c = y.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            let support: Vec<usize> = (0..r).filter(|&jj| y[jj].abs() > 1.0e-6 * max_c).collect();
            // Self-verification: the reported group must be a REAL relation
            // — reconstruct the dependent row from the support and reject
            // the group when the reconstruction residual is not at roundoff
            // (an ill-conditioned basis must never surface a fake group).
            let residual: f64 = (0..n_cols)
                .map(|rr| {
                    let recon: f64 = support.iter().map(|&jj| y[jj] * rows[jj][rr]).sum();
                    (recon - dep_row[rr]).abs()
                })
                .fold(0.0_f64, f64::max);
            let scale = dep_row.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            if !le(residual, 1.0e-9 * scale) {
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
        let z = solve_normal_equations(&j_sub, &f_sub, n_cols);

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

    // Soft gauge directions no component touches are unpinned by definition.
    soft_free += soft_seen.iter().filter(|seen| !**seen).count();
    d.rank = total_rank.max(0) as usize;
    d.dof = free_total as i64 - removed as i64 - soft_free as i64 - total_rank;
    d
}

/// Number of soft gauge directions that this component touches and that no
/// row of its Jacobian pins — freedom of the representation, not of the
/// geometry.
///
/// `j_raw` is the component Jacobian over its touched columns (`col_of`
/// maps global parameter index to local column). `seen` is flagged for
/// every direction the component touches, so the caller can count the
/// directions no component touches at all. The count is `k − rank(J·G)` for
/// the `k` touched directions `G`: the dimension of their span inside the
/// null space.
fn unpinned_gauge_dims(
    j_raw: &[Vec<f64>],
    col_of: &HashMap<usize, usize>,
    gauges: &[SparseRow],
    seen: &mut [bool],
) -> usize {
    let mut local: Vec<SparseRow> = Vec::new();
    for (gauge, flag) in gauges.iter().zip(seen.iter_mut()) {
        let cols: SparseRow = gauge
            .iter()
            .filter_map(|&(p, coeff)| col_of.get(&p).map(|&lc| (lc, coeff)))
            .collect();
        if !cols.is_empty() {
            *flag = true;
            local.push(cols);
        }
    }
    if local.is_empty() {
        return 0;
    }

    // M = J·G, rows scaled to a unit largest Jacobian entry like the rank
    // analysis. Against such a row a unit gauge direction is O(1) when the
    // row pins it and roundoff when it does not.
    let mut m = vec![vec![0.0_f64; local.len()]; j_raw.len()];
    for (m_row, j_row) in m.iter_mut().zip(j_raw) {
        let scale = j_row.iter().fold(0.0_f64, |mx, v| mx.max(v.abs()));
        if !gt(scale, 0.0) {
            continue;
        }
        for (entry, cols) in m_row.iter_mut().zip(&local) {
            let v = cols.iter().map(|&(lc, c)| j_row[lc] * c).sum::<f64>() / scale;
            if v.abs() > RANK_REL_TOL {
                *entry = v;
            }
        }
    }
    let (_, rank) = householder_qr_column_pivoting(&mut m, RANK_REL_TOL);
    local.len() - rank
}
