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
use std::sync::OnceLock;

use crate::alias::{Aliases, ClassSpace, Classes};
use crate::diagnosis::{Diagnosis, CONSIST_TOL, PARTIAL_TOL, RANK_REL_TOL};
use crate::numerics::{
    self, gt, householder_qr_column_pivoting, le, solve_normal_equations, AdditiveSpace,
    CoreOptions, RowSpaceGram, SolveResult, SolveStatus, SparseJacobian, StepSpace, LAMBDA,
};
use crate::structure;

/// A sparse row or direction: `(parameter or column index, value)` pairs.
pub(crate) type SparseRow = Vec<(usize, f64)>;

/// An equation whose residual is linear in at most two parameters with
/// unit coefficients (see [`Equations::alias`]).
#[derive(Debug, Clone, Copy)]
pub(crate) enum Alias {
    /// The residual is `params[param] − params[base] − offset`.
    Offset {
        param: usize,
        base: usize,
        offset: f64,
    },
    /// The residual is `params[param] − value`.
    Value { param: usize, value: f64 },
}

/// `Sync` when the connected components of one system may be solved on
/// several threads (feature `rayon`), nothing otherwise: a family that
/// caches per constraint then pays for a lock only where one is needed.
#[cfg(feature = "rayon")]
pub(crate) trait SharedAcrossThreads: Sync {}
#[cfg(feature = "rayon")]
impl<T: Sync + ?Sized> SharedAcrossThreads for T {}
/// See the `rayon` version.
#[cfg(not(feature = "rayon"))]
pub(crate) trait SharedAcrossThreads {}
#[cfg(not(feature = "rayon"))]
impl<T: ?Sized> SharedAcrossThreads for T {}

/// A list of constraints, each contributing scalar equations over a flat
/// parameter vector.
///
/// With feature `rayon` it is `Sync` (see [`SharedAcrossThreads`]): a
/// component evaluates only its own constraints, so whatever an
/// implementation caches per constraint is never contended.
pub(crate) trait Equations: SharedAcrossThreads {
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

    /// Whether [`Equations::step_space`] is the additive one, column `i`
    /// being parameter `i`: then a component can be solved in parts (see
    /// [`crate::structure`]), each over some of its parameters.
    fn steps_are_additive(&self) -> bool {
        true
    }

    /// The size of the quantity that column `column` steps, for judging
    /// what a rounding of it is worth: the magnitude of the parameter by
    /// default.
    fn column_magnitude(&self, column: usize, params: &[f64]) -> f64 {
        params[column].abs()
    }

    /// The size of the geometry the constraints `group` act on, if every
    /// equation of the family is a length: a component smaller than one
    /// unit is then held to the tolerance in proportion to that size,
    /// rather than to an absolute number that would let its equations stop
    /// far short of holding relative to it. `None` (the default) keeps the
    /// tolerance absolute.
    fn length_scale(&self, _group: &[usize], _params: &[f64]) -> Option<f64> {
        None
    }

    /// Equation `row` of constraint `index`, if its residual is exactly
    /// one parameter less another and a constant, or one parameter less a
    /// constant (see [`Alias`]): such an equation is eliminated before the
    /// solve rather than iterated on. `None` (the default) for any other.
    /// Only a family whose steps are additive may return one.
    fn alias(&self, _index: usize, _row: usize) -> Option<Alias> {
        None
    }

    /// Whether a solve cross-checks the Jacobian against finite differences
    /// in debug builds (two residual evaluations per free column).
    #[cfg_attr(not(debug_assertions), allow(dead_code))]
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
    eval_residuals_into(eqs, indices, params, &mut f);
    f
}

/// [`eval_residuals`] into the caller's buffer, which has one entry per
/// equation of `indices`.
pub(crate) fn eval_residuals_into<E: Equations + ?Sized>(
    eqs: &E,
    indices: &[usize],
    params: &[f64],
    f: &mut [f64],
) {
    let mut row = 0;
    for &i in indices {
        let count = eqs.equation_count(i);
        eqs.residuals(i, params, &mut f[row..row + count]);
        row += count;
    }
    debug_assert_eq!(row, f.len());
}

/// Sentinel of [`Assembly`]: an entry whose column is held, which the
/// Jacobian has no column for.
const HELD: usize = usize::MAX;

/// How the Jacobian of one component is assembled from the per-constraint
/// rows of [`Equations::jacobian`], laid out once per solve.
///
/// Each equation's row of the [`SparseJacobian`] holds the free columns its
/// constraint lists, and every entry the constraint lists is mapped to its
/// slot in the values when the layout is made. A build then writes each
/// derivative straight to its slot, with no lookup of the column. (The
/// trait promises that a constraint lists the same columns whatever their
/// values, in a fixed order; an entry that arrives in another order is
/// still placed, by a search of its row.)
struct Assembly {
    /// The free ones among them, ascending: column `k` of the Jacobian is
    /// `free[k]`.
    free: Vec<usize>,
    /// Per column of the family, its Jacobian column, or [`HELD`].
    col_of: Vec<usize>,
    /// The column `col_of[0]` is for: the lowest one touched.
    col_base: usize,
    /// Per listed row, its row of the pattern, or [`HELD`] for a row the
    /// [`Shape`] leaves out.
    pattern_row: Vec<usize>,
    /// Row `r` lists entries `entry_starts[r]..entry_starts[r + 1]`.
    entry_starts: Vec<usize>,
    /// The column of each listed entry, as the layout saw it.
    listed: Vec<usize>,
    /// The slot of each listed entry in the values, or [`HELD`].
    slots: Vec<usize>,
}

/// What an [`Assembly`] lays out, beyond one column per free parameter
/// and one row per equation.
#[derive(Clone, Copy, Default)]
struct Shape<'a> {
    /// One column per unknown class of these, each member's entries summed
    /// into it; the parameters of a known class are held.
    classes: Option<&'a Classes>,
    /// Rows laid out with no entries (they still have residuals): the
    /// absorbed ones.
    empty: Option<&'a [bool]>,
    /// The rows laid out at all; the others are not part of the Jacobian.
    keep: Option<&'a [bool]>,
    /// The only columns that may be free, by key (the class root, or the
    /// parameter), ascending; any other is held.
    columns: Option<&'a [usize]>,
}

/// The columns each equation of a component lists, the first pass of an
/// [`Assembly`].
struct Listing {
    /// Row `r` lists entries `entry_starts[r]..entry_starts[r + 1]`.
    entry_starts: Vec<usize>,
    /// The column of each listed entry.
    listed: Vec<usize>,
    /// The columns listed, ascending, each once.
    touched: Vec<usize>,
}

impl Listing {
    fn new<E: Equations + ?Sized>(
        eqs: &E,
        indices: &[usize],
        params: &[f64],
        scratch: &mut Vec<SparseRow>,
    ) -> Self {
        let mut entry_starts = vec![0];
        let mut listed: Vec<usize> = Vec::new();
        for &i in indices {
            for row in jacobian_rows(eqs, i, params, scratch) {
                listed.extend(row.iter().map(|&(column, _)| column));
                entry_starts.push(listed.len());
            }
        }
        let mut touched = listed.clone();
        touched.sort_unstable();
        touched.dedup();
        Self {
            entry_starts,
            listed,
            touched,
        }
    }
}

impl Assembly {
    /// Lay out the Jacobian of `listing` over the columns `fixed` does not
    /// hold, in the [`Shape`] `shape`.
    fn lay_out(listing: Listing, fixed: &[bool], shape: Shape<'_>) -> (Self, SparseJacobian) {
        let Listing {
            entry_starts,
            listed,
            touched,
        } = listing;
        let n_listed_rows = entry_starts.len() - 1;
        let empty = |r: usize| shape.empty.is_some_and(|empty| empty[r]);
        let kept = |r: usize| shape.keep.is_none_or(|keep| keep[r]);
        // The key of a free column: its class root, or itself.
        let key = |c: usize| match shape.classes {
            Some(classes) if classes.is_known(c) => None,
            Some(classes) => Some(classes.root_of(c)),
            None => Some(c),
        };
        let allowed = |k: usize| {
            shape
                .columns
                .is_none_or(|cols| cols.binary_search(&k).is_ok())
        };
        let free_key = |c: usize| {
            if fixed[c] {
                return None;
            }
            key(c).filter(|&k| allowed(k))
        };

        // Over the span of the touched columns only: a component of a large
        // system touches a few columns far from zero.
        let col_base = touched.first().copied().unwrap_or(0);
        let mut col_of = vec![HELD; touched.last().map_or(0, |&c| c + 1 - col_base)];
        let free: Vec<usize> = if shape.classes.is_none() && shape.columns.is_none() {
            let free: Vec<usize> = touched.iter().copied().filter(|&c| !fixed[c]).collect();
            for (k, &c) in free.iter().enumerate() {
                col_of[c - col_base] = k;
            }
            free
        } else {
            // The keys that a laid-out row touches.
            let mut keys: Vec<usize> = Vec::new();
            for (r, span) in entry_starts.windows(2).enumerate() {
                if empty(r) || !kept(r) {
                    continue;
                }
                keys.extend(listed[span[0]..span[1]].iter().filter_map(|&c| free_key(c)));
            }
            keys.sort_unstable();
            keys.dedup();
            for &c in &touched {
                if let Some(k) = free_key(c) {
                    if let Ok(at) = keys.binary_search(&k) {
                        col_of[c - col_base] = at;
                    }
                }
            }
            keys
        };

        // Per laid-out row, its free columns ascending and each once,
        // written straight into the pattern's column list.
        let mut pattern_row = vec![HELD; n_listed_rows];
        let mut starts = Vec::with_capacity(entry_starts.len());
        let mut cols: Vec<usize> = Vec::with_capacity(listed.len());
        starts.push(0);
        for (r, span) in entry_starts.windows(2).enumerate() {
            if !kept(r) {
                continue;
            }
            pattern_row[r] = starts.len() - 1;
            let begin = cols.len();
            if !empty(r) {
                cols.extend(
                    listed[span[0]..span[1]]
                        .iter()
                        .map(|&c| col_of[c - col_base])
                        .filter(|&k| k != HELD),
                );
            }
            cols[begin..].sort_unstable();
            let mut unique = begin;
            for at in begin..cols.len() {
                if unique == begin || cols[at] != cols[unique - 1] {
                    cols[unique] = cols[at];
                    unique += 1;
                }
            }
            cols.truncate(unique);
            starts.push(cols.len());
        }
        let pattern = SparseJacobian::from_parts(free.len(), starts, cols);
        let mut assembly = Self {
            free,
            col_of,
            col_base,
            entry_starts,
            listed,
            pattern_row,
            slots: Vec::new(),
        };
        let mut slots = Vec::with_capacity(assembly.listed.len());
        for r in 0..n_listed_rows {
            for e in assembly.entry_starts[r]..assembly.entry_starts[r + 1] {
                slots.push(if empty(r) {
                    HELD
                } else {
                    assembly.slot(&pattern, r, assembly.listed[e])
                });
            }
        }
        assembly.slots = slots;
        (assembly, pattern)
    }

    /// The slot of column `column` in listed row `r`, or [`HELD`].
    fn slot(&self, pattern: &SparseJacobian, r: usize, column: usize) -> usize {
        let row = self.pattern_row[r];
        let k = column
            .checked_sub(self.col_base)
            .and_then(|at| self.col_of.get(at))
            .copied()
            .unwrap_or(HELD);
        if k == HELD || row == HELD {
            return HELD;
        }
        match pattern.row_cols(row).binary_search(&k) {
            Ok(at) => pattern.row_start(row) + at,
            Err(_) => {
                // Outside the layout: the constraint listed a column it did
                // not list when the layout was made (or the row was laid
                // out empty).
                debug_assert!(
                    pattern.row_cols(row).is_empty(),
                    "row {r} lists column {column} outside its pattern"
                );
                HELD
            }
        }
    }

    /// Write the Jacobian of the constraints `indices` (those it was laid
    /// out for) at `params` into the values of `j`, which arrive zeroed.
    /// `scratch` is reused between builds.
    fn assemble<E: Equations + ?Sized>(
        &self,
        eqs: &E,
        indices: &[usize],
        params: &[f64],
        j: &mut SparseJacobian,
        scratch: &mut Vec<SparseRow>,
    ) {
        let mut r = 0;
        for &i in indices {
            for row in jacobian_rows(eqs, i, params, scratch) {
                let base = self.entry_starts[r];
                let laid_out = &self.listed[base..self.entry_starts[r + 1]];
                for (k, &(column, v)) in row.iter().enumerate() {
                    let slot = if laid_out.get(k) == Some(&column) {
                        self.slots[base + k]
                    } else {
                        self.slot(j, r, column)
                    };
                    if slot != HELD {
                        j.values_mut()[slot] += v;
                    }
                }
                r += 1;
            }
        }
    }
}

/// The Jacobian rows of constraint `index` at `params`, written into the
/// first rows of `scratch`, which keeps its rows (and their capacity) from
/// one call to the next.
fn jacobian_rows<'s, E: Equations + ?Sized>(
    eqs: &E,
    index: usize,
    params: &[f64],
    scratch: &'s mut Vec<SparseRow>,
) -> &'s [SparseRow] {
    let count = eqs.equation_count(index);
    if scratch.len() < count {
        scratch.resize_with(count, Vec::new);
    }
    let rows = &mut scratch[..count];
    for row in rows.iter_mut() {
        row.clear();
    }
    eqs.jacobian(index, params, rows);
    rows
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

/// Per equation of the constraints `indices`, the size of the numbers it is
/// computed from: `Σ|∂f/∂q|·|q|` over every column the equation touches,
/// held or free. One rounding of those numbers moves the residual by about
/// `ε` times this.
pub(crate) fn row_magnitudes<E: Equations + ?Sized>(
    eqs: &E,
    indices: &[usize],
    params: &[f64],
    out: &mut [f64],
) {
    let mut rows: Vec<SparseRow> = Vec::new();
    let mut row = 0;
    for &i in indices {
        let count = eqs.equation_count(i);
        rows.clear();
        rows.resize_with(count, Vec::new);
        eqs.jacobian(i, params, &mut rows);
        for (k, entries) in rows.iter().enumerate() {
            out[row + k] = entries
                .iter()
                .map(|&(column, v)| v.abs() * eqs.column_magnitude(column, params))
                .sum();
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
    solve_parts(
        eqs,
        params,
        fixed,
        components,
        None,
        only,
        options,
        after_solve,
    )
}

/// [`solve_components`] over the components of `plans`, keeping the plan of
/// each component it solves there for the next solve (see [`Plan`]).
pub(crate) fn solve_planned<E, F>(
    eqs: &E,
    params: &mut [f64],
    fixed: &[bool],
    plans: &Plans,
    only: Option<&[usize]>,
    options: &CoreOptions,
    after_solve: F,
) -> SolveResult
where
    E: Equations + ?Sized,
    F: FnOnce(&mut [f64], &[usize]),
{
    solve_parts(
        eqs,
        params,
        fixed,
        plans.components(),
        Some(plans),
        only,
        options,
        after_solve,
    )
}

#[allow(clippy::too_many_arguments)]
fn solve_parts<E, F>(
    eqs: &E,
    params: &mut [f64],
    fixed: &[bool],
    components: &[Vec<usize>],
    plans: Option<&Plans>,
    only: Option<&[usize]>,
    options: &CoreOptions,
    after_solve: F,
) -> SolveResult
where
    E: Equations + ?Sized,
    F: FnOnce(&mut [f64], &[usize]),
{
    let groups: Vec<&[usize]> = components
        .iter()
        .enumerate()
        .filter(|(ordinal, _)| only.is_none_or(|selected| selected.contains(ordinal)))
        .map(|(_, group)| group.as_slice())
        .collect();

    #[cfg(feature = "rayon")]
    let solved = {
        let runs = parallel_runs(eqs, &groups);
        if runs.len() >= 2 {
            solve_in_parallel(eqs, params, fixed, plans, &runs, options)
        } else {
            solve_in_turn(eqs, params, fixed, plans, &groups, options)
        }
    };
    #[cfg(not(feature = "rayon"))]
    let solved = solve_in_turn(eqs, params, fixed, plans, &groups, options);

    let mut iterations = 0usize;
    let mut converged = true;
    let mut status = SolveStatus::Converged;
    let mut solved_touched: Vec<usize> = Vec::new();
    for one in solved {
        iterations += one.result.iterations;
        converged &= one.result.converged;
        status = status.max(one.result.status);
        // Components share no parameters, so the union stays duplicate-free.
        solved_touched.extend(one.touched);
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

/// What the solve of one component left behind.
struct Solved {
    result: SolveResult,
    /// The columns its constraints touch, held or free.
    touched: Vec<usize>,
    /// The stored parameters the solve may have changed (those of its free
    /// columns; a rigid body stores seven numbers for six columns).
    #[cfg(feature = "rayon")]
    moved: Vec<usize>,
}

/// Solve the components `groups` one after the other on `params`.
fn solve_in_turn<E: Equations + ?Sized>(
    eqs: &E,
    params: &mut [f64],
    fixed: &[bool],
    plans: Option<&Plans>,
    groups: &[&[usize]],
    options: &CoreOptions,
) -> Vec<Solved> {
    groups
        .iter()
        .map(|group| solve_one(eqs, group, plans, params, fixed, options))
        .collect()
}

/// Solve runs of components on the threads of the current rayon pool, one
/// task per run (see [`parallel_runs`]).
///
/// Each task works on its own copy of the parameters and solves its run in
/// turn on it; afterwards every component's moved parameters are copied
/// back. A component reads only the parameters its constraints touch and
/// writes only its own, so it computes exactly what it would in
/// [`solve_in_turn`]: the results are bitwise the same, whatever the number
/// of threads.
#[cfg(feature = "rayon")]
fn solve_in_parallel<E: Equations + ?Sized>(
    eqs: &E,
    params: &mut [f64],
    fixed: &[bool],
    plans: Option<&Plans>,
    runs: &[&[&[usize]]],
    options: &CoreOptions,
) -> Vec<Solved> {
    use rayon::prelude::*;

    let parts: Vec<(Vec<f64>, Vec<Solved>)> = {
        let start: &[f64] = params;
        runs.par_iter()
            .map(|run| {
                let mut local = start.to_vec();
                let done = solve_in_turn(eqs, &mut local, fixed, plans, run, options);
                (local, done)
            })
            .collect()
    };
    let mut solved = Vec::new();
    for (local, done) in parts {
        for one in done {
            for &i in &one.moved {
                params[i] = local[i];
            }
            solved.push(one);
        }
    }
    solved
}

/// The components `groups` cut into runs of consecutive components, the
/// tasks of [`solve_in_parallel`]: one per [`PARALLEL_RUN_EQUATIONS`]
/// equations, at most one per thread of the pool, each holding about the
/// same number of equations. A single run means there is too little to
/// share out, and the components are solved in turn.
///
/// A task has to be worth waking a thread for, and the slowest task sets
/// the time. Measured on this crate's fixtures (24 threads): one task per
/// component cost 20–65 µs more than solving in turn on systems of 2–10
/// small components, broke even at 20 rectangles and gained only 30% at
/// 200. Tasks of at least 64, 96 and 128 equations, the work shared out
/// evenly, all made 200 rectangles about twice as fast and 40 and 80
/// grounded bodies 10–17% faster; 20 rectangles (360 equations) gained 21%,
/// 25% and 8% — at 128 they make two tasks — so 96 it is. Below 192
/// equations nothing is shared out.
#[cfg(feature = "rayon")]
fn parallel_runs<'g, E: Equations + ?Sized>(
    eqs: &E,
    groups: &'g [&'g [usize]],
) -> Vec<&'g [&'g [usize]]> {
    let sizes: Vec<usize> = groups
        .iter()
        .map(|group| equation_total(eqs, group))
        .collect();
    let total: usize = sizes.iter().sum();
    let tasks = (total / PARALLEL_RUN_EQUATIONS).min(rayon::current_num_threads());
    if tasks < 2 {
        return vec![groups];
    }
    // Cut where the running total passes each multiple of an equal share.
    let mut runs = Vec::with_capacity(tasks);
    let (mut begin, mut held) = (0, 0);
    for (k, size) in sizes.iter().enumerate() {
        held += size;
        if runs.len() + 1 < tasks && held * tasks >= (runs.len() + 1) * total {
            runs.push(&groups[begin..=k]);
            begin = k + 1;
        }
    }
    if begin < groups.len() {
        runs.push(&groups[begin..]);
    }
    runs
}

/// The equations per task of the parallel solve (see [`parallel_runs`]).
#[cfg(feature = "rayon")]
const PARALLEL_RUN_EQUATIONS: usize = 96;

/// What the solve of one component needs that depends on the structure of
/// the system alone — which constraints are active, which parameters they
/// list, which are held — and not on any value: kept by [`Plans`] from one
/// solve to the next while the topology holds.
pub(crate) struct Plan {
    /// The columns the constraints touch, held or free, ascending.
    touched: Vec<usize>,
    /// The layout over every free parameter: made with the plan for a
    /// component without classes, and only if a solve needs it (it falls
    /// back to the whole system) for one with them.
    whole: OnceLock<(Assembly, SparseJacobian)>,
    /// With classes of aliased parameters (see [`crate::alias`]): which
    /// rows they absorb, the layout over the classes and its step space.
    classes: Option<ClassPlan>,
    /// The blocks the component falls into (see [`crate::structure`]),
    /// when there are at least two and the family steps additively.
    blocks: Option<Vec<BlockPlan>>,
}

/// The part of a [`Plan`] for a component solved over classes.
struct ClassPlan {
    /// Per equation, whether its classes absorb it (the same at every
    /// solve: it depends on the structure only, not on the offsets).
    absorbed: Vec<bool>,
    assembly: Assembly,
    pattern: SparseJacobian,
    space: ClassSpace,
}

impl Plan {
    /// The plan of the constraints `group`, and the classes found at
    /// `params` while making it (for the solve at hand to use).
    fn new<E: Equations + ?Sized>(
        eqs: &E,
        group: &[usize],
        params: &[f64],
        fixed: &[bool],
        scratch: &mut Vec<SparseRow>,
    ) -> (Self, Option<Classes>) {
        let listing = Listing::new(eqs, group, params, scratch);
        let touched = listing.touched.clone();
        let classes = Aliases::classes_of(eqs, group, &touched, params, fixed);
        let whole = OnceLock::new();
        let plan = match &classes {
            Some(classes) => {
                let (assembly, pattern) = Assembly::lay_out(
                    listing,
                    fixed,
                    Shape {
                        classes: Some(classes),
                        empty: Some(&classes.absorbed),
                        ..Shape::default()
                    },
                );
                let space = ClassSpace::new(classes, &assembly.free, fixed);
                Some(ClassPlan {
                    absorbed: classes.absorbed.clone(),
                    assembly,
                    pattern,
                    space,
                })
            }
            None => {
                let _ = whole.set(Assembly::lay_out(listing, fixed, Shape::default()));
                None
            }
        };
        let blocks = if eqs.steps_are_additive() {
            match (&plan, &classes) {
                (Some(class_plan), Some(classes)) => block_plans(
                    eqs,
                    group,
                    params,
                    fixed,
                    Some(classes),
                    &class_plan.pattern,
                    &class_plan.assembly.free,
                    &classes.absorbed,
                    scratch,
                ),
                _ => whole.get().and_then(|(assembly, pattern)| {
                    block_plans(
                        eqs,
                        group,
                        params,
                        fixed,
                        None,
                        pattern,
                        &assembly.free,
                        &vec![false; pattern.n_rows()],
                        scratch,
                    )
                }),
            }
        } else {
            None
        };
        (
            Self {
                touched,
                whole,
                classes: plan,
                blocks,
            },
            classes,
        )
    }

    /// The layout over every free parameter of the constraints `group`
    /// (listed afresh if the plan does not have it: what an equation lists
    /// does not depend on `params`).
    fn whole<E: Equations + ?Sized>(
        &self,
        eqs: &E,
        group: &[usize],
        params: &[f64],
        fixed: &[bool],
        scratch: &mut Vec<SparseRow>,
    ) -> &(Assembly, SparseJacobian) {
        self.whole.get_or_init(|| {
            Assembly::lay_out(
                Listing::new(eqs, group, params, scratch),
                fixed,
                Shape::default(),
            )
        })
    }

    /// [`Plan::whole`], taken out of a plan that is not kept.
    fn into_whole<E: Equations + ?Sized>(
        self,
        eqs: &E,
        group: &[usize],
        params: &[f64],
        fixed: &[bool],
        scratch: &mut Vec<SparseRow>,
    ) -> (Assembly, SparseJacobian) {
        match self.whole.into_inner() {
            Some(whole) => whole,
            None => Assembly::lay_out(
                Listing::new(eqs, group, params, scratch),
                fixed,
                Shape::default(),
            ),
        }
    }
}

/// The connected components of a system and the [`Plan`] of each, made at
/// a component's first solve. Valid while the structure of the system does
/// not change; the owner keys it by its topology revision.
pub(crate) struct Plans {
    components: Vec<Vec<usize>>,
    /// Per constraint, the component it belongs to (`usize::MAX` for none).
    component_of: Vec<usize>,
    plans: Vec<OnceLock<Plan>>,
}

impl Plans {
    /// The components of `eqs` at `params`, with no plan made yet.
    pub(crate) fn new<E: Equations + ?Sized>(eqs: &E, params: &[f64]) -> Self {
        let components = components(eqs, params);
        let mut component_of = vec![usize::MAX; eqs.constraint_count()];
        for (ordinal, group) in components.iter().enumerate() {
            for &i in group {
                component_of[i] = ordinal;
            }
        }
        let plans = components.iter().map(|_| OnceLock::new()).collect();
        Self {
            components,
            component_of,
            plans,
        }
    }

    /// The components, as [`components`] gives them.
    pub(crate) fn components(&self) -> &[Vec<usize>] {
        &self.components
    }

    /// The slot of the plan of the component `group` belongs to.
    fn slot(&self, group: &[usize]) -> Option<&OnceLock<Plan>> {
        let ordinal = *self.component_of.get(*group.first()?)?;
        self.plans.get(ordinal)
    }
}

/// One block of a component (see [`crate::structure`]), laid out for its
/// own solve: its equations over its own columns, the rest held.
struct BlockPlan {
    /// The constraints its equations belong to, in the order of the group.
    sub: Vec<usize>,
    /// Per equation of `sub`, whether it is one of the block's.
    keep: Vec<bool>,
    /// How many are.
    kept: usize,
    assembly: Assembly,
    pattern: SparseJacobian,
    space: BlockSpace,
}

/// The step space of a block: over classes, or over plain parameters.
enum BlockSpace {
    Classes(ClassSpace),
    Plain(AdditiveSpace),
}

impl BlockSpace {
    fn space(&self) -> &dyn StepSpace {
        match self {
            BlockSpace::Classes(space) => space,
            BlockSpace::Plain(space) => space,
        }
    }
}

/// The blocks of the component `group`, laid out on `pattern` (over the
/// columns `keys`; rows `skip` marks are left out), or `None` when it does
/// not fall into at least two.
#[allow(clippy::too_many_arguments)]
fn block_plans<E: Equations + ?Sized>(
    eqs: &E,
    group: &[usize],
    params: &[f64],
    fixed: &[bool],
    classes: Option<&Classes>,
    pattern: &SparseJacobian,
    keys: &[usize],
    skip: &[bool],
    scratch: &mut Vec<SparseRow>,
) -> Option<Vec<BlockPlan>> {
    let blocks = structure::blocks(pattern, keys.len(), skip);
    if blocks.len() < 2 {
        return None;
    }
    // The first row of each constraint of the group, and the constraint
    // (by position in the group) of each row.
    let mut first = Vec::with_capacity(group.len() + 1);
    let mut owner = Vec::new();
    first.push(0);
    for (at, &i) in group.iter().enumerate() {
        let count = eqs.equation_count(i);
        owner.extend(std::iter::repeat_n(at, count));
        first.push(first[at] + count);
    }
    let members = classes.map(|classes| ClassSpace::members(classes, fixed));
    let plans = blocks
        .iter()
        .map(|block| {
            let block_keys: Vec<usize> = block.cols.iter().map(|&k| keys[k]).collect();
            let mut positions: Vec<usize> = block.rows.iter().map(|&r| owner[r]).collect();
            positions.dedup();
            let sub: Vec<usize> = positions.iter().map(|&at| group[at]).collect();
            let mut keep = Vec::with_capacity(equation_total(eqs, &sub));
            for &at in &positions {
                for r in first[at]..first[at + 1] {
                    keep.push(block.rows.binary_search(&r).is_ok());
                }
            }
            let listing = Listing::new(eqs, &sub, params, scratch);
            let shape = Shape {
                classes,
                empty: None,
                keep: Some(&keep),
                columns: Some(&block_keys),
            };
            let (assembly, pattern) = Assembly::lay_out(listing, fixed, shape);
            let space = match &members {
                Some(members) => BlockSpace::Classes(ClassSpace::of_roots(members, &assembly.free)),
                None => BlockSpace::Plain(AdditiveSpace::new(assembly.free.clone())),
            };
            BlockPlan {
                sub,
                kept: block.rows.len(),
                keep,
                assembly,
                pattern,
                space,
            }
        })
        .collect();
    Some(plans)
}

/// Solve the blocks `blocks` of a component one after the other, each for
/// its own columns with the others held, each from where it stands. Stops
/// at the first that does not converge and reports it; the caller puts
/// the parameters back.
fn solve_blocks<E: Equations + ?Sized>(
    eqs: &E,
    group: &[usize],
    blocks: &[BlockPlan],
    params: &mut [f64],
    unit: f64,
    scratch: &mut Vec<SparseRow>,
    options: &CoreOptions,
) -> SolveResult {
    let mut iterations = 0;
    let mut status = SolveStatus::Converged;
    for block in blocks {
        let result = match solve_small_block(eqs, block, params, unit, scratch, options) {
            Some(result) => result,
            None => {
                let space = block.space.space();
                let mut columns = |p: &[f64], out: &mut [f64]| match &block.space {
                    BlockSpace::Classes(space) => {
                        out.fill(0.0);
                        for (&member, &k) in space.storage().iter().zip(space.columns()) {
                            out[k] = out[k].max(eqs.column_magnitude(member, p));
                        }
                    }
                    BlockSpace::Plain(space) => {
                        for (magnitude, &column) in out.iter_mut().zip(space.storage()) {
                            *magnitude = eqs.column_magnitude(column, p);
                        }
                    }
                };
                newton_over(
                    eqs,
                    &block.sub,
                    Some(&block.keep),
                    params,
                    space,
                    &block.assembly,
                    block.pattern.clone(),
                    &mut columns,
                    unit,
                    scratch,
                    options,
                )
            }
        };
        iterations += result.iterations;
        status = status.max(result.status);
        if !result.converged {
            return SolveResult {
                iterations,
                ..result
            };
        }
    }
    SolveResult {
        converged: true,
        residual: numerics::rms(&eval_residuals(eqs, group, params)),
        iterations,
        status,
    }
}

/// The largest square block solved by [`solve_small_block`].
const SMALL_BLOCK: usize = 3;

/// A square block of at most [`SMALL_BLOCK`] unknowns solved by plain
/// Newton steps on its dense Jacobian, each halved until it lowers the
/// residual, from where the block stands — so it settles on the root
/// nearest to that, as the block's equations see it. Converged, it takes
/// further steps while they strictly lower the residual (up to the polish
/// count), as the general loop does. `None` when it does not apply or does
/// not get there (a singular step, no step that helps, the iteration limit,
/// an equation that cannot be met to the tolerance where its numbers
/// stand): the general loop then takes the block from where it started.
fn solve_small_block<E: Equations + ?Sized>(
    eqs: &E,
    block: &BlockPlan,
    params: &mut [f64],
    unit: f64,
    scratch: &mut Vec<SparseRow>,
    options: &CoreOptions,
) -> Option<SolveResult> {
    let space = block.space.space();
    let n = space.dim();
    if n == 0 || n > SMALL_BLOCK || block.kept != n || options.norm != numerics::ResidualNorm::Rms {
        return None;
    }
    let tolerance = options.tolerance * unit;
    let storage = space.storage();
    let start: Vec<f64> = storage.iter().map(|&p| params[p]).collect();
    let mut all = vec![0.0_f64; block.keep.len()];
    let mut eval = |p: &[f64], f: &mut [f64; SMALL_BLOCK]| {
        eval_residuals_into(eqs, &block.sub, p, &mut all);
        let mut at = 0;
        for (v, &k) in all.iter().zip(&block.keep) {
            if k {
                f[at] = *v;
                at += 1;
            }
        }
        (f[..n].iter().map(|v| v * v).sum::<f64>() / n as f64).sqrt()
    };
    let give_up = |params: &mut [f64]| {
        for (&p, &v) in storage.iter().zip(&start) {
            params[p] = v;
        }
        None
    };
    let mut j = block.pattern.clone();
    let mut f = [0.0_f64; SMALL_BLOCK];
    let mut residual = eval(params, &mut f);
    let mut iterations = 0;
    let mut polish = 0;
    loop {
        if !residual.is_finite() {
            return give_up(params);
        }
        let converged = residual < tolerance;
        if converged && (residual == 0.0 || polish >= options.polish_steps) {
            break;
        }
        if !converged && iterations >= options.max_iterations {
            return give_up(params);
        }
        // The dense Jacobian of the block.
        for v in j.values_mut() {
            *v = 0.0;
        }
        block
            .assembly
            .assemble(eqs, &block.sub, params, &mut j, scratch);
        let mut a = [[0.0_f64; SMALL_BLOCK]; SMALL_BLOCK];
        for (r, row) in a.iter_mut().enumerate().take(n) {
            for (&c, &v) in j.row_cols(r).iter().zip(j.row_values(r)) {
                row[c] = v;
            }
        }
        let Some(step) = solve_dense(&mut a, &f, n) else {
            if converged {
                break;
            }
            return give_up(params);
        };
        let saved: Vec<f64> = storage.iter().map(|&p| params[p]).collect();
        let mut taken = None;
        let mut scale = 1.0;
        for _ in 0..if converged { 1 } else { 20 } {
            let trial: Vec<f64> = step[..n].iter().map(|h| h * scale).collect();
            space.apply(params, &saved, &trial);
            let mut f_new = [0.0_f64; SMALL_BLOCK];
            let r_new = eval(params, &mut f_new);
            if r_new < residual {
                taken = Some((f_new, r_new));
                break;
            }
            scale *= 0.5;
        }
        match taken {
            Some((f_new, r_new)) => {
                f = f_new;
                residual = r_new;
                if converged {
                    polish += 1;
                } else {
                    iterations += 1;
                }
            }
            None => {
                for (&p, &v) in storage.iter().zip(&saved) {
                    params[p] = v;
                }
                if converged {
                    break;
                }
                return give_up(params);
            }
        }
    }
    Some(SolveResult {
        converged: true,
        residual,
        iterations,
        status: SolveStatus::Converged,
    })
}

/// Solve `a·h = −f` for the leading `n × n` of `a` by Gaussian elimination
/// with partial pivoting; `None` if a pivot is not above `PIVOT_TOL`
/// relative to the largest entry.
fn solve_dense(
    a: &mut [[f64; SMALL_BLOCK]; SMALL_BLOCK],
    f: &[f64; SMALL_BLOCK],
    n: usize,
) -> Option<[f64; SMALL_BLOCK]> {
    let mut b = [0.0_f64; SMALL_BLOCK];
    for i in 0..n {
        b[i] = -f[i];
    }
    let scale = a[..n]
        .iter()
        .flat_map(|row| row[..n].iter())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    if scale == 0.0 || !scale.is_finite() {
        return None;
    }
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&x, &y| a[x][k].abs().total_cmp(&a[y][k].abs()))
            .expect("rows left");
        if a[pivot][k].abs() <= numerics::PIVOT_TOL * scale {
            return None;
        }
        a.swap(k, pivot);
        b.swap(k, pivot);
        for i in k + 1..n {
            let m = a[i][k] / a[k][k];
            for c in k..n {
                a[i][c] -= m * a[k][c];
            }
            b[i] -= m * b[k];
        }
    }
    let mut h = [0.0_f64; SMALL_BLOCK];
    for k in (0..n).rev() {
        let mut v = b[k];
        for c in k + 1..n {
            v -= a[k][c] * h[c];
        }
        h[k] = v / a[k][k];
    }
    h.iter().all(|v| v.is_finite()).then_some(h)
}

/// Solve one connected subsystem: only the free parameters its constraints
/// touch may move.
///
/// The equations that only alias one parameter to another or to a value
/// are eliminated first (see [`crate::alias`]): the parameters are put where
/// those equations say, and the solve runs over one column per class of
/// parameters still unknown. If that solve does not converge, the
/// parameters are put back and the component is solved over every
/// parameter, as without the elimination.
///
/// What depends on the structure alone comes from the component's [`Plan`]
/// in `plans`, made now if it is not there yet; without `plans` it is made
/// for this solve and dropped.
fn solve_one<E: Equations + ?Sized>(
    eqs: &E,
    group: &[usize],
    plans: Option<&Plans>,
    params: &mut [f64],
    fixed: &[bool],
    options: &CoreOptions,
) -> Solved {
    let mut scratch: Vec<SparseRow> = Vec::new();
    let unit = equation_unit(eqs.length_scale(group, params));
    // A plan that is kept, or one made for this solve alone (which the
    // solve may then take apart instead of copying from).
    let mut local: Option<Plan> = None;
    let mut fresh: Option<Classes> = None;
    let kept: Option<&Plan> = match plans.and_then(|plans| plans.slot(group)) {
        Some(slot) => Some(match slot.get() {
            Some(plan) => plan,
            None => {
                let (made, classes) = Plan::new(eqs, group, params, fixed, &mut scratch);
                fresh = classes;
                // Components are solved once each per call, so nobody else
                // fills this slot meanwhile; if somebody did, theirs is the
                // same plan.
                let _ = slot.set(made);
                slot.get().expect("the plan was just set")
            }
        }),
        None => {
            let (made, classes) = Plan::new(eqs, group, params, fixed, &mut scratch);
            fresh = classes;
            local = Some(made);
            None
        }
    };
    let plan: &Plan = kept.unwrap_or_else(|| local.as_ref().expect("made above"));
    let touched_owned = plan.touched.clone();
    let touched = &touched_owned;

    // One-shot integration check (debug builds): the analytic assembly
    // must agree with central differences at the starting point, over
    // every parameter and before anything moves. The exhaustive
    // per-variant gates live with the constraint definitions; this catches
    // per-call assembly bugs (column mapping, mask filtering).
    #[cfg(debug_assertions)]
    if equation_total(eqs, group) > 0 && eqs.debug_check_jacobian() {
        let (assembly, pattern) = plan.whole(eqs, group, params, fixed, &mut scratch);
        let free = &assembly.free;
        if !free.is_empty() {
            let space = eqs.step_space(free);
            let mut build = |p: &[f64], j: &mut SparseJacobian| {
                assembly.assemble(eqs, group, p, j, &mut scratch);
            };
            debug_check_jacobian(eqs, group, params, free, &*space, pattern, &mut build);
        }
    }

    // A start that is not a number is reported as it stands (`Diverged`),
    // not mended by putting the parameters where their classes say.
    let finite = touched.iter().all(|&p| params[p].is_finite());
    if let (true, Some(class_plan)) = (finite, &plan.classes) {
        // The offsets and the known values follow the dimensions, so the
        // classes are found afresh; what they absorb does not change.
        let classes = fresh.or_else(|| Aliases::classes_of(eqs, group, touched, params, fixed));
        if let Some(classes) = classes {
            debug_assert_eq!(classes.absorbed, class_plan.absorbed);
            let movable: Vec<usize> = touched.iter().copied().filter(|&p| !fixed[p]).collect();
            let start: Vec<f64> = movable.iter().map(|&p| params[p]).collect();
            classes.snap(params, fixed);
            if let Some(blocks) = &plan.blocks {
                let result = solve_blocks(eqs, group, blocks, params, unit, &mut scratch, options);
                if result.converged {
                    return Solved {
                        result,
                        #[cfg(feature = "rayon")]
                        moved: movable,
                        touched: touched_owned,
                    };
                }
                for (&p, &v) in movable.iter().zip(&start) {
                    params[p] = v;
                }
                classes.snap(params, fixed);
            }
            let space = &class_plan.space;
            let mut columns = |p: &[f64], out: &mut [f64]| {
                out.fill(0.0);
                for (&member, &k) in space.storage().iter().zip(space.columns()) {
                    out[k] = out[k].max(eqs.column_magnitude(member, p));
                }
            };
            let result = newton_over(
                eqs,
                group,
                None,
                params,
                space,
                &class_plan.assembly,
                class_plan.pattern.clone(),
                &mut columns,
                unit,
                &mut scratch,
                options,
            );
            if result.converged {
                return Solved {
                    result,
                    #[cfg(feature = "rayon")]
                    moved: movable,
                    touched: touched_owned,
                };
            }
            for (&p, &v) in movable.iter().zip(&start) {
                params[p] = v;
            }
        }
    }

    // Only this component's own free columns move: the Newton loop's free
    // list is exactly its movable coordinates.
    if let (true, None, Some(blocks)) = (finite, &plan.classes, &plan.blocks) {
        let movable: Vec<usize> = touched.iter().copied().filter(|&p| !fixed[p]).collect();
        let start: Vec<f64> = movable.iter().map(|&p| params[p]).collect();
        let result = solve_blocks(eqs, group, blocks, params, unit, &mut scratch, options);
        if result.converged {
            return Solved {
                result,
                #[cfg(feature = "rayon")]
                moved: movable,
                touched: touched_owned,
            };
        }
        for (&p, &v) in movable.iter().zip(&start) {
            params[p] = v;
        }
    }
    let owned: Assembly;
    let (assembly, pattern) = match (kept, local.take()) {
        (Some(plan), _) => {
            let (assembly, pattern) = plan.whole(eqs, group, params, fixed, &mut scratch);
            (assembly, pattern.clone())
        }
        (None, Some(plan)) => {
            let (assembly, pattern) = plan.into_whole(eqs, group, params, fixed, &mut scratch);
            owned = assembly;
            (&owned, pattern)
        }
        (None, None) => unreachable!("a plan is kept or made"),
    };
    let free = &assembly.free;
    let space = eqs.step_space(free);
    let mut columns = |p: &[f64], out: &mut [f64]| {
        for (magnitude, &column) in out.iter_mut().zip(free) {
            *magnitude = eqs.column_magnitude(column, p);
        }
    };
    let result = newton_over(
        eqs,
        group,
        None,
        params,
        &*space,
        assembly,
        pattern,
        &mut columns,
        unit,
        &mut scratch,
        options,
    );
    Solved {
        result,
        #[cfg(feature = "rayon")]
        moved: space.storage().to_vec(),
        touched: touched_owned,
    }
}

/// The Newton solve of the constraints `group` over the columns of
/// `space`, whose Jacobian `assembly` builds into `pattern`; with `keep`,
/// of only the equations it marks (one flag per equation of `group`).
#[allow(clippy::too_many_arguments)]
fn newton_over<E: Equations + ?Sized>(
    eqs: &E,
    group: &[usize],
    keep: Option<&[bool]>,
    params: &mut [f64],
    space: &dyn StepSpace,
    assembly: &Assembly,
    pattern: SparseJacobian,
    columns: &mut dyn FnMut(&[f64], &mut [f64]),
    unit: f64,
    scratch: &mut Vec<SparseRow>,
    options: &CoreOptions,
) -> SolveResult {
    let mut build = |p: &[f64], j: &mut SparseJacobian| {
        assembly.assemble(eqs, group, p, j, scratch);
    };
    match keep {
        None => {
            let n_eq = equation_total(eqs, group);
            let mut rows = |p: &[f64], out: &mut [f64]| row_magnitudes(eqs, group, p, out);
            numerics::newton_solve_in(
                params,
                space,
                n_eq,
                |p, f| eval_residuals_into(eqs, group, p, f),
                numerics::Jacobian::Analytic {
                    pattern,
                    build: &mut build,
                },
                Some(numerics::Magnitudes {
                    columns,
                    rows: &mut rows,
                    unit,
                }),
                options,
            )
        }
        Some(keep) => {
            // Every equation of the constraints is evaluated; the kept ones
            // are passed on.
            let n_eq = keep.iter().filter(|&&k| k).count();
            let mut all = vec![0.0_f64; keep.len()];
            let mut all_rows = vec![0.0_f64; keep.len()];
            let compact = |all: &[f64], out: &mut [f64]| {
                let mut at = 0;
                for (v, &k) in all.iter().zip(keep) {
                    if k {
                        out[at] = *v;
                        at += 1;
                    }
                }
            };
            let mut rows = |p: &[f64], out: &mut [f64]| {
                row_magnitudes(eqs, group, p, &mut all_rows);
                compact(&all_rows, out);
            };
            numerics::newton_solve_in(
                params,
                space,
                n_eq,
                |p, f| {
                    eval_residuals_into(eqs, group, p, &mut all);
                    compact(&all, f);
                },
                numerics::Jacobian::Analytic {
                    pattern,
                    build: &mut build,
                },
                Some(numerics::Magnitudes {
                    columns,
                    rows: &mut rows,
                    unit,
                }),
                options,
            )
        }
    }
}

/// The unit a component's equations are measured in, from its
/// [`Equations::length_scale`]: that size when it is below one, and one
/// otherwise — a model larger than a unit keeps the absolute tolerance,
/// which the precision of its numbers already bounds (see
/// [`numerics::PRECISION_ULPS`]). A size of zero (every point in one
/// place) or none at all gives one too.
pub(crate) fn equation_unit(length_scale: Option<f64>) -> f64 {
    match length_scale {
        Some(size) if size > 0.0 && size < 1.0 => size,
        _ => 1.0,
    }
}

/// The debug-build cross-check of [`solve_one`]: the analytic Jacobian of a
/// component against central differences of its residuals.
///
/// A difference is only evidence where it can be trusted, and it cannot
/// always: a step of 1e-7 is a tenth of a feature one micron across
/// (truncation), and it is a single float spacing a billion units from the
/// origin (rounding). So a mismatch is put to a second difference at half
/// the step. The two extrapolate to a better estimate (the truncation error
/// falls by four when the step halves); an analytic entry that agrees with
/// that is right. Otherwise, if the two differences agree closely with each
/// other, they are resolving the derivative and the analytic entry is what
/// is wrong: this panics. If they do not, the entry is left unjudged.
#[cfg(debug_assertions)]
fn debug_check_jacobian<E: Equations + ?Sized>(
    eqs: &E,
    group: &[usize],
    params: &mut [f64],
    free: &[usize],
    space: &dyn StepSpace,
    pattern: &SparseJacobian,
    build: &mut numerics::JacobianBuilder<'_>,
) {
    let n_eq = equation_total(eqs, group);
    let mut eval = |p: &[f64], f: &mut [f64]| eval_residuals_into(eqs, group, p, f);
    let mut built = pattern.clone();
    build(params, &mut built);
    let j_analytic = built.to_dense();
    let mut j_fd = vec![vec![0.0_f64; free.len()]; n_eq];
    numerics::numerical_jacobian_for(params, space, n_eq, &mut eval, &mut j_fd);

    // Rounding of the residuals, blown up by the division by the step.
    let mut size = vec![0.0_f64; n_eq];
    row_magnitudes(eqs, group, params, &mut size);
    let allowed = |row: usize, fd: f64| {
        1.0e-4 * (1.0 + fd.abs()) + 4.0 * f64::EPSILON * size[row] / (2.0 * numerics::FD_H)
    };
    let additive = space.additive_params().is_some();

    for ci in 0..free.len() {
        let off: Vec<usize> = (0..n_eq)
            .filter(|&ri| (j_analytic[ri][ci] - j_fd[ri][ci]).abs() > allowed(ri, j_fd[ri][ci]))
            .collect();
        if off.is_empty() {
            continue;
        }
        // A second opinion at half the step, for plain parameters. (The
        // rigid-pose family, whose columns are not parameters, does not run
        // this check.)
        let second: Option<Vec<f64>> = additive.then(|| {
            let gi = free[ci];
            let orig = params[gi];
            let h = 0.5 * numerics::fd_step(orig);
            params[gi] = orig + h;
            let above = params[gi];
            let f_plus = eval_residuals(eqs, group, params);
            params[gi] = orig - h;
            let below = params[gi];
            let f_minus = eval_residuals(eqs, group, params);
            params[gi] = orig;
            f_plus
                .iter()
                .zip(&f_minus)
                .map(|(a, b)| (a - b) / (above - below))
                .collect()
        });
        for ri in off {
            let (va, vf) = (j_analytic[ri][ci], j_fd[ri][ci]);
            let wrong = match &second {
                None => true,
                Some(half) => {
                    let extrapolated = (4.0 * half[ri] - vf) / 3.0;
                    let limit = allowed(ri, vf);
                    (va - extrapolated).abs() > limit && (half[ri] - vf).abs() <= 0.25 * limit
                }
            };
            assert!(
                !wrong,
                "analytic Jacobian drifts from finite differences at row {ri} col {ci}: {va} vs {vf}"
            );
        }
    }
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

        // What the thresholds are relative to: the larger of the largest
        // residual and the unit the equations are measured in (one, or the
        // size of a component smaller than one; see `equation_unit`). With
        // a floor of one, a contradiction of 0.2% on a triangle 1e-6 across
        // read as partially redundant.
        let unit = equation_unit(eqs.length_scale(group, params));

        // Nothing in this component can move: any unsatisfied residual is a
        // genuine conflict.
        let flag_unsatisfied = |d: &mut Diagnosis, f_c: &[f64]| {
            let scale = f_c.iter().fold(unit, |m, v| m.max(v.abs()));
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

        let scale = f_c.iter().fold(unit, |m, v| m.max(v.abs()));
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

#[cfg(all(test, feature = "rayon"))]
pub(crate) mod parallel_tests {
    use super::*;
    use crate::constraint::Constraint;
    use crate::entity::PointRef;
    use crate::solver::SketchEquations;
    use crate::Sketch;

    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|x| x.to_bits()).collect()
    }

    /// Solve `groups` in turn and in parallel, one run per component (the
    /// most threads the parallel path can use), from the same start; both
    /// must leave the same parameters, bit for bit, and report the same.
    pub(crate) fn assert_paths_agree<E: Equations + ?Sized>(
        eqs: &E,
        start: &[f64],
        fixed: &[bool],
        groups: &[&[usize]],
    ) {
        let options = CoreOptions::from(&crate::SolveOptions::default());
        let mut in_turn = start.to_vec();
        let a = solve_in_turn(eqs, &mut in_turn, fixed, None, groups, &options);
        let runs: Vec<&[&[usize]]> = groups.chunks(1).collect();
        let mut parallel = start.to_vec();
        let b = solve_in_parallel(eqs, &mut parallel, fixed, None, &runs, &options);
        assert_ne!(bits(&in_turn), bits(start), "the start must need solving");
        assert_eq!(bits(&in_turn), bits(&parallel));
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().zip(&b) {
            assert_eq!(a.result.status, b.result.status);
            assert_eq!(a.result.iterations, b.result.iterations);
            assert_eq!(a.result.residual.to_bits(), b.result.residual.to_bits());
            assert_eq!(a.touched, b.touched);
        }
    }

    /// `count` independent anchored segments of six equations each, knocked
    /// off their solution.
    fn segments(count: usize) -> Sketch {
        let mut sk = Sketch::new();
        for i in 0..count {
            let (x, y) = ((i % 4) as f64 * 10.0, (i / 4) as f64 * 10.0);
            let origin = sk.add_point(x, y);
            let side = sk.add_line(x, y, x + 3.0, y);
            sk.add_constraint(Constraint::fix_point(origin, x, y));
            sk.add_constraint(Constraint::coincident(
                PointRef::LineStart(side),
                PointRef::Point(origin),
            ));
            sk.add_constraint(Constraint::Horizontal(side));
            sk.add_constraint(Constraint::LineLength {
                line: side,
                length: 3.0 + i as f64,
            });
        }
        sk.drift_free_params(0.1, 3);
        sk
    }

    #[test]
    fn parallel_components_solve_bitwise_as_in_turn_2d() {
        let sk = segments(12);
        let eqs = SketchEquations::new(&sk.entities, &sk.constraints);
        let components = components(&eqs, &sk.params);
        assert_eq!(components.len(), 12);
        let groups: Vec<&[usize]> = components.iter().map(Vec::as_slice).collect();
        assert_paths_agree(&eqs, &sk.params, &sk.fixed, &groups);
    }

    /// One task per `PARALLEL_RUN_EQUATIONS` equations and at most one per
    /// thread, the components kept in order and shared out evenly; a
    /// system too small for two tasks stays in one run.
    #[test]
    fn parallel_runs_share_the_equations_out_evenly() {
        let threads = rayon::current_num_threads();
        for count in [1, 15, 16, 31, 32, 33, 50, 97, 400] {
            let sk = segments(count);
            let eqs = SketchEquations::new(&sk.entities, &sk.constraints);
            let components = components(&eqs, &sk.params);
            let groups: Vec<&[usize]> = components.iter().map(Vec::as_slice).collect();
            let runs = parallel_runs(&eqs, &groups);

            let total = 6 * count;
            let tasks = (total / PARALLEL_RUN_EQUATIONS).min(threads).max(1);
            assert_eq!(runs.len(), tasks, "{count} components");
            let rejoined: Vec<&[usize]> = runs.iter().flat_map(|run| run.iter().copied()).collect();
            assert_eq!(rejoined, groups, "{count} components");
            for run in &runs {
                let held: usize = run.iter().map(|group| equation_total(&eqs, group)).sum();
                assert!(
                    held.abs_diff(total / tasks) <= 6,
                    "{count} components: a run of {held} equations for a share of {}",
                    total / tasks
                );
            }
        }
    }
}

#[cfg(test)]
mod small_block_tests {
    use super::{solve_dense, SMALL_BLOCK};

    #[test]
    fn a_small_system_is_solved_with_pivoting() {
        // a·h = −f with a zero leading entry (needs a row exchange).
        let mut a = [[0.0, 2.0, 1.0], [1.0, 1.0, 0.0], [3.0, 0.0, 1.0]];
        let h = [1.0, -2.0, 0.5];
        let mut f = [0.0_f64; SMALL_BLOCK];
        for (fi, row) in f.iter_mut().zip(&a) {
            *fi = -(row[0] * h[0] + row[1] * h[1] + row[2] * h[2]);
        }
        let got = solve_dense(&mut a, &f, 3).expect("regular");
        for (g, e) in got.iter().zip(h) {
            assert!((g - e).abs() < 1e-14, "{got:?}");
        }
    }

    #[test]
    fn a_singular_or_empty_system_has_no_step() {
        let mut a = [[1.0, 2.0, 0.0], [2.0, 4.0, 0.0], [0.0, 0.0, 0.0]];
        assert!(solve_dense(&mut a, &[1.0, 1.0, 0.0], 2).is_none());
        let mut zero = [[0.0; SMALL_BLOCK]; SMALL_BLOCK];
        assert!(solve_dense(&mut zero, &[1.0, 0.0, 0.0], 1).is_none());
    }
}
