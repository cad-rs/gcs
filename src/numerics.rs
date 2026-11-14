//! The shared Newton-Raphson least-squares core.
//!
//! All three families run the same loop: evaluate residuals, build the
//! Jacobian (analytic rows supplied by the family; central finite
//! differences for callers of [`newton_solve`] that have none), solve the
//! regularised normal equations `(JᵀJ + λI) Δx = Jᵀ(−F)`, and take a step
//! inside a trust radius. The residual and Jacobian evaluation differ per
//! family, so they are injected as closures; the rigid-pose family also
//! supplies how a step is applied to a pose.
//!
//! The Moore-Penrose pseudo-inverse J⁺ is realised through the normal
//! equations with Tikhonov regularisation λ: under-constrained systems get
//! the minimum-norm step instead of diverging.
//!
//! Step control is Powell's dogleg. The Gauss–Newton step is taken whole
//! while it fits inside the trust radius, which starts out unbounded, so a
//! well-behaved solve is plain Gauss–Newton. A step is judged by its gain:
//! the reduction of the squared residual it achieved over the reduction the
//! linear model promised. A step with no gain is not taken; the radius
//! shrinks and the next trial bends from the Gauss–Newton direction toward
//! steepest descent, which costs a residual evaluation and no new
//! factorisation. A solve that ends without converging reports why
//! ([`SolveStatus`]) and, by default, puts the starting parameters back
//! ([`FailurePolicy`]).

use std::cmp::Ordering;

/// Result returned by every [`solve`][`crate::solver::solve`] entry point.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub struct SolveResult {
    /// `true` if every solved component ended with its RMS constraint
    /// residual below the tolerance ([`SolveOptions::tolerance`],
    /// [`RESIDUAL_TOL`] by default), or with every equation as exact as the
    /// numbers it is computed from allow
    /// ([`SolveStatus::ConvergedToPrecision`]). A component with no free
    /// parameter is judged on the residual it has.
    pub converged: bool,
    /// RMS residual of the parameters the solve left behind. After a failed
    /// solve under [`FailurePolicy::Restore`] this is the residual of the
    /// restored starting state.
    pub residual: f64,
    /// Number of Newton iterations performed.
    pub iterations: usize,
    /// Why the solve ended: [`SolveStatus::Converged`] or
    /// [`SolveStatus::ConvergedToPrecision`] exactly when `converged` is
    /// `true`.
    pub status: SolveStatus,
}

/// Why a solve ended.
///
/// When several connected components are solved, the reported status is the
/// most severe one among them (the variants are declared in increasing
/// severity).
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SolveStatus {
    /// The residual is below the tolerance.
    Converged,
    /// The residual is above the tolerance, but every equation holds to
    /// within the rounding of the numbers it is computed from: at this
    /// model's coordinates the tolerance asks for more than double
    /// precision can represent. (A point a billion units from the origin is
    /// stored to about 1e-7; no equation on it can be met to 1e-10.) The
    /// solve counts as converged. To see this status less often, keep
    /// coordinates near the origin or pass a tolerance that fits the model.
    ConvergedToPrecision,
    /// The iteration limit was reached while the residual was still above
    /// the tolerance and still falling.
    MaxIterations,
    /// The iteration cannot move from here: no step, however short, lowers
    /// the residual (a stationary point of the residual norm that is not a
    /// solution), or the unsatisfied equations have no free parameter at
    /// all. This is how a conflicting system ends: at its least-squares
    /// compromise.
    Stalled,
    /// The regularised normal equations were numerically singular; no step
    /// could be computed.
    Singular,
    /// A residual or a step became NaN or infinite.
    Diverged,
}

/// What a failed solve leaves in the parameter vector.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FailurePolicy {
    /// Put back the parameters the solve started from, so a failed solve
    /// leaves the geometry untouched. Applies per connected component: a
    /// component that converged keeps its solution.
    #[default]
    Restore,
    /// Keep the last iterate (what 0.1.0 always did) — useful to inspect
    /// how far a failed solve got.
    KeepLast,
}

/// Options for the 2D-sketch and 3D-geometric solver families.
///
/// The rigid-pose family has its own [`crate::rigid3d::PoseSolveOptions`]:
/// its tolerance is an L2 norm and it has two more knobs.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, PartialEq)]
pub struct SolveOptions {
    /// Convergence tolerance on the RMS constraint residual.
    pub tolerance: f64,
    /// Maximum Newton iterations per connected component.
    pub max_iterations: usize,
    /// What a failed solve leaves in the parameter vector.
    pub on_failure: FailurePolicy,
}

impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            tolerance: RESIDUAL_TOL,
            max_iterations: MAX_ITER,
            on_failure: FailurePolicy::default(),
        }
    }
}

/// Which norm of the residual vector a tolerance is measured in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ResidualNorm {
    /// Root mean square (the 2D and 3D-geometric families).
    Rms,
    /// Euclidean norm (the rigid-pose family).
    L2,
}

/// What the Newton loop is told, whichever family is calling.
///
/// The public [`SolveOptions`] belongs to the 2D and 3D-geometric families
/// and the rigid-pose family has its own option type (an L2 tolerance, a
/// step scale, a base regularisation); both map onto this.
#[derive(Debug, Clone)]
pub(crate) struct CoreOptions {
    /// Convergence tolerance on the residual, measured in `norm`.
    pub(crate) tolerance: f64,
    /// The norm `tolerance` refers to.
    pub(crate) norm: ResidualNorm,
    /// Maximum Newton iterations.
    pub(crate) max_iterations: usize,
    /// What a failed solve leaves in the parameters.
    pub(crate) on_failure: FailurePolicy,
    /// Base Tikhonov regularisation of the normal equations.
    pub(crate) lambda: f64,
    /// The factor the Gauss–Newton step is scaled by (1 is the full step).
    pub(crate) initial_step: f64,
    /// Upper bound on convergence-polish steps (see [`POLISH_STEPS`]).
    pub(crate) polish_steps: usize,
    /// Whether the convergence polish falls back to a regularised step when
    /// the unregularised one does not exist (an under-constrained
    /// component). Off for the 2D and 3D-geometric families, which leave
    /// such a component at the first iterate under the tolerance.
    pub(crate) polish_singular: bool,
}

impl From<&SolveOptions> for CoreOptions {
    fn from(options: &SolveOptions) -> Self {
        Self {
            tolerance: options.tolerance,
            norm: ResidualNorm::Rms,
            max_iterations: options.max_iterations,
            on_failure: options.on_failure,
            lambda: LAMBDA,
            initial_step: 1.0,
            polish_steps: POLISH_STEPS,
            polish_singular: false,
        }
    }
}

/// `a > b`, and `false` when either side is NaN. Spelled through
/// `partial_cmp` so a guard written as `!gt(a, b)` visibly takes the NaN
/// case too.
#[inline]
pub(crate) fn gt(a: f64, b: f64) -> bool {
    a.partial_cmp(&b) == Some(Ordering::Greater)
}

/// `a <= b`, and `false` when either side is NaN (see [`gt`]).
#[inline]
pub(crate) fn le(a: f64, b: f64) -> bool {
    matches!(a.partial_cmp(&b), Some(Ordering::Less | Ordering::Equal))
}

/// Convergence tolerance on the RMS constraint residual.
pub const RESIDUAL_TOL: f64 = 1e-10;
/// Maximum Newton iterations.
pub const MAX_ITER: usize = 100;
/// Step size for numerical Jacobian (central differences).
pub const FD_H: f64 = 1e-7;
/// Tikhonov regularisation coefficient of the normal equations the step is
/// solved from. Where those are ill-conditioned or rank-deficient (their
/// smallest pivot below `1e-8` of the largest) the step is instead the
/// minimum-norm solution of a rank-truncated QR of the Jacobian, and this
/// coefficient does not enter it.
pub const LAMBDA: f64 = 1e-10;
/// Pivot threshold for Gaussian elimination.
pub const PIVOT_TOL: f64 = 1e-14;
/// Normalization guard threshold — prevents division by near-zero vectors.
pub const SOLVER_NORM_TOL: f64 = 1e-10;
/// How many rounding units of its own numbers an equation is allowed
/// before it counts as unmet, when the tolerance is finer than that: an
/// equation whose terms have magnitude `size = Σ|∂f/∂q|·|q|` is held to
/// `max(tolerance, PRECISION_ULPS · ε · size)`.
///
/// One rounding of the terms is worth `ε·size`. About two of them is what
/// can be reached at best: each parameter is stored to half a unit, and
/// evaluating the equation rounds a few more times. Four leaves a factor
/// of two over that. A larger allowance gives accuracy away where a tight
/// tolerance is still reachable (with 8, the assembly corpus, which asks
/// for 1e-12 at coordinates of some tens, came back five times less
/// exactly from its round trips; with 4, half again as loosely as with no
/// allowance at all); a smaller one risks never settling.
pub const PRECISION_ULPS: f64 = 4.0;
/// How often the equation weights are taken afresh when a solve ends
/// without converging (the geometry may have moved to larger numbers).
const REWEIGH_PASSES: usize = 2;
/// A step whose gain (achieved over promised reduction of the squared
/// residual) is above this widens the trust radius.
const TRUST_WIDEN_ABOVE: f64 = 0.75;
/// A step whose gain is at or below this narrows the trust radius.
const TRUST_NARROW_BELOW: f64 = 0.25;
/// How far a good step lets the radius grow, as a multiple of its length.
const TRUST_GROWTH: f64 = 3.0;
/// Trial steps per iteration before the iterate is declared stationary.
/// Every rejection halves the radius, so the last trial is some twelve
/// orders of magnitude shorter than the first.
const TRUST_ATTEMPTS: usize = 40;
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
///
/// Uses the default [`SolveOptions`]; a solve that does not converge puts
/// the starting parameters back (see [`FailurePolicy`]).
pub fn newton_solve<F>(
    params: &mut [f64],
    fixed: &[bool],
    n_eq: usize,
    eval_residuals: F,
) -> SolveResult
where
    F: FnMut(&[f64], &mut [f64]),
{
    newton_solve_with(
        params,
        fixed,
        n_eq,
        eval_residuals,
        &SolveOptions::default(),
    )
}

/// [`newton_solve`] with explicit [`SolveOptions`].
pub fn newton_solve_with<F>(
    params: &mut [f64],
    fixed: &[bool],
    n_eq: usize,
    eval_residuals: F,
    options: &SolveOptions,
) -> SolveResult
where
    F: FnMut(&[f64], &mut [f64]),
{
    newton_solve_with_jacobian(
        params,
        fixed,
        n_eq,
        eval_residuals,
        Jacobian::FiniteDifference,
        options,
    )
}

/// Analytic Jacobian assembly callback: writes the derivatives at the given
/// parameters into the values of the [`SparseJacobian`] whose pattern it was
/// handed with (zeroed by the loop before each call).
pub(crate) type JacobianBuilder<'a> = dyn FnMut(&[f64], &mut SparseJacobian) + 'a;

/// Where the Newton loop sources its Jacobian each iteration.
pub(crate) enum Jacobian<'a> {
    /// Central finite differences over the residual closure (for callers
    /// without analytic derivatives), over a dense pattern.
    FiniteDifference,
    /// Analytic assembly (solver roadmap M2) into a fixed pattern, see
    /// [`JacobianBuilder`].
    Analytic {
        /// The entries every build may write: one row per equation, one
        /// column per free coordinate.
        pattern: SparseJacobian,
        build: &'a mut JacobianBuilder<'a>,
    },
}

/// A Jacobian in compressed sparse rows: per equation, the columns it may
/// have a derivative in (ascending) and the derivatives.
///
/// A constraint touches a handful of parameters, so a row has a handful of
/// entries however large the component is, and which ones is fixed for a
/// solve. The pattern is laid out once; each Jacobian rewrites the values.
/// An entry of the pattern may hold an exact zero (a derivative that
/// vanishes at this point); every product skips those, as a dense matrix
/// skipped its zeros, so the arithmetic is the same.
#[derive(Clone)]
pub(crate) struct SparseJacobian {
    n_cols: usize,
    /// Row `r` is entries `starts[r]..starts[r + 1]`.
    starts: Vec<usize>,
    cols: Vec<usize>,
    vals: Vec<f64>,
}

impl SparseJacobian {
    /// The pattern of rows `rows` (columns ascending, each listed once), over
    /// `n_cols` columns, with every value zero.
    pub(crate) fn from_rows<R: AsRef<[usize]>>(n_cols: usize, rows: &[R]) -> Self {
        let mut starts = Vec::with_capacity(rows.len() + 1);
        let mut cols = Vec::new();
        starts.push(0);
        for row in rows {
            cols.extend_from_slice(row.as_ref());
            starts.push(cols.len());
        }
        Self::from_parts(n_cols, starts, cols)
    }

    /// The pattern whose row `r` is `cols[starts[r]..starts[r + 1]]`
    /// (columns ascending, each listed once), with every value zero.
    pub(crate) fn from_parts(n_cols: usize, starts: Vec<usize>, cols: Vec<usize>) -> Self {
        debug_assert_eq!(starts.first(), Some(&0));
        debug_assert_eq!(starts.last(), Some(&cols.len()));
        debug_assert!(starts
            .windows(2)
            .all(|span| cols[span[0]..span[1]].windows(2).all(|w| w[0] < w[1])));
        debug_assert!(cols.iter().all(|&c| c < n_cols));
        let vals = vec![0.0; cols.len()];
        Self {
            n_cols,
            starts,
            cols,
            vals,
        }
    }

    /// Every column of every row: the pattern of a Jacobian whose entries
    /// are not known in advance (finite differences).
    pub(crate) fn dense(n_rows: usize, n_cols: usize) -> Self {
        let all: Vec<usize> = (0..n_cols).collect();
        Self::from_rows(n_cols, &vec![all; n_rows])
    }

    /// The non-zero entries of a dense matrix.
    #[cfg(test)]
    pub(crate) fn from_dense(j: &[Vec<f64>], n_cols: usize) -> Self {
        let rows: Vec<Vec<usize>> = j
            .iter()
            .map(|row| (0..n_cols).filter(|&c| row[c] != 0.0).collect())
            .collect();
        let mut sparse = Self::from_rows(n_cols, &rows);
        sparse.copy_from_dense(j);
        sparse
    }

    pub(crate) fn n_rows(&self) -> usize {
        self.starts.len() - 1
    }

    /// The first entry of row `r`: the slot of the `k`-th column of the
    /// row's pattern is `row_start(r) + k`.
    pub(crate) fn row_start(&self, r: usize) -> usize {
        self.starts[r]
    }

    /// The columns of row `r`'s pattern.
    pub(crate) fn row_cols(&self, r: usize) -> &[usize] {
        &self.cols[self.starts[r]..self.starts[r + 1]]
    }

    /// The values of row `r`, one per entry of [`SparseJacobian::row_cols`].
    pub(crate) fn row_values(&self, r: usize) -> &[f64] {
        &self.vals[self.starts[r]..self.starts[r + 1]]
    }

    /// The values, slot by slot (see [`SparseJacobian::row_start`]).
    pub(crate) fn values_mut(&mut self) -> &mut [f64] {
        &mut self.vals
    }

    /// The non-zero entries of row `r`, as `(column, value)`, columns
    /// ascending.
    fn row(&self, r: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let span = self.starts[r]..self.starts[r + 1];
        self.cols[span.clone()]
            .iter()
            .zip(&self.vals[span])
            .filter(|(_, v)| **v != 0.0)
            .map(|(c, v)| (*c, *v))
    }

    /// The non-zero entries of every row, row by row.
    fn rows(&self) -> impl Iterator<Item = impl Iterator<Item = (usize, f64)> + '_> + '_ {
        (0..self.n_rows()).map(|r| self.row(r))
    }

    fn clear_values(&mut self) {
        self.vals.fill(0.0);
    }

    /// Take the values of the pattern's entries from a dense matrix.
    fn copy_from_dense(&mut self, j: &[Vec<f64>]) {
        for (r, row) in j.iter().enumerate() {
            let span = self.starts[r]..self.starts[r + 1];
            for (v, &c) in self.vals[span.clone()].iter_mut().zip(&self.cols[span]) {
                *v = row[c];
            }
        }
    }

    /// The matrix with every row written out in full.
    #[cfg(debug_assertions)]
    pub(crate) fn to_dense(&self) -> Vec<Vec<f64>> {
        let mut dense = vec![vec![0.0_f64; self.n_cols]; self.n_rows()];
        for (r, row) in dense.iter_mut().enumerate() {
            let span = self.starts[r]..self.starts[r + 1];
            for (&c, &v) in self.cols[span.clone()].iter().zip(&self.vals[span]) {
                row[c] = v;
            }
        }
        dense
    }
}

/// How a step in the free coordinates is applied to the stored parameters.
///
/// The Newton loop works in free coordinates: one per Jacobian column. For
/// most parameters a coordinate *is* a stored parameter and a step is an
/// addition ([`AdditiveSpace`]). A parameter that lives on a manifold — the
/// orientation quaternion of a rigid body — has fewer coordinates than
/// stored numbers and is stepped by a retraction instead (rotate, then
/// renormalise). The loop only ever asks for "the saved state moved by this
/// step", so it does not need to know which kind it is driving.
pub(crate) trait StepSpace {
    /// Number of free coordinates (Jacobian columns).
    fn dim(&self) -> usize;

    /// The stored parameters a step may change, each listed once.
    fn storage(&self) -> &[usize];

    /// Write `saved ⊕ step` into the stored parameters: `saved` holds the
    /// values of [`StepSpace::storage`] the step starts from, `step` has
    /// [`StepSpace::dim`] entries.
    fn apply(&self, params: &mut [f64], saved: &[f64], step: &[f64]);

    /// The free parameters themselves, when every coordinate is a stored
    /// parameter and steps are additions.
    fn additive_params(&self) -> Option<&[usize]> {
        None
    }
}

/// The step space of plain parameters: coordinate `k` is the stored
/// parameter `free[k]`, and a step adds to it.
pub(crate) struct AdditiveSpace {
    free: Vec<usize>,
}

impl AdditiveSpace {
    /// `free` lists the parameters that may move.
    pub(crate) fn new(free: Vec<usize>) -> Self {
        Self { free }
    }

    /// Every parameter not held by the `fixed` mask.
    pub(crate) fn from_mask(fixed: &[bool]) -> Self {
        Self::new((0..fixed.len()).filter(|&i| !fixed[i]).collect())
    }
}

impl StepSpace for AdditiveSpace {
    fn dim(&self) -> usize {
        self.free.len()
    }

    fn storage(&self) -> &[usize] {
        &self.free
    }

    fn apply(&self, params: &mut [f64], saved: &[f64], step: &[f64]) {
        for (k, &i) in self.free.iter().enumerate() {
            params[i] = saved[k] + step[k];
        }
    }

    fn additive_params(&self) -> Option<&[usize]> {
        Some(&self.free)
    }
}

/// The values of the stored parameters `storage` in `params`.
fn save(params: &[f64], storage: &[usize]) -> Vec<f64> {
    storage.iter().map(|&i| params[i]).collect()
}

/// Put `saved` back into the stored parameters `storage`.
fn restore(params: &mut [f64], storage: &[usize], saved: &[f64]) {
    for (&i, &v) in storage.iter().zip(saved) {
        params[i] = v;
    }
}

/// The damped Newton loop over the parameters not held by `fixed`, each
/// stepped additively.
pub(crate) fn newton_solve_with_jacobian<F>(
    params: &mut [f64],
    fixed: &[bool],
    n_eq: usize,
    eval_residuals: F,
    jacobian: Jacobian<'_>,
    options: &SolveOptions,
) -> SolveResult
where
    F: FnMut(&[f64], &mut [f64]),
{
    let space = AdditiveSpace::from_mask(&fixed[..params.len()]);
    newton_solve_in(
        params,
        &space,
        n_eq,
        eval_residuals,
        jacobian,
        None,
        &CoreOptions::from(options),
    )
}

/// How large the numbers of a system are, for judging what a rounding of
/// them is worth (see [`PRECISION_ULPS`]).
pub(crate) struct Magnitudes<'a> {
    /// Per free coordinate of the step space: the size of the quantity it
    /// moves.
    pub(crate) columns: &'a mut dyn FnMut(&[f64], &mut [f64]),
    /// Per equation: `Σ|∂f/∂q|·|q|` over everything the equation touches,
    /// held or free.
    pub(crate) rows: &'a mut dyn FnMut(&[f64], &mut [f64]),
    /// The unit the equations are measured in, at most one: the tolerance
    /// applies to `residual / unit`. A family whose equations are all
    /// lengths passes the size of a component smaller than one unit, so
    /// that the tolerance is relative to it; one passes 1 otherwise.
    pub(crate) unit: f64,
}

/// The weight of an equation whose terms have magnitude `size`: `unit`,
/// unless the tolerance is finer than those terms resolve, in which case
/// the equation is measured in units of that resolution instead. One
/// rounding of the terms is worth `ε·size`; the weight is
/// `PRECISION_ULPS·ε·size / tolerance` once that exceeds `unit`.
pub(crate) fn weight_of(size: f64, tolerance: f64, unit: f64) -> f64 {
    let weight = PRECISION_ULPS * f64::EPSILON * size / tolerance;
    // A tolerance of zero, or a size that is not a number, weighs nothing.
    if weight.is_finite() {
        weight.max(unit)
    } else {
        unit
    }
}

/// The unit of [`Magnitudes`], if any: one for a family that has none.
fn unit_of(magnitudes: &Option<Magnitudes<'_>>) -> f64 {
    magnitudes.as_ref().map_or(1.0, |m| m.unit)
}

/// Move residuals weighted by `old` onto the weights `new`.
fn reweigh(f: &mut [f64], old: &[f64], new: &[f64]) {
    for ((fi, old), new) in f.iter_mut().zip(old).zip(new) {
        *fi = *fi * old / new;
    }
}

/// The weights from the Jacobian the loop has just built: `size` is taken
/// over the free columns only, from rows the loop needs anyway, so an
/// ordinary solve pays nothing for it. (What the held columns add is found
/// by [`fresh_weights`] if the solve turns out to need it.) Returns
/// whether any weight is not one.
fn weights_from_rows(
    weights: &mut [f64],
    j: &SparseJacobian,
    column: &[f64],
    tolerance: f64,
    unit: f64,
) -> bool {
    let mut any = false;
    for (weight, row) in weights.iter_mut().zip(j.rows()) {
        let size: f64 = row.map(|(col, v)| v.abs() * column[col]).sum();
        *weight = weight_of(size, tolerance, unit);
        any |= *weight != 1.0;
    }
    any
}

/// The weights at `params`, from the full magnitudes of the `n_eq`
/// equations.
fn weights_at(
    n_eq: usize,
    params: &[f64],
    magnitudes: &mut Magnitudes<'_>,
    tolerance: f64,
) -> Vec<f64> {
    let mut size = vec![0.0_f64; n_eq];
    (magnitudes.rows)(params, &mut size);
    size.iter()
        .map(|s| weight_of(*s, tolerance, magnitudes.unit))
        .collect()
}

/// The weights at `params` ([`weights_at`]) — if they differ enough from
/// `weights` to be worth another pass: some weight at least twice what it
/// was.
fn fresh_weights(
    weights: &[f64],
    params: &[f64],
    magnitudes: &mut Option<Magnitudes<'_>>,
    tolerance: f64,
) -> Option<Vec<f64>> {
    let magnitudes = magnitudes.as_mut()?;
    let fresh = weights_at(weights.len(), params, magnitudes, tolerance);
    fresh
        .iter()
        .zip(weights)
        .any(|(new, old)| *new >= 2.0 * old)
        .then_some(fresh)
}

/// Divide each residual by the weight of its equation.
fn weigh(f: &mut [f64], weights: &[f64]) {
    for (fi, weight) in f.iter_mut().zip(weights) {
        *fi /= weight;
    }
}

/// Divide each Jacobian row by the weight of its equation.
fn weigh_rows(j: &mut SparseJacobian, weights: &[f64]) {
    for (r, weight) in weights.iter().enumerate() {
        if *weight != 1.0 {
            let span = j.starts[r]..j.starts[r + 1];
            for v in j.vals[span].iter_mut() {
                *v /= weight;
            }
        }
    }
}

/// Per column, the scale the step is computed in when the equations are
/// weighted.
///
/// Dividing a row by a weight of 1e7 divides its share of `JᵀJ` by 1e14. A
/// column that only such rows touch is left with entries of that size,
/// while the regularisation `λ` and the pivot threshold of the elimination
/// are small absolute numbers sized for entries of order one: the column
/// would read as singular, or `λ` would swallow it and its step vanish. So
/// each column is scaled back up by the smallest weight among the rows
/// that touch it ([`scale_columns`]), the step is found in those scaled
/// coordinates, and scaled back before it is applied. A column touched by
/// an unweighted row has scale one.
fn column_scales(j: &SparseJacobian, weights: &[f64], n: usize) -> Vec<f64> {
    let mut scales = vec![f64::INFINITY; n];
    for (row, weight) in j.rows().zip(weights) {
        for (col, _) in row {
            if *weight < scales[col] {
                scales[col] = *weight;
            }
        }
    }
    // A column no row touches is left alone.
    for scale in scales.iter_mut() {
        if !scale.is_finite() {
            *scale = 1.0;
        }
    }
    scales
}

/// Multiply each column of the Jacobian by its scale: the Jacobian with
/// respect to the scaled coordinates `y`, where the step is `h = scale·y`.
fn scale_columns(j: &mut SparseJacobian, scales: &[f64]) {
    for (v, &col) in j.vals.iter_mut().zip(&j.cols) {
        *v *= scales[col];
    }
}

/// The step in the parameters for a step `y` in scaled coordinates.
fn unscaled(y: &[f64], scales: &[f64]) -> Vec<f64> {
    y.iter().zip(scales).map(|(y, scale)| y * scale).collect()
}

/// RMS of the residuals as the equations state them, from the weighted
/// ones.
fn unweighted_rms(f: &[f64], weights: &[f64]) -> f64 {
    let sum: f64 = f
        .iter()
        .zip(weights)
        .map(|(fi, weight)| (fi * weight) * (fi * weight))
        .sum();
    (sum / f.len().max(1) as f64).sqrt()
}

/// The trust-region Newton loop over the free coordinates of `space`, with
/// an explicit Jacobian source.
///
/// With `magnitudes`, every equation is weighed by what its numbers can
/// resolve (see [`weight_of`]). Near the origin every weight is one, and
/// the weighting is skipped altogether. Far from it, an equation on
/// coordinates of a billion cannot be met to 1e-10, and neither should the
/// least-squares step trade the rounding of such an equation against one
/// that is exact: the weights put both in the same units, and a solve whose
/// equations are all as exact as their numbers allow ends as
/// [`SolveStatus::ConvergedToPrecision`].
pub(crate) fn newton_solve_in<F>(
    params: &mut [f64],
    space: &dyn StepSpace,
    n_eq: usize,
    mut eval_residuals: F,
    jacobian: Jacobian<'_>,
    mut magnitudes: Option<Magnitudes<'_>>,
    options: &CoreOptions,
) -> SolveResult
where
    F: FnMut(&[f64], &mut [f64]),
{
    let n_free = space.dim();
    let storage = space.storage();

    if n_eq == 0 {
        return SolveResult {
            converged: true,
            residual: 0.0,
            iterations: 0,
            status: SolveStatus::Converged,
        };
    }

    // The loop measures the RMS residual. An L2 tolerance over `n_eq`
    // equations is the same test with the tolerance scaled.
    let tolerance = match options.norm {
        ResidualNorm::Rms => options.tolerance,
        ResidualNorm::L2 => options.tolerance / (n_eq as f64).sqrt(),
    };

    // One per equation; `weighted` says whether any is not one. Until the
    // first Jacobian gives the precision of each equation, every weight is
    // the unit of the equations (see `Magnitudes::unit`).
    let unit = unit_of(&magnitudes);
    let mut weights = vec![unit; n_eq];
    let mut weighted = unit != 1.0;

    if n_free == 0 {
        // Equations but nothing that may move: the state is what it is,
        // satisfied or not. Report which (0.1.0 called both converged).
        let f = eval_vec(&mut eval_residuals, params, n_eq);
        let residual = rms(&f);
        // Judged in the units of the equations and as far as their numbers
        // go: met as stated, met to the precision of the numbers, or unmet.
        let mut weighed = f;
        if let Some(magnitudes) = magnitudes.as_mut() {
            weigh(
                &mut weighed,
                &weights_at(n_eq, params, magnitudes, tolerance),
            );
        }
        let judged = rms(&weighed);
        let status = if judged < tolerance {
            if residual < tolerance {
                SolveStatus::Converged
            } else {
                SolveStatus::ConvergedToPrecision
            }
        } else if !residual.is_finite() {
            SolveStatus::Diverged
        } else {
            SolveStatus::Stalled
        };
        return SolveResult {
            converged: status <= SolveStatus::ConvergedToPrecision,
            residual,
            iterations: 0,
            status,
        };
    }

    // The state a failed solve is rolled back to.
    let initial = save(params, storage);
    // Why the loop stopped short of convergence. Only reported when the
    // final residual is not under the tolerance.
    let mut failure = SolveStatus::MaxIterations;

    let mut residual;
    let mut iters = 0;
    // The Jacobian source is matched once (it is consumed by value); the
    // per-iteration branch checks the option. The pattern is laid out here
    // and its values rewritten by every Jacobian (the analytic builder
    // receives them zeroed; differences are taken into a dense scratch
    // matrix and copied over).
    let (mut j, mut analytic, mut fd_scratch) = match jacobian {
        Jacobian::FiniteDifference => (
            SparseJacobian::dense(n_eq, n_free),
            None,
            vec![vec![0.0_f64; n_free]; n_eq],
        ),
        Jacobian::Analytic { pattern, build } => (pattern, Some(build), Vec::new()),
    };
    debug_assert_eq!(j.n_rows(), n_eq);
    debug_assert_eq!(j.n_cols, n_free);

    // The residual vector at the current iterate (weighted, once there are
    // weights). An accepted trial step hands its residuals over, so each
    // point is evaluated once.
    let mut f = eval_vec(&mut eval_residuals, params, n_eq);
    if weighted {
        weigh(&mut f, &weights);
    }

    // The trust radius, in the norm of the step vector. Unbounded until a
    // step has been judged poor.
    let mut radius = f64::INFINITY;

    // The weights are taken from the first Jacobian and then held while the
    // loop runs, so that the residuals it compares are in the same units.
    // If it ends without converging, the geometry may have travelled to
    // where its numbers are larger than they were at the start, or the
    // held parameters may be the large ones; the weights are then taken
    // afresh, in full, and the loop goes on from where it stands.
    let mut primed = false;
    let mut reweighs = 0;
    loop {
        while iters < options.max_iterations {
            residual = rms(&f);
            if residual < tolerance {
                break;
            }
            if !residual.is_finite() {
                failure = SolveStatus::Diverged;
                break;
            }

            // Jacobian J (n_eq × n_free).
            if let Some(build) = analytic.as_mut() {
                j.clear_values();
                build(params, &mut j);
            } else {
                numerical_jacobian_for(params, space, n_eq, &mut eval_residuals, &mut fd_scratch);
                j.copy_from_dense(&fd_scratch);
            }

            if !primed {
                primed = true;
                if let Some(magnitudes) = magnitudes.as_mut() {
                    let mut column = vec![0.0_f64; n_free];
                    (magnitudes.columns)(params, &mut column);
                    let before = weighted.then(|| weights.clone());
                    weighted = weights_from_rows(&mut weights, &j, &column, tolerance, unit);
                    match &before {
                        Some(before) => reweigh(&mut f, before, &weights),
                        None if weighted => weigh(&mut f, &weights),
                        None => {}
                    }
                    if rms(&f) < tolerance {
                        break;
                    }
                }
            }
            let scales = weighted.then(|| column_scales(&j, &weights, n_free));
            if let Some(scales) = &scales {
                weigh_rows(&mut j, &weights);
                scale_columns(&mut j, scales);
            }

            // The linear model f + J·h at this iterate, and a dogleg step
            // inside the trust radius. A trial that gains nothing is undone
            // and the radius halved; the next one is shorter and closer to
            // the steepest-descent direction. The model is not rebuilt for
            // it. (With weights, the model, its steps and the radius are in
            // the scaled coordinates.)
            let saved = save(params, storage);
            let model = LinearModel::new(&j, &f, n_free, options);
            // The residuals of the trial step that was accepted, if one was.
            let mut accepted: Option<Vec<f64>> = None;
            for _ in 0..TRUST_ATTEMPTS {
                let Some(step) = model.step(radius) else {
                    break;
                };
                let length = norm(&step);
                match &scales {
                    None => space.apply(params, &saved, &step),
                    Some(scales) => space.apply(params, &saved, &unscaled(&step, scales)),
                }
                let mut f_new = eval_vec(&mut eval_residuals, params, n_eq);
                if weighted {
                    weigh(&mut f_new, &weights);
                }
                let gain = model.gain(&step, &f_new);
                if gt(gain, TRUST_WIDEN_ABOVE) {
                    radius = radius.max(TRUST_GROWTH * length);
                } else if !gt(gain, TRUST_NARROW_BELOW) {
                    radius = 0.5 * radius.min(length);
                }
                if gt(gain, 0.0) {
                    accepted = Some(f_new);
                    break;
                }
            }
            let stuck = model.failure();
            match accepted {
                Some(f_new) => f = f_new,
                None => {
                    // No step of any length lowered the residual: put the
                    // iterate back and stop — the next iteration would
                    // rebuild the same model.
                    restore(params, storage, &saved);
                    failure = stuck;
                    break;
                }
            }

            iters += 1;
        }

        residual = rms(&f);
        if residual.is_finite() && residual >= tolerance && reweighs < REWEIGH_PASSES {
            if let Some(fresh) = fresh_weights(&weights, params, &mut magnitudes, tolerance) {
                reweighs += 1;
                weights = fresh;
                weighted = true;
                primed = true;
                f = eval_vec(&mut eval_residuals, params, n_eq);
                weigh(&mut f, &weights);
                radius = f64::INFINITY;
                failure = SolveStatus::MaxIterations;
                continue;
            }
        }
        break;
    }

    // `f` belongs to the iterate the loop left: an accepted step handed its
    // residuals over, and a rejected one was put back.
    let mut f_final = f;
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
    //
    // The Jacobian is built and the normal equations factored once, at the
    // first polish step; later steps reuse them (a chord iteration). The
    // iterate has already converged, so the Jacobian barely moves between
    // polish steps, and each further step costs a back-substitution
    // instead of a Jacobian and a factorisation: the polish used to be two
    // thirds of the factorisations of an easy solve.
    let mut polish_steps = 0usize;
    let mut chord: Option<Chord> = None;
    while residual > 0.0 && residual < tolerance && polish_steps < options.polish_steps {
        if chord.is_none() {
            if let Some(build) = analytic.as_mut() {
                j.clear_values();
                build(params, &mut j);
            } else {
                numerical_jacobian_for(params, space, n_eq, &mut eval_residuals, &mut fd_scratch);
                j.copy_from_dense(&fd_scratch);
            }
            let scales = weighted.then(|| column_scales(&j, &weights, n_free));
            if let Some(scales) = &scales {
                weigh_rows(&mut j, &weights);
                scale_columns(&mut j, scales);
            }
            // Unregularised step: see the block comment above. A singular
            // system has none; a family that asks for it gets the
            // regularised step there instead, so an under-constrained
            // component is polished too (toward its minimum-norm
            // correction). Otherwise there is no polish.
            let (normal, b) = j.normal_equations(&f_final);
            // The matrix is copied only where a regularised retry may
            // need it.
            let factors = if options.polish_singular {
                LuFactors::factor(normal.clone()).or_else(|| {
                    let mut normal = normal;
                    normal.add_diagonal(options.lambda);
                    LuFactors::factor(normal)
                })
            } else {
                LuFactors::factor(normal)
            };
            let Some(factors) = factors else {
                break;
            };
            chord = Some(Chord {
                factors,
                scales,
                first: Some(b),
            });
        }
        let Some(Chord {
            factors,
            scales,
            first,
        }) = &mut chord
        else {
            break;
        };
        // b = Jᵀ(−f). The first step has it from the normal equations; a
        // later one forms it over the same Jacobian (still in `j`), in the
        // order they sum it.
        let mut b = match first.take() {
            Some(b) => b,
            None => {
                let mut b = vec![0.0_f64; n_free];
                for (row, fi) in j.rows().zip(&f_final) {
                    for (col, v) in row {
                        b[col] -= v * fi;
                    }
                }
                b
            }
        };
        let delta = factors.solve(&mut b);
        if delta.iter().any(|d| !d.is_finite()) {
            break;
        }
        let delta = match scales {
            None => delta,
            Some(scales) => unscaled(&delta, scales),
        };
        let saved = save(params, storage);
        space.apply(params, &saved, &delta);
        let mut f_polish = eval_vec(&mut eval_residuals, params, n_eq);
        if weighted {
            weigh(&mut f_polish, &weights);
        }
        let residual_polish = rms(&f_polish);
        if residual_polish < residual {
            residual = residual_polish;
            f_final = f_polish;
            polish_steps += 1;
        } else {
            // Worsened (or NaN, or no strict gain) — reject the polish
            // step, restore, and stop.
            restore(params, storage, &saved);
            break;
        }
    }

    // What is reported is the residual of the equations as stated.
    let stated = if weighted {
        unweighted_rms(&f_final, &weights)
    } else {
        residual
    };
    let status = if residual < tolerance {
        if stated < tolerance {
            SolveStatus::Converged
        } else {
            SolveStatus::ConvergedToPrecision
        }
    } else if !residual.is_finite() {
        SolveStatus::Diverged
    } else {
        failure
    };
    let converged = status <= SolveStatus::ConvergedToPrecision;
    residual = stated;
    if !converged && options.on_failure == FailurePolicy::Restore {
        restore(params, storage, &initial);
        residual = rms(&eval_vec(&mut eval_residuals, params, n_eq));
    }

    SolveResult {
        converged,
        residual,
        iterations: iters,
        status,
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

/// What the convergence polish keeps from its first step.
struct Chord {
    /// The factored normal matrix.
    factors: LuFactors,
    /// The column scales, when the equations are weighted.
    scales: Option<Vec<f64>>,
    /// The right-hand side of the first step, until it is taken.
    first: Option<Vec<f64>>,
}

/// Columns whose remaining norm falls to this fraction of the largest
/// column norm or below end the rank of a [`min_norm_step`]: such a
/// direction is too weakly determined to step along. (The Tikhonov term of
/// the normal equations did the same job less sharply, by damping every
/// direction whose squared singular value is near `λ`.)
pub(crate) const STEP_RCOND: f64 = 1e-12;

/// When the smallest pivot of the factored normal equations falls below
/// this fraction of the largest, the Gauss–Newton step is taken from
/// [`min_norm_step`] instead. The normal equations square the condition
/// number, so at this ratio they have lost about half the digits; and in a
/// direction whose pivot is near `λ` the regularisation, not the geometry,
/// sets the step (rank deficiency, or a direction as weakly determined as
/// the side of a needle-thin triangle).
pub(crate) const QR_TRIGGER: f64 = 1e-8;

/// The Gauss–Newton step `h` of `J·h ≈ −f`: from the regularised normal
/// equations `(JᵀJ + λI)·h = −Jᵀf` while they are well conditioned (the
/// fast path, bitwise what elimination gives), otherwise from a
/// rank-truncated QR of `J` (see [`QR_TRIGGER`]).
fn gauss_newton_step(j: &SparseJacobian, f: &[f64], lambda: f64) -> Option<Vec<f64>> {
    let (mut a, mut b) = j.normal_equations(f);
    a.add_diagonal(lambda);
    match LuFactors::factor(a) {
        Some(factors) if factors.pivot_ratio() >= QR_TRIGGER => Some(factors.solve(&mut b)),
        _ => min_norm_step(j, f),
    }
}

/// Apply the Householder reflector `I − 2·v·vᵀ / (vᵀv)` to `x`.
fn reflect(x: &mut [f64], v: &[f64], v_norm2: f64) {
    let beta = 2.0 * dot(x, v) / v_norm2;
    for (xi, vi) in x.iter_mut().zip(v) {
        *xi -= beta * vi;
    }
}

/// The reflector that takes `x` to `(α, 0, …, 0)`: returns `(v, vᵀv, α)`,
/// with the sign of `α` chosen against `x[0]` so that no cancellation
/// occurs. `None` for a zero `x`.
fn householder(x: &[f64]) -> Option<(Vec<f64>, f64, f64)> {
    let sigma = norm(x);
    if !gt(sigma, 0.0) {
        return None;
    }
    let alpha = if x[0] > 0.0 { -sigma } else { sigma };
    let mut v = x.to_vec();
    v[0] -= alpha;
    let v_norm2 = dot(&v, &v);
    Some((v, v_norm2, alpha))
}

/// The minimum-norm least-squares solution of `J·h = −f`, by a complete
/// orthogonal decomposition (as LAPACK's `xGELSY`): a column-pivoted
/// Householder QR of `J` decides the rank — columns are taken while the
/// remaining norm exceeds [`STEP_RCOND`] times the largest — and, when that
/// rank is short of the column count, a second QR of the leading rows of
/// `R` gives the shortest of the solutions. Nothing is squared, so the step
/// is as accurate as `κ(J)` allows rather than `κ(J)²`.
///
/// `None` when `J` has no column above the threshold. A `J` or `f` that is
/// not finite gives a step of NaN, for the caller to report.
pub(crate) fn min_norm_step(j: &SparseJacobian, f: &[f64]) -> Option<Vec<f64>> {
    let (m, n) = (f.len(), j.n_cols);
    if f.iter().any(|v| !v.is_finite()) || j.vals.iter().any(|v| !v.is_finite()) {
        return Some(vec![f64::NAN; n]);
    }
    // The columns of J, so that the reflectors run over contiguous memory.
    let mut cols = vec![vec![0.0_f64; m]; n];
    for (i, row) in j.rows().enumerate() {
        for (c, v) in row {
            cols[c][i] += v;
        }
    }
    let mut b: Vec<f64> = f.iter().map(|v| -v).collect();
    let mut perm: Vec<usize> = (0..n).collect();
    // Remaining column norms, downdated as rows are eliminated and taken
    // afresh where the downdate has cancelled too far (as LAPACK does).
    let mut norms: Vec<f64> = cols.iter().map(|c| norm(c)).collect();
    let mut exact = norms.clone();
    let largest = norms.iter().fold(0.0_f64, |a, &b| a.max(b));
    if !gt(largest, 0.0) {
        return None;
    }
    let mut rank = 0;
    for k in 0..m.min(n) {
        let p = (k..n)
            .max_by(|&a, &c| norms[a].total_cmp(&norms[c]))
            .expect("k < n");
        cols.swap(k, p);
        perm.swap(k, p);
        norms.swap(k, p);
        exact.swap(k, p);
        let Some((v, v_norm2, alpha)) = householder(&cols[k][k..]) else {
            break;
        };
        if !gt(alpha.abs(), STEP_RCOND * largest) {
            break;
        }
        for col in cols[k + 1..].iter_mut() {
            reflect(&mut col[k..], &v, v_norm2);
        }
        reflect(&mut b[k..], &v, v_norm2);
        cols[k][k] = alpha;
        for x in cols[k][k + 1..].iter_mut() {
            *x = 0.0;
        }
        for j in k + 1..n {
            if norms[j] > 0.0 {
                let t = cols[j][k] / norms[j];
                let shrink = (1.0 - t * t).max(0.0);
                if shrink * (norms[j] / exact[j]).powi(2) <= f64::EPSILON.sqrt() {
                    norms[j] = norm(&cols[j][k + 1..]);
                    exact[j] = norms[j];
                } else {
                    norms[j] *= shrink.sqrt();
                }
            }
        }
        rank += 1;
    }
    if rank == 0 {
        return None;
    }
    // R = [R11 R12] in the first `rank` rows: R[i][j] = cols[j][i], j ≥ i.
    let c = &b[..rank];
    let mut y = vec![0.0_f64; n];
    if rank == n {
        for i in (0..n).rev() {
            let mut s = c[i];
            for j in i + 1..n {
                s -= cols[j][i] * y[j];
            }
            y[i] = s / cols[i][i];
        }
    } else {
        // Rᵀ (n × rank) = Z·[T; 0]: the columns of Rᵀ are the rows of R.
        let mut rt: Vec<Vec<f64>> = (0..rank)
            .map(|i| {
                (0..n)
                    .map(|j| if j >= i { cols[j][i] } else { 0.0 })
                    .collect()
            })
            .collect();
        let mut reflectors = Vec::with_capacity(rank);
        for i in 0..rank {
            let (v, v_norm2, alpha) = householder(&rt[i][i..])?;
            for col in rt[i + 1..].iter_mut() {
                reflect(&mut col[i..], &v, v_norm2);
            }
            rt[i][i] = alpha;
            reflectors.push((v, v_norm2));
        }
        // R = [Tᵀ 0]·Zᵀ, so R·y = c is Tᵀ·w = c with y = Z·[w; 0]; T[j][i]
        // (j ≤ i) is rt[i][j].
        let mut w = vec![0.0_f64; rank];
        for i in 0..rank {
            let mut s = c[i];
            for j in 0..i {
                s -= rt[i][j] * w[j];
            }
            w[i] = s / rt[i][i];
        }
        y[..rank].copy_from_slice(&w);
        for (i, (v, v_norm2)) in reflectors.iter().enumerate().rev() {
            reflect(&mut y[i..], v, *v_norm2);
        }
    }
    let mut h = vec![0.0_f64; n];
    for (j, &p) in perm.iter().enumerate() {
        h[p] = y[j];
    }
    Some(h)
}

/// The linear model `f + J·h` of the residuals around one iterate, with the
/// two steps a dogleg is drawn between.
struct LinearModel<'a> {
    j: &'a SparseJacobian,
    f: &'a [f64],
    /// `‖f‖²`.
    cost: f64,
    /// `Jᵀf`: the squared residual falls fastest along `−g`.
    g: Vec<f64>,
    /// The Gauss–Newton step: the regularised normal equations, scaled by
    /// the caller's step factor. `None` when they are numerically singular.
    newton: Option<Vec<f64>>,
    /// The Cauchy step: where the model is smallest along `−g`. `None` at a
    /// stationary point, where `g` vanishes.
    cauchy: Option<Vec<f64>>,
    /// Whether the normal equations gave something that is not a number.
    non_finite: bool,
}

impl<'a> LinearModel<'a> {
    fn new(j: &'a SparseJacobian, f: &'a [f64], n_free: usize, options: &CoreOptions) -> Self {
        let mut g = vec![0.0_f64; n_free];
        for (row, fi) in j.rows().zip(f) {
            for (col, v) in row {
                g[col] += v * fi;
            }
        }
        let g_sq = dot(&g, &g);
        let jg_sq: f64 = j
            .rows()
            .map(|row| {
                let jg: f64 = row.map(|(col, v)| v * g[col]).sum();
                jg * jg
            })
            .sum();
        let cauchy = (gt(g_sq, 0.0) && gt(jg_sq, 0.0)).then(|| {
            let alpha = g_sq / jg_sq;
            g.iter().map(|gk| -alpha * gk).collect()
        });

        let solved = gauss_newton_step(j, f, options.lambda);
        let non_finite = solved
            .as_ref()
            .is_some_and(|h| h.iter().any(|d| !d.is_finite()));
        let newton = solved
            .filter(|_| !non_finite)
            .map(|h| h.iter().map(|d| options.initial_step * d).collect());

        Self {
            j,
            f,
            cost: dot(f, f),
            g,
            newton,
            cauchy,
            non_finite,
        }
    }

    /// The dogleg step of length at most `radius`: the Gauss–Newton step if
    /// it fits; otherwise the point where the path from the Cauchy step to
    /// the Gauss–Newton step leaves the radius, or, if even the Cauchy step
    /// is too long, a step of length `radius` straight down the gradient.
    /// `None` when there is no step to take.
    fn step(&self, radius: f64) -> Option<Vec<f64>> {
        let step: Vec<f64> = match (&self.newton, &self.cauchy) {
            (Some(newton), _) if norm(newton) <= radius => newton.clone(),
            (_, Some(cauchy)) if norm(cauchy) >= radius => {
                let scale = radius / norm(&self.g);
                self.g.iter().map(|gk| -scale * gk).collect()
            }
            (Some(newton), Some(cauchy)) => {
                // cauchy + β (newton − cauchy), with β such that the length
                // is the radius.
                let across: Vec<f64> = newton.iter().zip(cauchy).map(|(b, a)| b - a).collect();
                let (dd, ad, aa) = (
                    dot(&across, &across),
                    dot(cauchy, &across),
                    dot(cauchy, cauchy),
                );
                let root = (ad * ad + dd * (radius * radius - aa)).max(0.0).sqrt();
                let beta = if ad <= 0.0 {
                    (-ad + root) / dd
                } else {
                    (radius * radius - aa) / (ad + root)
                };
                cauchy
                    .iter()
                    .zip(&across)
                    .map(|(a, d)| a + beta * d)
                    .collect()
            }
            (Some(newton), None) => {
                let scale = radius / norm(newton);
                newton.iter().map(|b| scale * b).collect()
            }
            (None, Some(cauchy)) => cauchy.clone(),
            (None, None) => return None,
        };
        (step.iter().all(|d| d.is_finite()) && step.iter().any(|&d| d != 0.0)).then_some(step)
    }

    /// The gain of a step: the reduction of the squared residual it
    /// achieved, over the reduction `‖f‖² − ‖f + J·h‖²` the model promised.
    /// Negative when the model promised nothing or a number went missing.
    fn gain(&self, step: &[f64], f_new: &[f64]) -> f64 {
        let modelled: f64 = self
            .j
            .rows()
            .zip(self.f)
            .map(|(row, fi)| {
                let jh: f64 = row.map(|(col, v)| v * step[col]).sum();
                let r = fi + jh;
                r * r
            })
            .sum();
        let promised = self.cost - modelled;
        let achieved = self.cost - dot(f_new, f_new);
        if gt(promised, 0.0) && achieved.is_finite() {
            achieved / promised
        } else {
            -1.0
        }
    }

    /// Why no step could be accepted from this iterate.
    fn failure(&self) -> SolveStatus {
        if self.non_finite {
            SolveStatus::Diverged
        } else if self.newton.is_none() && self.cauchy.is_none() && gt(dot(&self.g, &self.g), 0.0) {
            // A gradient, and still no step: the linear algebra failed.
            SolveStatus::Singular
        } else {
            SolveStatus::Stalled
        }
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

/// Jacobian by central finite differences over the free coordinates of
/// `space`, written into the caller's buffer. `params` is perturbed in place
/// and put back.
pub(crate) fn numerical_jacobian_for<F>(
    params: &mut [f64],
    space: &dyn StepSpace,
    n_eq: usize,
    eval_residuals: &mut F,
    out: &mut [Vec<f64>],
) where
    F: FnMut(&[f64], &mut [f64]),
{
    if let Some(free) = space.additive_params() {
        numerical_jacobian_into(params, free, n_eq, eval_residuals, out);
        return;
    }

    let storage = space.storage();
    let saved = save(params, storage);
    let mut step = vec![0.0_f64; space.dim()];
    for col in 0..space.dim() {
        step[col] = FD_H;
        space.apply(params, &saved, &step);
        let f_plus = eval_vec(eval_residuals, params, n_eq);

        step[col] = -FD_H;
        space.apply(params, &saved, &step);
        let f_minus = eval_vec(eval_residuals, params, n_eq);

        step[col] = 0.0;
        for (row, o) in out.iter_mut().enumerate() {
            o[col] = (f_plus[row] - f_minus[row]) / (2.0 * FD_H);
        }
    }
    restore(params, storage, &saved);
}

/// The spacing of floats near a parameter, in units of which a
/// finite-difference step is kept at least this wide.
const FD_MIN_ULPS: f64 = 4096.0;

/// The finite-difference step for a parameter of value `x`: [`FD_H`], widened
/// for a parameter so large that `FD_H` is only a few float spacings there
/// (beyond about 1e5), where a difference taken with it would be rounding
/// noise. The step is deliberately not made relative to `x`: geometry far
/// from the origin is no larger for it, and a step that grows with the
/// offset soon exceeds the features it is meant to resolve. Nor does it
/// shrink for small geometry, which it cannot know about: a caller whose
/// features are below about 1e-4 needs analytic rows.
pub(crate) fn fd_step(x: f64) -> f64 {
    FD_H.max(FD_MIN_ULPS * f64::EPSILON * x.abs())
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
        let h = fd_step(orig);

        p[gi] = orig + h;
        let above = p[gi];
        let f_plus = eval_vec(eval_residuals, &p, n_eq);

        p[gi] = orig - h;
        let below = p[gi];
        let f_minus = eval_vec(eval_residuals, &p, n_eq);

        p[gi] = orig;

        // Divide by the step that was actually taken.
        let taken = above - below;
        for (row, o) in out.iter_mut().enumerate() {
            o[col] = (f_plus[row] - f_minus[row]) / taken;
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
    if !gt(max0, 0.0) {
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
///
/// Dense reference implementation, compiled for tests only: production code
/// uses the sparse variant in [`crate::mobility`].
#[cfg(test)]
pub(crate) fn null_space_support(j: &[Vec<f64>], rel_tol: f64) -> Vec<bool> {
    let n_cols = j.first().map_or(0, |row| row.len());
    let mut movable = vec![false; n_cols];
    for v in null_space_basis(j, rel_tol) {
        for (m, x) in movable.iter_mut().zip(&v) {
            *m |= x.abs() > rel_tol;
        }
    }
    movable
}

/// A basis of the null space of `j` (rows of equal length), from the same
/// reduced row echelon form as [`null_space_support`]: one vector per free
/// column, 1 there and the negated RREF entry at every pivot column.
/// Entries at or below `rel_tol` (relative to the equilibrated rows) are
/// left as computed; an empty `j` has every column free.
pub(crate) fn null_space_basis(j: &[Vec<f64>], rel_tol: f64) -> Vec<Vec<f64>> {
    let n_cols = j.first().map_or(0, |row| row.len());
    if n_cols == 0 {
        return Vec::new();
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
        // The row with the largest entry in the pivot column becomes the
        // pivot row (the column was chosen over all the rows left; the row
        // standing at the pivot position may have a zero there).
        let best_row = (pivot_row..n_rows)
            .max_by(|&a, &b| r[a][pcol].abs().total_cmp(&r[b][pcol].abs()))
            .expect("rows remain");
        r.swap(pivot_row, best_row);
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
    remaining.sort_unstable();
    remaining
        .iter()
        .map(|&free| {
            let mut v = vec![0.0_f64; n_cols];
            v[free] = 1.0;
            for (row, &pivot) in pivot_cols.iter().enumerate() {
                v[pivot] = -r[row][free];
            }
            v
        })
        .collect()
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
    fn a_pivot_column_missing_from_the_pivot_row_is_found_below() {
        // Column 1 has the largest norm but no entry in row 0: x0 pinned,
        // then x1 − x0 + x2 = 0 and x1 = 0 pin x1 and x2 too.
        let j = vec![
            vec![1.0, 0.0, 0.0],
            vec![-1.0, 2.0, 1.0],
            vec![0.0, 2.0, 0.0],
        ];
        assert_eq!(null_space_support(&j, 1e-9), vec![false, false, false]);
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
    try_solve_normal_equations(j, f, n, lambda).unwrap_or_else(|| vec![0.0; n])
}

/// [`solve_normal_equations_with_lambda`] that reports a singular system as
/// `None` instead of a zero step, so the Newton loop can tell "no step
/// exists" from "the step is zero".
pub(crate) fn try_solve_normal_equations(
    j: &[Vec<f64>],
    f: &[f64],
    n: usize,
    lambda: f64,
) -> Option<Vec<f64>> {
    let (mut a, mut b) = normal_equations(j, f, n);
    // Tikhonov regularisation
    for i in 0..n {
        a[i][i] += lambda;
    }

    gaussian_elimination(&mut a, &mut b)
}

/// The normal-equation pair `A = JᵀJ` (n×n) and `b = Jᵀ(−F)`.
///
/// A constraint row touches a handful of parameters, so `J` is mostly
/// zeros even in its dense storage. Only the non-zero entries of each row
/// are visited: the products that are skipped are exact zeros, so every sum
/// accumulates the same terms in the same order as the plain triple loop
/// and the result is bitwise identical — at a cost proportional to the
/// squared row density instead of `n²` per row.
pub(crate) fn normal_equations(j: &[Vec<f64>], f: &[f64], n: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    let mut a = vec![vec![0.0_f64; n]; n];
    let mut b = vec![0.0_f64; n];
    let mut nonzero: Vec<usize> = Vec::new();

    for (row, &fi) in j.iter().zip(f) {
        nonzero.clear();
        nonzero.extend((0..n).filter(|&col| row[col] != 0.0));
        for &col in &nonzero {
            b[col] -= row[col] * fi;
            for &k in &nonzero {
                a[col][k] += row[col] * row[k];
            }
        }
    }
    (a, b)
}

impl SparseJacobian {
    /// The normal-equation pair `A = JᵀJ` and `b = Jᵀ(−F)`: the sums of
    /// [`normal_equations`], term for term and in the same order, at a cost
    /// of the squared entries of each row and nothing per empty column.
    pub(crate) fn normal_equations(&self, f: &[f64]) -> (Square, Vec<f64>) {
        let n = self.n_cols;
        let mut a = Square::zeros(n);
        let mut b = vec![0.0_f64; n];
        let mut nonzero: Vec<(usize, f64)> = Vec::new();

        for (r, &fi) in f.iter().enumerate() {
            nonzero.clear();
            nonzero.extend(self.row(r));
            for &(col, v) in &nonzero {
                b[col] -= v * fi;
                let a_col = a.row_mut(col);
                for &(k, w) in &nonzero {
                    a_col[k] += v * w;
                }
            }
        }
        (a, b)
    }
}

/// Gaussian elimination with partial pivoting; `None` when a pivot falls
/// below [`PIVOT_TOL`] (singular system).
///
/// Non-finite input does not panic: a NaN entry is picked as the pivot (it
/// orders above every finite magnitude) and propagates into the result,
/// which the caller checks for finiteness.
pub(crate) fn gaussian_elimination(a: &mut [Vec<f64>], b: &mut [f64]) -> Option<Vec<f64>> {
    let n = b.len();
    debug_assert_eq!(a.len(), n);

    for col in 0..n {
        // Find pivot row (max absolute value in column)
        let pivot_row = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
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

/// LU factorisation with partial pivoting of a square matrix, kept so that
/// several right-hand sides can be solved against it.
///
/// Solving through the factorisation performs exactly the arithmetic of
/// [`gaussian_elimination`] on the right-hand side (same pivots, same
/// multipliers, same order), so the result is bitwise the one a fresh
/// elimination would give.
pub(crate) struct LuFactors {
    /// The factors in one matrix (rows in pivoted order): the upper
    /// triangle left by the elimination on and above the diagonal, and
    /// below it the elimination factors, each moved with its row by the
    /// later exchanges (as LAPACK stores them).
    lu: Square,
    /// The row exchanged with `col` at elimination step `col`.
    swaps: Vec<usize>,
}

/// A square matrix stored row after row in one buffer: one allocation
/// however many rows, where a vector per row cost one each.
#[derive(Clone)]
pub(crate) struct Square {
    n: usize,
    a: Vec<f64>,
}

impl Square {
    pub(crate) fn zeros(n: usize) -> Self {
        Self {
            n,
            a: vec![0.0; n * n],
        }
    }

    /// The matrix with the rows `rows` (each as long as there are rows).
    #[cfg(test)]
    pub(crate) fn from_rows(rows: &[Vec<f64>]) -> Self {
        let n = rows.len();
        let mut square = Self::zeros(n);
        for (target, row) in square.a.chunks_exact_mut(n.max(1)).zip(rows) {
            target.copy_from_slice(row);
        }
        square
    }

    #[cfg(test)]
    pub(crate) fn row(&self, i: usize) -> &[f64] {
        &self.a[i * self.n..(i + 1) * self.n]
    }

    pub(crate) fn row_mut(&mut self, i: usize) -> &mut [f64] {
        &mut self.a[i * self.n..(i + 1) * self.n]
    }

    /// Add `lambda` to every diagonal entry.
    pub(crate) fn add_diagonal(&mut self, lambda: f64) {
        for i in 0..self.n {
            self.a[i * self.n + i] += lambda;
        }
    }

    /// Exchange rows `i` and `j`.
    fn swap_rows(&mut self, i: usize, j: usize) {
        if i != j {
            let (lo, hi) = (i.min(j), i.max(j));
            let (first, second) = self.a.split_at_mut(hi * self.n);
            first[lo * self.n..(lo + 1) * self.n].swap_with_slice(&mut second[..self.n]);
        }
    }
}

impl LuFactors {
    /// Factor `a`; `None` when a pivot falls below [`PIVOT_TOL`].
    pub(crate) fn factor(mut a: Square) -> Option<Self> {
        let n = a.n;
        let mut swaps = Vec::with_capacity(n);

        for col in 0..n {
            let at = |i: usize| a.a[i * n + col].abs();
            let pivot_row = (col..n).max_by(|&i, &j| at(i).total_cmp(&at(j)))?;
            a.swap_rows(col, pivot_row);
            swaps.push(pivot_row);

            let pivot = a.a[col * n + col];
            if pivot.abs() < PIVOT_TOL {
                return None;
            }

            // The pivot row and the rows below as slices, so that the
            // inner loop runs without bounds checks: indexed, and with the
            // factors kept in a vector per column, the same arithmetic took
            // two to three times as long as `gaussian_elimination`. The
            // operations on the upper triangle and their order are those of
            // `gaussian_elimination`; the entry below the pivot, which it
            // drives to (about) zero and never reads again, keeps the
            // factor instead.
            let (above, below) = a.a.split_at_mut((col + 1) * n);
            let pivot_tail = &above[col * n + col + 1..];
            for row in below.chunks_exact_mut(n) {
                let factor = row[col] / pivot;
                for (x, p) in row[col + 1..].iter_mut().zip(pivot_tail) {
                    let v = p * factor;
                    *x -= v;
                }
                row[col] = factor;
            }
        }
        Some(Self { lu: a, swaps })
    }

    /// The smallest pivot over the largest, in magnitude: about the
    /// reciprocal of the condition number. NaN if a pivot is.
    pub(crate) fn pivot_ratio(&self) -> f64 {
        let (mut small, mut large) = (f64::INFINITY, 0.0_f64);
        for i in 0..self.lu.n {
            let p = self.lu.a[i * self.lu.n + i].abs();
            if p.is_nan() {
                return f64::NAN;
            }
            small = small.min(p);
            large = large.max(p);
        }
        small / large
    }

    /// Solve `A x = b` for the factored `A`. `b` is used as scratch space.
    pub(crate) fn solve(&self, b: &mut [f64]) -> Vec<f64> {
        let n = b.len();
        debug_assert_eq!(self.lu.n, n);
        let lu = &self.lu.a;

        // The factors moved with their rows, so the exchanges are applied
        // first; every entry then receives the updates the interleaved
        // elimination gives it, in the same order.
        for (col, &swap) in self.swaps.iter().enumerate() {
            b.swap(col, swap);
        }
        for col in 0..n {
            let pivot = b[col];
            let below = lu[(col + 1) * n..].chunks_exact(n);
            for (row, target) in below.zip(&mut b[col + 1..]) {
                let bv = pivot * row[col];
                *target -= bv;
            }
        }

        let mut x = vec![0.0_f64; n];
        for i in (0..n).rev() {
            let row = &lu[i * n..(i + 1) * n];
            x[i] = b[i];
            for j in (i + 1)..n {
                let v = row[j] * x[j];
                x[i] -= v;
            }
            x[i] /= row[i];
        }
        x
    }
}

/// The row-space Gram system of a set of independent rows `R`, factored once
/// and solved for many right-hand sides.
///
/// For a dependent row `d`, [`RowSpaceGram::coefficients`] returns the `y`
/// with `(R Rᵀ + λI) y = R d`: the combination of independent rows that
/// reproduces `d`. The diagnosis asks this for every dependent row; the
/// matrix is the same each time, so it is built and factored here once
/// instead of once per dependent row. Rows are sparse, so the products are
/// accumulated over non-zero entries only.
pub(crate) struct RowSpaceGram {
    /// For each column, the independent rows with a non-zero entry there.
    by_column: Vec<Vec<(usize, f64)>>,
    /// Number of independent rows.
    rank: usize,
    /// `None` when the regularised Gram matrix is numerically singular.
    factors: Option<LuFactors>,
}

impl RowSpaceGram {
    /// `rows` are the independent rows (each `n_cols` long).
    pub(crate) fn new(rows: &[Vec<f64>], n_cols: usize, lambda: f64) -> Self {
        let rank = rows.len();
        let mut by_column: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_cols];
        for (i, row) in rows.iter().enumerate() {
            for (col, &v) in row.iter().enumerate() {
                if v != 0.0 {
                    by_column[col].push((i, v));
                }
            }
        }

        let mut gram = Square::zeros(rank);
        for entries in &by_column {
            for &(i, vi) in entries {
                let gram_i = gram.row_mut(i);
                for &(j, vj) in entries {
                    gram_i[j] += vi * vj;
                }
            }
        }
        gram.add_diagonal(lambda);

        Self {
            by_column,
            rank,
            factors: LuFactors::factor(gram),
        }
    }

    /// Coefficients of `dependent` over the independent rows; all zeros when
    /// the Gram system is singular.
    pub(crate) fn coefficients(&self, dependent: &[f64]) -> Vec<f64> {
        let Some(factors) = &self.factors else {
            return vec![0.0_f64; self.rank];
        };
        let mut rhs = vec![0.0_f64; self.rank];
        for (entries, &d) in self.by_column.iter().zip(dependent) {
            if d != 0.0 {
                for &(i, vi) in entries {
                    rhs[i] += vi * d;
                }
            }
        }
        factors.solve(&mut rhs)
    }
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

    #[test]
    fn non_finite_residual_reports_diverged_and_restores() {
        // ln(x − 2) is NaN at the start x = 1.
        let mut params = vec![1.0_f64];
        let result = newton_solve(&mut params, &[false], 1, |p, f| {
            f[0] = (p[0] - 2.0).ln();
        });
        assert!(!result.converged);
        assert_eq!(result.status, SolveStatus::Diverged);
        assert_eq!(params[0], 1.0);
    }

    #[test]
    fn exact_stationary_point_without_a_root_reports_stalled() {
        // x² + 1 has no root, and at x = 0 its gradient is exactly zero:
        // every regularisation gives the zero step.
        let mut params = vec![0.0_f64];
        let result = newton_solve(&mut params, &[false], 1, |p, f| {
            f[0] = p[0] * p[0] + 1.0;
        });
        assert!(!result.converged);
        assert_eq!(result.status, SolveStatus::Stalled);
        assert_eq!(result.iterations, 0);
        assert_eq!(params[0], 0.0);
    }

    #[test]
    fn gaussian_elimination_survives_nan() {
        // 0.1.0 unwrapped a `partial_cmp` here and panicked.
        let mut a = vec![vec![f64::NAN, 1.0], vec![1.0, 2.0]];
        let mut b = vec![1.0, 2.0];
        if let Some(x) = gaussian_elimination(&mut a, &mut b) {
            assert!(x.iter().any(|v| !v.is_finite()), "NaN must propagate");
        }
    }

    /// A small dense, well-conditioned test matrix with a row exchange in
    /// its elimination (the largest first-column entry is not on top).
    fn sample_matrix() -> Vec<Vec<f64>> {
        vec![
            vec![0.5, 2.0, -1.0, 0.25],
            vec![4.0, 1.0, 0.5, -2.0],
            vec![-1.5, 0.75, 3.0, 1.0],
            vec![2.0, -0.5, 1.25, 5.0],
        ]
    }

    #[test]
    fn lu_factors_reproduce_gaussian_elimination_bitwise() {
        let factors =
            LuFactors::factor(Square::from_rows(&sample_matrix())).expect("regular matrix");
        for rhs in [
            vec![1.0, 2.0, 3.0, 4.0],
            vec![-0.3, 0.0, 7.5, 1.0e-3],
            vec![0.0, 0.0, 0.0, 0.0],
        ] {
            let mut a = sample_matrix();
            let mut b = rhs.clone();
            let direct = gaussian_elimination(&mut a, &mut b).expect("regular matrix");
            let through_factors = factors.solve(&mut rhs.clone());
            assert_eq!(
                direct.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                through_factors
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                "rhs {rhs:?}"
            );
        }

        // A larger, unsymmetric matrix, where many later exchanges move
        // factors stored by earlier steps.
        let mut seed = 1_u64;
        let mut random = || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
        };
        let n = 40;
        let a: Vec<Vec<f64>> = (0..n).map(|_| (0..n).map(|_| random()).collect()).collect();
        let rhs: Vec<f64> = (0..n).map(|_| random()).collect();
        let factors = LuFactors::factor(Square::from_rows(&a)).expect("regular matrix");
        let direct = gaussian_elimination(&mut a.clone(), &mut rhs.clone()).expect("regular");
        let through_factors = factors.solve(&mut rhs.clone());
        assert!(
            factors
                .swaps
                .iter()
                .enumerate()
                .filter(|(c, s)| *c != **s)
                .count()
                > n / 2
        );
        assert_eq!(
            direct.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            through_factors
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        );
    }

    /// Dense rows to the sparse rows `min_norm_step` takes.
    fn sparse(j: &[Vec<f64>], n: usize) -> SparseJacobian {
        SparseJacobian::from_dense(j, n)
    }

    fn close(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.iter()
            .zip(b)
            .all(|(x, y)| (x - y).abs() <= tol * (1.0 + y.abs()))
    }

    #[test]
    fn min_norm_step_solves_a_square_system() {
        let j = sample_matrix();
        let f = [1.0, -2.0, 0.5, 3.0];
        let h = min_norm_step(&sparse(&j, 4), &f).expect("regular");
        let mut b: Vec<f64> = f.iter().map(|v| -v).collect();
        let direct = gaussian_elimination(&mut j.clone(), &mut b).expect("regular");
        assert!(close(&h, &direct, 1e-13), "{h:?} vs {direct:?}");
    }

    #[test]
    fn min_norm_step_takes_the_shortest_of_many_solutions() {
        // Two equations, three unknowns, full row rank: h = −Jᵀ(JJᵀ)⁻¹f.
        let j = vec![vec![1.0, 2.0, -1.0], vec![0.5, -1.0, 3.0]];
        let f = [1.0, -2.0];
        let h = min_norm_step(&sparse(&j, 3), &f).expect("rank 2");
        let jjt = [
            [dot(&j[0], &j[0]), dot(&j[0], &j[1])],
            [dot(&j[1], &j[0]), dot(&j[1], &j[1])],
        ];
        let det = jjt[0][0] * jjt[1][1] - jjt[0][1] * jjt[1][0];
        let y = [
            (-f[0] * jjt[1][1] + f[1] * jjt[0][1]) / det,
            (-f[1] * jjt[0][0] + f[0] * jjt[1][0]) / det,
        ];
        let expect: Vec<f64> = (0..3).map(|c| j[0][c] * y[0] + j[1][c] * y[1]).collect();
        assert!(close(&h, &expect, 1e-13), "{h:?} vs {expect:?}");
    }

    #[test]
    fn min_norm_step_handles_rank_deficiency() {
        // Two identical columns: consistent, the step splits evenly.
        let h = min_norm_step(&sparse(&[vec![1.0, 1.0]], 2), &[-2.0]).unwrap();
        assert!(close(&h, &[1.0, 1.0], 1e-14), "{h:?}");
        // Two identical rows asking for different things: the least-squares
        // compromise (their mean), and of those the shortest.
        let j = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        let h = min_norm_step(&sparse(&j, 2), &[-1.0, -3.0]).unwrap();
        assert!(close(&h, &[1.0, 1.0], 1e-14), "{h:?}");
        // A dependent column among independent ones.
        let j = vec![
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 1.0],
            vec![1.0, 1.0, 2.0],
        ];
        let f = [-1.0, -2.0, -3.0];
        let h = min_norm_step(&sparse(&j, 3), &f).unwrap();
        // Solutions are h = (1, 2, 0) + t·(1, 1, −1); the shortest has t = −1.
        assert!(close(&h, &[0.0, 1.0, 1.0], 1e-13), "{h:?}");
    }

    #[test]
    fn min_norm_step_fits_an_overdetermined_system() {
        let j = vec![vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, 1.0]];
        let f = [-1.0, -1.0, -3.0];
        let h = min_norm_step(&sparse(&j, 2), &f).unwrap();
        // Normal equations: [[2, 1], [1, 2]]·h = [4, 4] → h = (4/3, 4/3).
        assert!(close(&h, &[4.0 / 3.0, 4.0 / 3.0], 1e-14), "{h:?}");
    }

    /// Two nearly parallel rows (`κ ≈ 1e7`): the step matches Cramer's rule
    /// to a few rounding units of `κ`, where the normal equations (`κ²`)
    /// would lose most of the digits.
    #[test]
    fn min_norm_step_keeps_the_digits_of_an_ill_conditioned_system() {
        let e = 1e-7;
        let j = vec![vec![1.0, e], vec![1.0, -e]];
        let f = [-1.0, -1.0 - 2e-9];
        let h = min_norm_step(&sparse(&j, 2), &f).unwrap();
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let expect = [
            (-f[0] * j[1][1] + f[1] * j[0][1]) / det,
            (-f[1] * j[0][0] + f[0] * j[1][0]) / det,
        ];
        assert!(close(&h, &expect, 1e-8), "{h:?} vs {expect:?}");
    }

    #[test]
    fn min_norm_step_reports_no_step_and_nan() {
        assert!(min_norm_step(&sparse(&[vec![0.0, 0.0]], 2), &[1.0]).is_none());
        let h = min_norm_step(&sparse(&[vec![f64::NAN, 1.0]], 2), &[1.0]).unwrap();
        assert!(h.iter().all(|v| v.is_nan()));
    }

    #[test]
    fn lu_factors_reject_a_singular_matrix() {
        let singular = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert!(LuFactors::factor(Square::from_rows(&singular)).is_none());
    }

    #[test]
    fn row_space_gram_matches_a_fresh_solve_per_row() {
        // Independent rows with zeros in them, and two rows to express.
        let rows = vec![
            vec![1.0, 0.0, 2.0, 0.0, -1.0],
            vec![0.0, 3.0, 0.0, 0.5, 0.0],
            vec![0.25, 0.0, 0.0, -2.0, 4.0],
        ];
        let n_cols = 5;
        let gram = RowSpaceGram::new(&rows, n_cols, LAMBDA);

        for dependent in [
            vec![1.25, 3.0, 2.0, -1.5, 3.0], // rows 0 + 1 + 2
            vec![0.0, -6.0, 0.0, -1.0, 0.0], // −2 × row 1
        ] {
            // The form the diagnosis used: rebuild and eliminate per row.
            let r = rows.len();
            let mut a = vec![vec![0.0_f64; r]; r];
            let mut b = vec![0.0_f64; r];
            for (i, ri) in rows.iter().enumerate() {
                b[i] = ri.iter().zip(&dependent).map(|(x, y)| x * y).sum::<f64>();
                for (j, rj) in rows.iter().enumerate() {
                    a[i][j] = ri.iter().zip(rj).map(|(x, y)| x * y).sum::<f64>();
                }
                a[i][i] += LAMBDA;
            }
            let reference = gaussian_elimination(&mut a, &mut b).expect("regular");

            let y = gram.coefficients(&dependent);
            for (got, want) in y.iter().zip(&reference) {
                assert!(
                    (got - want).abs() <= 1e-15 * want.abs().max(1.0),
                    "{got} vs {want}"
                );
            }
        }
    }

    #[test]
    fn normal_equations_match_the_dense_triple_loop_bitwise() {
        let j = vec![
            vec![1.5, 0.0, -2.0, 0.0],
            vec![0.0, 0.0, 3.0, 0.25],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![-1.0, 4.0, 0.0, 0.5],
        ];
        let f = vec![0.3, -1.2, 7.0, 2.5];
        let n = 4;

        let (a, b) = normal_equations(&j, &f, n);

        let mut a_ref = vec![vec![0.0_f64; n]; n];
        let mut b_ref = vec![0.0_f64; n];
        for i in 0..j.len() {
            for col in 0..n {
                b_ref[col] -= j[i][col] * f[i];
                for k in 0..n {
                    a_ref[col][k] += j[i][col] * j[i][k];
                }
            }
        }
        let bits =
            |m: &[Vec<f64>]| -> Vec<u64> { m.iter().flatten().map(|v| v.to_bits()).collect() };
        assert_eq!(bits(&a), bits(&a_ref));
        assert_eq!(
            b.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            b_ref.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );

        // The sparse Jacobian gives the same sums, whether its pattern holds
        // only the non-zero entries or every entry, zeros included.
        let mut every = SparseJacobian::dense(j.len(), n);
        every.copy_from_dense(&j);
        for sparse in [SparseJacobian::from_dense(&j, n), every] {
            let (a_sparse, b_sparse) = sparse.normal_equations(&f);
            let a_sparse: Vec<Vec<f64>> = (0..n).map(|i| a_sparse.row(i).to_vec()).collect();
            assert_eq!(bits(&a_sparse), bits(&a_ref));
            assert_eq!(
                b_sparse.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                b_ref.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
            );
        }
    }
}
