//! The state of every entity and constraint of a system, for a UI to show.
//!
//! [`Status`] gathers what [`crate::diagnosis::Diagnosis`] and the
//! per-parameter mobility say into one answer per entity and per
//! constraint. Each gets a single state, the most serious that applies,
//! so that a host can colour geometry by it directly; the facts behind it
//! stay available (`diagnose()`, `mobility()`).

use crate::diagnosis::Diagnosis;
use crate::numerics::{self, RESIDUAL_TOL};
use crate::system::{self, Equations};

/// The state of an entity, from the least to the most serious. An entity
/// gets the most serious one that applies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum EntityStatus {
    /// Nothing about it can move, and nothing on it is redundant.
    WellDefined,
    /// Some free parameter of it can still move without breaking any
    /// constraint (it is in the null space of the constraint Jacobian).
    UnderDefined,
    /// A constraint on it depends on others only where the geometry
    /// stands: at a configuration that is singular for it (an exact
    /// tangency), not in general (see [`ConstraintStatus::Singular`]). Its
    /// first-order motions there need not be real ones.
    Singular,
    /// A constraint on it is redundant (or partially redundant): removable
    /// without changing the solution.
    OverDefined,
    /// A constraint on it is not met where the geometry stands, in a
    /// component that is not at a solution (not solved yet, or a solve
    /// that failed).
    Unsolved,
    /// A constraint on it conflicts with the others.
    Conflicting,
}

/// The state of a constraint, from the least to the most serious. A
/// constraint gets the most serious one that applies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ConstraintStatus {
    /// Met, and independent of the others.
    Met,
    /// Met, and dependent on others only where the geometry stands: it
    /// could go without changing the rank here, but at a witness
    /// configuration — the geometry moved off, with the same geometric
    /// relations holding and every dimension at what it measures there —
    /// it is independent. The configuration is singular (an exact
    /// tangency, a dimension at a special value), not the constraint
    /// redundant. Only the 2D sketch looks for a witness; elsewhere such a
    /// constraint is [`ConstraintStatus::Redundant`].
    Singular,
    /// Dependent on others and consistent with them: removable without
    /// changing the solution. Every member of a consistent dependency group
    /// (see [`Diagnosis::dependency_groups`]) is: each alone could go, not
    /// all together.
    Redundant,
    /// Dependent on others and only nearly consistent with them (the gray
    /// band of [`Diagnosis::partially_redundant`]).
    PartiallyRedundant,
    /// Not met where the geometry stands, in a component that is not at a
    /// solution.
    Unmet,
    /// Inconsistent with the others: no geometry meets them all. Every
    /// member of a dependency group that holds a conflict is, not only the
    /// one [`Diagnosis::conflicting`] names.
    Conflicting,
    /// Suppressed or removed: not part of the system.
    Inactive,
}

/// One state per entity and per constraint (see [`EntityStatus`],
/// [`ConstraintStatus`]), with the degrees of freedom and the dependency
/// groups of the diagnosis it is built from.
#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Status {
    /// Degrees of freedom left: free parameters − rank (see
    /// [`Diagnosis::dof`]).
    pub dof: i64,
    /// Per entity id; `None` for an entity that was removed.
    pub entities: Vec<Option<EntityStatus>>,
    /// Per constraint index.
    pub constraints: Vec<ConstraintStatus>,
    /// The constraints that depend on each other, group by group (see
    /// [`Diagnosis::dependency_groups`]).
    pub dependency_groups: Vec<Vec<usize>>,
}

impl Status {
    /// The entities in state `state`, ascending.
    pub fn entities_in(&self, state: EntityStatus) -> Vec<usize> {
        (0..self.entities.len())
            .filter(|&e| self.entities[e] == Some(state))
            .collect()
    }

    /// The constraints in state `state`, ascending.
    pub fn constraints_in(&self, state: ConstraintStatus) -> Vec<usize> {
        (0..self.constraints.len())
            .filter(|&c| self.constraints[c] == state)
            .collect()
    }
}

/// What a family tells [`status_of`] about its entities.
pub(crate) struct Entities<'a> {
    /// Number of entities.
    pub(crate) count: usize,
    /// Of the constraints given (removable where the geometry stands, of a
    /// system of rank `.1` there), those that could still go at a witness
    /// configuration (see [`ConstraintStatus::Singular`]); `None` when the
    /// family has none or none was found.
    #[allow(clippy::type_complexity)]
    pub(crate) witness: Option<&'a dyn Fn(&[usize], usize) -> Option<Vec<usize>>>,
    /// Whether entity `e` was removed.
    pub(crate) removed: &'a dyn Fn(usize) -> bool,
    /// Whether some free parameter of entity `e` can move.
    pub(crate) movable: &'a dyn Fn(usize) -> bool,
    /// The entities constraint `c` refers to.
    pub(crate) of_constraint: &'a dyn Fn(usize) -> Vec<usize>,
}

/// The [`Status`] of the system `eqs` at `params` from its diagnosis.
pub(crate) fn status_of<E: Equations + ?Sized, G: system::Gauge>(
    eqs: &E,
    params: &[f64],
    fixed: &[bool],
    gauge: &G,
    diagnosis: Diagnosis,
    entities: Entities<'_>,
) -> Status {
    let n = eqs.constraint_count();
    let unmet = unmet_constraints(eqs, params, RESIDUAL_TOL);
    let mut constraints = vec![ConstraintStatus::Met; n];
    for (c, state) in constraints.iter_mut().enumerate() {
        if eqs.equation_count(c) == 0 {
            *state = ConstraintStatus::Inactive;
        } else if unmet[c] {
            *state = ConstraintStatus::Unmet;
        }
    }
    fn raise(constraints: &mut [ConstraintStatus], list: &[usize], to: ConstraintStatus) {
        for &c in list {
            if constraints[c] != ConstraintStatus::Inactive && constraints[c] < to {
                constraints[c] = to;
            }
        }
    }
    // A conflict involves every member of its dependency group, not only
    // the representative the pivoting named.
    raise(
        &mut constraints,
        &diagnosis.conflicting,
        ConstraintStatus::Conflicting,
    );
    let in_conflict = |group: &[usize]| group.iter().any(|c| diagnosis.conflicting.contains(c));
    for group in &diagnosis.dependency_groups {
        if in_conflict(group) {
            raise(&mut constraints, group, ConstraintStatus::Conflicting);
        }
    }
    // A redundant constraint is one that could go without changing the
    // solution: dependent rows are not enough (two point coincidences of a
    // rigid body share the distance along their line, yet neither can go),
    // so each candidate — a member of a consistent dependency, or a row the
    // diagnosis named — is taken out on its own and kept as redundant if
    // the rank of its component does not drop.
    let partial = |c: usize| {
        diagnosis.partially_redundant.contains(&c)
            || diagnosis.dependency_groups.iter().any(|g| {
                g.contains(&c) && g.iter().any(|m| diagnosis.partially_redundant.contains(m))
            })
    };
    let mut candidates: Vec<usize> = diagnosis
        .redundant
        .iter()
        .chain(&diagnosis.partially_redundant)
        .copied()
        .chain(
            diagnosis
                .dependency_groups
                .iter()
                .filter(|g| !in_conflict(g))
                .flat_map(|g| g.iter().copied()),
        )
        .filter(|&c| constraints[c] != ConstraintStatus::Conflicting)
        .collect();
    candidates.sort_unstable();
    candidates.dedup();
    let here = removable(eqs, params, fixed, gauge, &candidates);
    // Removable here but not at a witness: dependent only where the
    // geometry stands.
    let generic = match entities.witness {
        Some(witness) if !here.is_empty() => witness(&here, diagnosis.rank),
        _ => None,
    };
    for c in here {
        let to = if generic.as_ref().is_some_and(|g| !g.contains(&c)) {
            ConstraintStatus::Singular
        } else if partial(c) {
            ConstraintStatus::PartiallyRedundant
        } else {
            ConstraintStatus::Redundant
        };
        raise(&mut constraints, &[c], to);
    }

    let mut states: Vec<Option<EntityStatus>> = (0..entities.count)
        .map(|e| {
            if (entities.removed)(e) {
                None
            } else if (entities.movable)(e) {
                Some(EntityStatus::UnderDefined)
            } else {
                Some(EntityStatus::WellDefined)
            }
        })
        .collect();
    for (c, &state) in constraints.iter().enumerate() {
        let to = match state {
            ConstraintStatus::Met | ConstraintStatus::Inactive => continue,
            ConstraintStatus::Singular => EntityStatus::Singular,
            ConstraintStatus::Redundant | ConstraintStatus::PartiallyRedundant => {
                EntityStatus::OverDefined
            }
            ConstraintStatus::Unmet => EntityStatus::Unsolved,
            ConstraintStatus::Conflicting => EntityStatus::Conflicting,
        };
        for e in (entities.of_constraint)(c) {
            if let Some(Some(current)) = states.get_mut(e) {
                if *current < to {
                    *current = to;
                }
            }
        }
    }
    Status {
        dof: diagnosis.dof,
        entities: states,
        constraints,
        dependency_groups: diagnosis.dependency_groups,
    }
}

/// The constraints of `eqs` that `active` marks; the others have no
/// equations.
struct Only<'a, E: ?Sized> {
    eqs: &'a E,
    active: Vec<bool>,
}

impl<E: Equations + ?Sized> Equations for Only<'_, E> {
    fn constraint_count(&self) -> usize {
        self.eqs.constraint_count()
    }
    fn equation_count(&self, index: usize) -> usize {
        if self.active[index] {
            self.eqs.equation_count(index)
        } else {
            0
        }
    }
    fn residuals(&self, index: usize, params: &[f64], out: &mut [f64]) {
        if self.active[index] {
            self.eqs.residuals(index, params, out);
        }
    }
    fn jacobian(&self, index: usize, params: &[f64], out: &mut [system::SparseRow]) {
        if self.active[index] {
            self.eqs.jacobian(index, params, out);
        }
    }
    fn column_count(&self, params: &[f64]) -> usize {
        self.eqs.column_count(params)
    }
    fn step_space(&self, free: &[usize]) -> Box<dyn crate::numerics::StepSpace + '_> {
        self.eqs.step_space(free)
    }
    fn steps_are_additive(&self) -> bool {
        self.eqs.steps_are_additive()
    }
    fn column_magnitude(&self, column: usize, params: &[f64]) -> f64 {
        self.eqs.column_magnitude(column, params)
    }
    fn length_scale(&self, group: &[usize], params: &[f64]) -> Option<f64> {
        self.eqs.length_scale(group, params)
    }
    fn debug_check_jacobian(&self) -> bool {
        self.eqs.debug_check_jacobian()
    }
}

/// The constraints among `candidates` that can each be taken out alone
/// without lowering the rank of their component: the rank of the
/// component's Jacobian with and without it, from the same diagnosis (and
/// gauge) the status is built from.
pub(crate) fn removable<E: Equations + ?Sized, G: system::Gauge>(
    eqs: &E,
    params: &[f64],
    fixed: &[bool],
    gauge: &G,
    candidates: &[usize],
) -> Vec<usize> {
    if candidates.is_empty() {
        return Vec::new();
    }
    let n = eqs.constraint_count();
    let rank_of = |active: Vec<bool>, fixed: &[bool]| {
        system::diagnose(&Only { eqs, active }, params, fixed, gauge).rank
    };
    let mut out = Vec::new();
    for group in &system::components(eqs, params) {
        let mine: Vec<usize> = candidates
            .iter()
            .copied()
            .filter(|c| group.contains(c))
            .collect();
        if mine.is_empty() {
            continue;
        }
        // An anchor is a constraint of full rank on its own free columns
        // (a fixed point, a grounded body): it decides them alone, so with
        // them held the rank of the rest is the rank of the component less
        // their number — and a constraint other than an anchor can be
        // tested on the smaller pieces the rest falls into (a hub with
        // forty parts pinned to it falls into forty). Only without a
        // gauge: under one, full rank in the stored parameters is not the
        // same test.
        let anchors: Vec<usize> = if gauge.is_trivial() {
            group
                .iter()
                .copied()
                .filter(|&c| full_rank_on_own_columns(eqs, c, params, fixed))
                .collect()
        } else {
            Vec::new()
        };
        let mut held = fixed.to_vec();
        for &a in &anchors {
            for column in columns_of(eqs, a, params) {
                held[column] = true;
            }
        }
        // The pieces of the component with the anchors' columns held: the
        // other constraints joined by the columns they share and nobody
        // holds.
        let rest: Vec<usize> = group
            .iter()
            .copied()
            .filter(|c| !anchors.contains(c))
            .collect();
        let pieces = pieces_of(eqs, &rest, params, &held);
        let mut whole_rank = None;
        for c in mine {
            if anchors.contains(&c) || anchors.is_empty() {
                // Tested on the whole component.
                let mut active = vec![false; n];
                for &k in group {
                    active[k] = true;
                }
                let full = *whole_rank.get_or_insert_with(|| rank_of(active.clone(), fixed));
                active[c] = false;
                if rank_of(active, fixed) == full {
                    out.push(c);
                }
                continue;
            }
            let piece = pieces.iter().find(|p| p.contains(&c)).expect("in a piece");
            let mut active = vec![false; n];
            for &k in piece {
                active[k] = true;
            }
            let full = rank_of(active.clone(), &held);
            active[c] = false;
            if rank_of(active, &held) == full {
                out.push(c);
            }
        }
    }
    out
}

/// The free columns constraint `c` lists, ascending.
fn columns_of<E: Equations + ?Sized>(eqs: &E, c: usize, params: &[f64]) -> Vec<usize> {
    let mut rows: Vec<system::SparseRow> = vec![Vec::new(); eqs.equation_count(c)];
    eqs.jacobian(c, params, &mut rows);
    let mut columns: Vec<usize> = rows.iter().flatten().map(|&(column, _)| column).collect();
    columns.sort_unstable();
    columns.dedup();
    columns
}

/// Whether the equations of constraint `c` alone have full rank on the
/// free columns it lists (at least one).
fn full_rank_on_own_columns<E: Equations + ?Sized>(
    eqs: &E,
    c: usize,
    params: &[f64],
    fixed: &[bool],
) -> bool {
    let mut rows: Vec<system::SparseRow> = vec![Vec::new(); eqs.equation_count(c)];
    eqs.jacobian(c, params, &mut rows);
    let columns: Vec<usize> = columns_of(eqs, c, params)
        .into_iter()
        .filter(|&column| !fixed[column])
        .collect();
    if columns.is_empty() || rows.len() < columns.len() {
        return false;
    }
    // Aᵀ with one column per equation (row-equilibrated), as the diagnosis
    // takes its rank.
    let mut a = vec![vec![0.0_f64; rows.len()]; columns.len()];
    for (r, row) in rows.iter().enumerate() {
        let scale = row.iter().fold(0.0_f64, |m, &(_, v)| m.max(v.abs()));
        for &(column, v) in row {
            if let Ok(k) = columns.binary_search(&column) {
                a[k][r] = if scale > 0.0 { v / scale } else { 0.0 };
            }
        }
    }
    let (_, rank) =
        numerics::householder_qr_column_pivoting(&mut a, crate::diagnosis::RANK_REL_TOL);
    rank == columns.len()
}

/// The constraints `indices` grouped by the columns they share that
/// `held` does not hold (union–find over constraints).
fn pieces_of<E: Equations + ?Sized>(
    eqs: &E,
    indices: &[usize],
    params: &[f64],
    held: &[bool],
) -> Vec<Vec<usize>> {
    use std::collections::HashMap;
    let mut parent: Vec<usize> = (0..indices.len()).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let mut first_of_column: HashMap<usize, usize> = HashMap::new();
    for (k, &c) in indices.iter().enumerate() {
        for column in columns_of(eqs, c, params) {
            if held[column] {
                continue;
            }
            match first_of_column.get(&column) {
                Some(&other) => {
                    let (a, b) = (find(&mut parent, k), find(&mut parent, other));
                    if a != b {
                        parent[a] = b;
                    }
                }
                None => {
                    first_of_column.insert(column, k);
                }
            }
        }
    }
    let mut by_root: HashMap<usize, Vec<usize>> = HashMap::new();
    for (k, &c) in indices.iter().enumerate() {
        let root = find(&mut parent, k);
        by_root.entry(root).or_default().push(c);
    }
    by_root.into_values().collect()
}

/// The motions the constraints leave one entity, each as an orthonormal
/// basis over the entity's own columns (its parameters, or a rigid body's
/// six coordinates), scaled as [`motion_spaces`] says.
#[derive(Debug, Clone, Default)]
pub(crate) struct MotionSpace {
    /// The motions of the entity while the rest may move along: the null
    /// space of the constraint Jacobian restricted to the entity.
    pub(crate) follow: Vec<Vec<f64>>,
    /// The motions of the entity with every other column held: the null
    /// space of the Jacobian's columns of the entity alone. Left empty
    /// unless asked for.
    pub(crate) own: Vec<Vec<f64>>,
}

/// Adds `v` to the orthonormal set `basis` (over `dim` numbers) when it is
/// not already in its span, by more than `tol`.
fn extend_basis(basis: &mut Vec<Vec<f64>>, dim: usize, mut v: Vec<f64>, tol: f64) {
    if basis.len() >= dim {
        return;
    }
    for _ in 0..2 {
        for u in basis.iter() {
            let d: f64 = v.iter().zip(u).map(|(a, b)| a * b).sum();
            for (a, b) in v.iter_mut().zip(u) {
                *a -= d * b;
            }
        }
    }
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > tol {
        for x in &mut v {
            *x /= norm;
        }
        basis.push(v);
    }
}

/// Per entity (`count` of them, entity `e` stepped by the columns
/// `params_of(e)`: its parameters, or a rigid body's six coordinates), the
/// motions of it that no constraint resists and that are physical, not
/// just another way of storing it ([`MotionSpace`]). A motion is a
/// direction of the null space of the constraint Jacobian — over the
/// effective columns of `gauge`, which leave out what a representation
/// stores twice (a plane's normal length) — mapped back to the stored
/// parameters, scaled so that its largest entry is one, restricted to the
/// entity and with the entity's own soft gauge directions taken out (a
/// cylinder's axis point sliding along the axis); it counts when something
/// of it is still above the rank tolerance. A free parameter that no
/// constraint touches moves on its own. `own` also gives each entity's
/// motions with all other columns held.
pub(crate) fn motion_spaces<E: Equations + ?Sized, G: system::Gauge>(
    eqs: &E,
    params: &[f64],
    fixed: &[bool],
    gauge: &G,
    count: usize,
    params_of: &dyn Fn(usize) -> std::ops::Range<usize>,
    own: bool,
) -> Vec<MotionSpace> {
    use std::collections::HashMap;
    let tol = crate::diagnosis::RANK_REL_TOL;
    // Over the columns of the family (a rigid body stores seven numbers
    // and moves in six).
    let n_cols = eqs.column_count(params);
    let mut owner = vec![usize::MAX; n_cols];
    for e in 0..count {
        for p in params_of(e) {
            owner[p] = e;
        }
    }
    let soft = gauge.soft_directions(params, fixed);
    let mut spaces = vec![MotionSpace::default(); count];
    // The part of `motion` (stored parameter → value, largest entry one)
    // on entity `e`, its own soft gauge directions (unit vectors with
    // disjoint supports) taken out; `None` when nothing of it is left.
    let local = |e: usize, motion: &HashMap<usize, f64>| -> Option<Vec<f64>> {
        let range = params_of(e);
        let mut v: Vec<f64> = range
            .clone()
            .map(|p| motion.get(&p).copied().unwrap_or(0.0))
            .collect();
        for g in soft
            .iter()
            .filter(|g| g.iter().all(|&(p, _)| owner[p] == e))
        {
            let along: f64 = g.iter().map(|&(p, w)| v[p - range.start] * w).sum();
            for &(p, w) in g {
                v[p - range.start] -= along * w;
            }
        }
        v.iter().any(|x| x.abs() > tol).then_some(v)
    };
    let mut touched_any = vec![false; n_cols];
    for group in system::components(eqs, params) {
        let mut touched = system::touched_params(eqs, &group, params);
        for &p in &touched {
            touched_any[p] = true;
        }
        touched.retain(|&p| !fixed[p]);
        if touched.is_empty() {
            continue;
        }
        let n_eq = system::equation_total(eqs, &group);
        let col_of: HashMap<usize, usize> =
            touched.iter().enumerate().map(|(c, &p)| (p, c)).collect();
        let mut j_raw = vec![vec![0.0_f64; touched.len()]; n_eq];
        system::analytic_jacobian(eqs, &group, params, &col_of, &mut j_raw);
        // Effective columns as combinations of the local ones (the local
        // ones themselves without a gauge).
        let columns: Vec<Vec<(usize, f64)>> = gauge
            .effective_columns(params, fixed, &touched)
            .unwrap_or_else(|| (0..touched.len()).map(|c| vec![(c, 1.0)]).collect());
        let j: Vec<Vec<f64>> = j_raw
            .iter()
            .map(|row| {
                columns
                    .iter()
                    .map(|col| col.iter().map(|&(lc, w)| row[lc] * w).sum())
                    .collect()
            })
            .collect();
        let null_space = |j: &[Vec<f64>], width: usize| -> Vec<Vec<f64>> {
            if n_eq == 0 {
                (0..width)
                    .map(|c| {
                        let mut v = vec![0.0; width];
                        v[c] = 1.0;
                        v
                    })
                    .collect()
            } else {
                numerics::null_space_basis(j, tol)
            }
        };
        // A null vector over the effective columns `cols` as a motion of
        // the stored parameters, largest entry one.
        let to_motion = |b: &[f64], cols: &[usize]| -> HashMap<usize, f64> {
            let mut motion: HashMap<usize, f64> = HashMap::new();
            for (k, &coeff) in b.iter().enumerate() {
                if coeff != 0.0 {
                    for &(lc, w) in &columns[cols[k]] {
                        *motion.entry(touched[lc]).or_insert(0.0) += coeff * w;
                    }
                }
            }
            let scale = motion.values().fold(0.0_f64, |m, v| m.max(v.abs()));
            if scale > 0.0 {
                for v in motion.values_mut() {
                    *v /= scale;
                }
            }
            motion
        };
        let mut present: Vec<usize> = touched
            .iter()
            .map(|&p| owner[p])
            .filter(|&e| e != usize::MAX)
            .collect();
        present.sort_unstable();
        present.dedup();
        let all: Vec<usize> = (0..columns.len()).collect();
        for b in null_space(&j, columns.len()) {
            let motion = to_motion(&b, &all);
            for &e in &present {
                if let Some(v) = local(e, &motion) {
                    let dim = params_of(e).len();
                    extend_basis(&mut spaces[e].follow, dim, v, tol);
                }
            }
        }
        if !own {
            continue;
        }
        for &e in &present {
            let mine: Vec<usize> = (0..columns.len())
                .filter(|&c| columns[c].iter().all(|&(lc, _)| owner[touched[lc]] == e))
                .collect();
            if mine.is_empty() {
                continue;
            }
            let j_mine: Vec<Vec<f64>> = j
                .iter()
                .map(|row| mine.iter().map(|&c| row[c]).collect())
                .collect();
            let dim = params_of(e).len();
            for b in null_space(&j_mine, mine.len()) {
                if let Some(v) = local(e, &to_motion(&b, &mine)) {
                    extend_basis(&mut spaces[e].own, dim, v, tol);
                }
            }
        }
    }
    // A free parameter no constraint touches can move as it likes.
    for (p, &e) in owner.iter().enumerate() {
        if e != usize::MAX && !fixed[p] && !touched_any[p] {
            let range = params_of(e);
            let mut v = vec![0.0; range.len()];
            v[p - range.start] = 1.0;
            extend_basis(&mut spaces[e].follow, range.len(), v.clone(), tol);
            if own {
                extend_basis(&mut spaces[e].own, range.len(), v, tol);
            }
        }
    }
    spaces
}

/// Per entity (`count` of them, entity `e` stepped by the columns
/// `params_of(e)`), whether some motion of it that no constraint resists
/// is physical (see [`motion_spaces`]). Without a gauge this is "some
/// parameter of it is movable" (see [`crate::mobility`]).
pub(crate) fn movable_entities<E: Equations + ?Sized, G: system::Gauge>(
    eqs: &E,
    params: &[f64],
    fixed: &[bool],
    gauge: &G,
    count: usize,
    params_of: &dyn Fn(usize) -> std::ops::Range<usize>,
) -> Vec<bool> {
    motion_spaces(eqs, params, fixed, gauge, count, params_of, false)
        .iter()
        .map(|m| !m.follow.is_empty())
        .collect()
}

/// Per constraint, whether it is not met where the geometry stands: some
/// equation of it is off by more than `tolerance` in its own units (the
/// size of a small component, the precision of its numbers — as a solve
/// judges it), in a component that, judged as a solve judges it, is not at
/// a solution. A component at a solution has none, whatever its rounding.
pub(crate) fn unmet_constraints<E: Equations + ?Sized>(
    eqs: &E,
    params: &[f64],
    tolerance: f64,
) -> Vec<bool> {
    let mut unmet = vec![false; eqs.constraint_count()];
    for group in system::components(eqs, params) {
        let n_eq = system::equation_total(eqs, &group);
        if n_eq == 0 {
            continue;
        }
        let unit = system::equation_unit(eqs.length_scale(&group, params));
        let f = system::eval_residuals(eqs, &group, params);
        let mut size = vec![0.0_f64; n_eq];
        system::row_magnitudes(eqs, &group, params, &mut size);
        let weighed: Vec<f64> = f
            .iter()
            .zip(&size)
            .map(|(r, s)| r / numerics::weight_of(*s, tolerance, unit))
            .collect();
        let rms = numerics::rms(&weighed);
        if rms < tolerance {
            continue;
        }
        let mut row = 0;
        for &c in &group {
            let count = eqs.equation_count(c);
            unmet[c] = weighed[row..row + count]
                .iter()
                .any(|v| v.is_nan() || v.abs() > tolerance);
            row += count;
        }
    }
    unmet
}
