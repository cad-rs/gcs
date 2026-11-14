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
    /// Whether entity `e` was removed.
    pub(crate) removed: &'a dyn Fn(usize) -> bool,
    /// Whether some free parameter of entity `e` can move.
    pub(crate) movable: &'a dyn Fn(usize) -> bool,
    /// The entities constraint `c` refers to.
    pub(crate) of_constraint: &'a dyn Fn(usize) -> Vec<usize>,
}

/// The [`Status`] of the system `eqs` at `params` from its diagnosis.
pub(crate) fn status_of<E: Equations + ?Sized>(
    eqs: &E,
    params: &[f64],
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
    let mut raise = |list: &[usize], to: ConstraintStatus| {
        for &c in list {
            if constraints[c] != ConstraintStatus::Inactive && constraints[c] < to {
                constraints[c] = to;
            }
        }
    };
    raise(&diagnosis.redundant, ConstraintStatus::Redundant);
    raise(
        &diagnosis.partially_redundant,
        ConstraintStatus::PartiallyRedundant,
    );
    raise(&diagnosis.conflicting, ConstraintStatus::Conflicting);
    // The lists above name one representative per dependency, whichever
    // the pivoting met last; a dependency group is the whole of one. Every
    // member of a group takes the state of the group: conflicting if it
    // holds a conflicting constraint, otherwise (partially) redundant —
    // each member alone could go without changing the solution.
    for group in &diagnosis.dependency_groups {
        let has = |list: &[usize]| group.iter().any(|c| list.contains(c));
        let state = if has(&diagnosis.conflicting) {
            ConstraintStatus::Conflicting
        } else if has(&diagnosis.partially_redundant) {
            ConstraintStatus::PartiallyRedundant
        } else {
            ConstraintStatus::Redundant
        };
        raise(group, state);
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

/// Per entity (`count` of them, entity `e` stored in the parameters
/// `params_of(e)`), whether some motion of it that no constraint resists
/// is physical, not just another way of storing it: a direction of the
/// null space of the constraint Jacobian — over the effective columns of
/// `gauge`, which leave out what a representation stores twice (a plane's
/// normal length) — that, restricted to the entity's parameters and with
/// the entity's own soft gauge directions taken out (a cylinder's axis
/// point sliding along the axis), still moves something. A free parameter
/// that no constraint touches makes its entity movable outright. Without a
/// gauge this is "some parameter of it is movable" (see
/// [`crate::mobility`]).
pub(crate) fn movable_entities<E: Equations + ?Sized, G: system::Gauge>(
    eqs: &E,
    params: &[f64],
    fixed: &[bool],
    gauge: &G,
    count: usize,
    params_of: &dyn Fn(usize) -> std::ops::Range<usize>,
) -> Vec<bool> {
    use std::collections::HashMap;
    let tol = crate::diagnosis::RANK_REL_TOL;
    let mut owner = vec![usize::MAX; params.len()];
    for e in 0..count {
        for p in params_of(e) {
            owner[p] = e;
        }
    }
    let soft = gauge.soft_directions(params, fixed);
    let mut movable = vec![false; count];
    let mut touched_any = vec![false; params.len()];
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
        let basis = if n_eq == 0 {
            (0..columns.len())
                .map(|c| {
                    let mut v = vec![0.0; columns.len()];
                    v[c] = 1.0;
                    v
                })
                .collect()
        } else {
            numerics::null_space_basis(&j, tol)
        };
        for b in basis {
            // The motion in the stored parameters.
            let mut motion: HashMap<usize, f64> = HashMap::new();
            for (c, &coeff) in b.iter().enumerate() {
                if coeff != 0.0 {
                    for &(lc, w) in &columns[c] {
                        *motion.entry(touched[lc]).or_insert(0.0) += coeff * w;
                    }
                }
            }
            let scale = motion.values().fold(0.0_f64, |m, v| m.max(v.abs()));
            if scale == 0.0 {
                continue;
            }
            let mut by_entity: HashMap<usize, Vec<(usize, f64)>> = HashMap::new();
            for (&p, &v) in &motion {
                if owner[p] != usize::MAX {
                    by_entity.entry(owner[p]).or_default().push((p, v / scale));
                }
            }
            for (e, mut part) in by_entity {
                if movable[e] {
                    continue;
                }
                // Take out the entity's own soft gauge directions (unit
                // vectors with disjoint supports).
                for g in soft
                    .iter()
                    .filter(|g| g.iter().all(|&(p, _)| owner[p] == e))
                {
                    let along: f64 = g
                        .iter()
                        .map(|&(p, w)| {
                            part.iter()
                                .find(|(q, _)| *q == p)
                                .map_or(0.0, |(_, v)| v * w)
                        })
                        .sum();
                    for &(p, w) in g {
                        if let Some(entry) = part.iter_mut().find(|(q, _)| *q == p) {
                            entry.1 -= along * w;
                        }
                    }
                }
                if part.iter().any(|(_, v)| v.abs() > tol) {
                    movable[e] = true;
                }
            }
        }
    }
    // A free parameter no constraint touches can move as it likes.
    for (p, &e) in owner.iter().enumerate() {
        if e != usize::MAX && !fixed[p] && !touched_any[p] {
            movable[e] = true;
        }
    }
    movable
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
