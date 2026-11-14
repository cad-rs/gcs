//! The main [`Sketch`] API.

use std::sync::Arc;

use glam::DVec2;

use crate::constraint::Constraint;
use crate::edit::{ConstraintState, SketchError};
use crate::entity::{Entity, EntityId, EntityKind, PointRef};
use crate::numerics::CoreOptions;
use crate::solver::{SketchEquations, SolveOptions, SolveResult};
use crate::system::{self, NoGauge};

/// Detailed DOF analysis result.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub struct DofReport {
    /// Number of free (non-fixed) parameters.
    pub free_params: usize,
    /// Total number of constraint equations.
    pub equations: usize,
    /// Rank of the constraint Jacobian over the free parameters at the
    /// current configuration.
    #[cfg_attr(feature = "serde", serde(default))]
    pub rank: usize,
    /// True degrees of freedom = `free_params − rank`. Redundant equations
    /// consume none, so this is never negative.
    pub dof: i64,
    /// Count-based estimate = `free_params − equations`. Negative when there
    /// are more equations than free parameters.
    #[cfg_attr(feature = "serde", serde(default))]
    pub structural_dof: i64,
    /// Indices of free parameters in the sketch's parameter vector.
    pub free_param_indices: Vec<usize>,
}

/// A 2D parametric sketch.
///
/// Holds geometric entities (points, lines, circles, arcs) and constraints
/// between them.  Call [`solve`][Sketch::solve] to run the GCS solver.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Sketch {
    /// Flat parameter vector shared by all entities.
    pub(crate) params: Vec<f64>,
    /// `fixed[i] == true` means `params[i]` is held constant during solving.
    pub(crate) fixed: Vec<bool>,
    /// Entity metadata.
    pub entities: Vec<Entity>,
    /// Constraint list. Slots are never reused: see `states`.
    pub(crate) constraints: Vec<Constraint>,
    /// The state of each constraint slot (a missing entry means active).
    #[cfg_attr(feature = "serde", serde(default))]
    states: Vec<ConstraintState>,
    /// Per entity: removed (a missing entry means present).
    #[cfg_attr(feature = "serde", serde(default))]
    removed: Vec<bool>,
    /// See [`Sketch::topology_revision`].
    #[cfg_attr(feature = "serde", serde(default))]
    revision: u64,
    /// What the last solve kept of the structure (see [`PlanCache`]).
    #[cfg_attr(feature = "serde", serde(skip))]
    plans: PlanCache,
}

/// The components of a sketch and what each one's solve needs that depends
/// on the structure alone (`system::Plans`), kept from one solve to the
/// next while the structure holds: while the topology revision is the same
/// and so are the entities (a host may edit the public entity list
/// directly, which the revision does not see). A clone shares them.
#[derive(Clone, Default)]
struct PlanCache {
    key: Option<(u64, u64)>,
    plans: Option<Arc<system::Plans>>,
}

impl std::fmt::Debug for PlanCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PlanCache")
            .field("key", &self.key)
            .field("kept", &self.plans.is_some())
            .finish()
    }
}

impl PlanCache {
    /// The plans for `key`, made by `make` unless they are kept for it.
    fn get(&mut self, key: (u64, u64), make: impl FnOnce() -> system::Plans) -> Arc<system::Plans> {
        if self.key != Some(key) || self.plans.is_none() {
            self.plans = Some(Arc::new(make()));
            self.key = Some(key);
        }
        Arc::clone(self.plans.as_ref().expect("just made"))
    }
}

impl Sketch {
    /// Create an empty sketch.
    pub fn new() -> Self {
        Self {
            params: Vec::new(),
            fixed: Vec::new(),
            entities: Vec::new(),
            constraints: Vec::new(),
            states: Vec::new(),
            removed: Vec::new(),
            revision: 0,
            plans: PlanCache::default(),
        }
    }

    /// The key of [`PlanCache`]: the topology revision, and a digest of
    /// what the structure is read from besides (the entity list, the sizes
    /// of the parameter, constraint and state lists).
    fn plan_key(&self) -> (u64, u64) {
        // FNV-1a over 64-bit words.
        let mut digest: u64 = 0xcbf2_9ce4_8422_2325;
        let mut mix = |word: u64| {
            digest ^= word;
            digest = digest.wrapping_mul(0x0000_0100_0000_01b3);
        };
        mix(self.params.len() as u64);
        mix(self.constraints.len() as u64);
        mix(self.states.len() as u64);
        for e in &self.entities {
            mix(e.kind as u64);
            mix(e.param_start as u64);
        }
        (self.revision, digest)
    }

    // ── Entity constructors ───────────────────────────────────────────────────

    /// Add a 2D point.  Returns its [`EntityId`].
    pub fn add_point(&mut self, x: f64, y: f64) -> EntityId {
        self.push_entity(EntityKind::Point, &[x, y])
    }

    /// Add a line defined by two endpoints (x1,y1)→(x2,y2).
    pub fn add_line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64) -> EntityId {
        self.push_entity(EntityKind::Line, &[x1, y1, x2, y2])
    }

    /// Add a circle with center (cx,cy) and radius r.
    pub fn add_circle(&mut self, cx: f64, cy: f64, r: f64) -> EntityId {
        self.push_entity(EntityKind::Circle, &[cx, cy, r])
    }

    /// Add an arc with center (cx,cy), radius r, from `start_angle` to
    /// `end_angle` (both in radians, counter-clockwise).
    pub fn add_arc(
        &mut self,
        cx: f64,
        cy: f64,
        r: f64,
        start_angle: f64,
        end_angle: f64,
    ) -> EntityId {
        self.push_entity(EntityKind::Arc, &[cx, cy, r, start_angle, end_angle])
    }

    fn push_entity(&mut self, kind: EntityKind, init: &[f64]) -> EntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        self.params.extend_from_slice(init);
        self.fixed.resize(self.fixed.len() + init.len(), false);
        self.entities.push(Entity::new(kind, param_start));
        self.revision += 1;
        id
    }

    // ── Parameter access ──────────────────────────────────────────────────────

    /// Read the current (x, y) of a point-like entity or point reference.
    pub fn point_coords(&self, p: PointRef) -> DVec2 {
        let (xi, yi) = p.param_indices(&self.entities);
        DVec2::new(self.params[xi], self.params[yi])
    }

    /// Read all parameters of an entity as a slice.
    pub fn entity_params(&self, id: EntityId) -> &[f64] {
        let e = &self.entities[id];
        &self.params[e.param_start..e.param_start + e.kind.param_count()]
    }

    /// Parameter index range owned by an entity.
    pub fn entity_param_range(&self, id: EntityId) -> std::ops::Range<usize> {
        let e = &self.entities[id];
        e.param_start..e.param_start + e.kind.param_count()
    }

    /// Per-parameter mobility at the current configuration (see
    /// [`crate::mobility::mobility`]). Suppressed and removed constraints
    /// hold nothing.
    pub fn mobility(&self) -> Vec<bool> {
        crate::mobility::mobility_of(&self.equations(), &self.params, &self.fixed)
    }

    /// Copy of the full flat parameter vector (persistence-contract and
    /// inspection use).
    pub fn params_snapshot(&self) -> Vec<f64> {
        self.params.clone()
    }

    /// Overwrite the full flat parameter vector (corpus drift probes and
    /// warm-start experiments). Length must match.
    pub fn set_params(&mut self, params: Vec<f64>) {
        assert_eq!(
            params.len(),
            self.params.len(),
            "param vector length mismatch"
        );
        self.params = params;
    }

    /// Apply a deterministic seeded random drift of amplitude `amp` to the
    /// FREE parameters only — masked (fixed) parameters are untouched, so
    /// the drift probe never silently re-positions fixed geometry.
    pub fn drift_free_params(&mut self, amp: f64, seed: u64) {
        let mut state = seed ^ 0x9E3779B97F4A7C15;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let x = state ^ (state >> 29);
            let x = x.wrapping_mul(0xBF58476D1CE4E5B9);
            (x >> 11) as f64 / (1u64 << 53) as f64
        };
        for (i, v) in self.params.iter_mut().enumerate() {
            if !self.fixed[i] {
                *v += (next() * 2.0 - 1.0) * amp;
            }
        }
    }

    /// Fix a single parameter (by absolute index) so the solver won't move it.
    pub fn fix_param(&mut self, param_idx: usize) {
        self.fixed[param_idx] = true;
        self.revision += 1;
    }

    /// Fix all parameters of an entity (equivalent to a rigid body constraint).
    pub fn fix_entity(&mut self, id: EntityId) {
        let e = &self.entities[id];
        for i in e.param_start..e.param_start + e.kind.param_count() {
            self.fixed[i] = true;
        }
        self.revision += 1;
    }

    // ── Constraint management ─────────────────────────────────────────────────

    /// Add a constraint. Returns its index, which stays valid for the life
    /// of the sketch: removing other constraints does not shift it.
    pub fn add_constraint(&mut self, c: Constraint) -> usize {
        self.constraints.push(c);
        self.states
            .resize(self.constraints.len(), ConstraintState::Active);
        self.revision += 1;
        self.constraints.len() - 1
    }

    /// Every constraint slot, in index order — including suppressed and
    /// removed ones, whose state [`Sketch::constraint_state`] tells.
    pub fn constraints(&self) -> &[Constraint] {
        &self.constraints
    }

    /// Number of constraint slots (their indices are the positions in
    /// [`Sketch::constraints`]; removed slots are counted).
    pub fn constraints_len(&self) -> usize {
        self.constraints.len()
    }

    /// Whether constraint `index` is enforced, suppressed or removed.
    pub fn constraint_state(&self, index: usize) -> ConstraintState {
        self.states.get(index).copied().unwrap_or_default()
    }

    /// Delete constraint `index`. Its slot stays (so other indices do not
    /// shift) but it is never enforced again. Returns `false` if it was
    /// already removed.
    pub fn remove_constraint(&mut self, index: usize) -> bool {
        self.set_state(index, ConstraintState::Removed)
    }

    /// Enforce constraint `index` (`true`) or keep it without enforcing it
    /// (`false`). A suppressed dimension is a reference dimension: nothing
    /// drives it, and [`Sketch::measure`] reads what it currently is.
    /// Returns `false`, changing nothing, for a removed constraint.
    pub fn set_constraint_enabled(&mut self, index: usize, enabled: bool) -> bool {
        if self.constraint_state(index) == ConstraintState::Removed {
            return false;
        }
        self.set_state(
            index,
            if enabled {
                ConstraintState::Active
            } else {
                ConstraintState::Suppressed
            },
        );
        true
    }

    /// Put a different constraint in slot `index` and make it active.
    pub fn replace_constraint(&mut self, index: usize, c: Constraint) {
        self.constraints[index] = c;
        self.states
            .resize(self.constraints.len(), ConstraintState::Active);
        self.states[index] = ConstraintState::Active;
        self.revision += 1;
    }

    /// The driving value of dimensional constraint `index` (see
    /// [`Constraint::value`]).
    pub fn constraint_value(&self, index: usize) -> Option<f64> {
        self.constraints[index].value()
    }

    /// Change the driving value of dimensional constraint `index` in place.
    /// The geometry is not touched until the next solve, which starts from
    /// where the sketch stands. Returns `false` for a constraint without a
    /// scalar value.
    ///
    /// This does not change which parameters any constraint touches, so
    /// [`Sketch::topology_revision`] stays the same.
    pub fn set_constraint_value(&mut self, index: usize, value: f64) -> bool {
        self.constraints[index].set_value(value)
    }

    /// What dimensional constraint `index` measures in the current
    /// geometry, whatever its state (see [`Constraint::measure`]).
    pub fn measure(&self, index: usize) -> Option<f64> {
        self.constraints[index].measure(&self.params, &self.entities)
    }

    fn set_state(&mut self, index: usize, state: ConstraintState) -> bool {
        self.states
            .resize(self.constraints.len(), ConstraintState::Active);
        if self.states[index] == state {
            return false;
        }
        self.states[index] = state;
        self.revision += 1;
        true
    }

    // ── Entity removal ───────────────────────────────────────────────────────

    /// Delete an entity together with every constraint that refers to it.
    /// Returns the indices of the constraints removed with it.
    ///
    /// The slot and the parameters of the entity stay in place (so other
    /// ids and parameter indices do not shift); the parameters are held and
    /// no longer count as freedoms.
    pub fn remove_entity(&mut self, id: EntityId) -> Vec<usize> {
        self.removed.resize(self.entities.len(), false);
        if self.removed[id] {
            return Vec::new();
        }
        self.removed[id] = true;
        for i in self.entity_param_range(id) {
            self.fixed[i] = true;
        }
        self.revision += 1;

        let dependents: Vec<usize> = (0..self.constraints.len())
            .filter(|&index| self.constraint_state(index) != ConstraintState::Removed)
            .filter(|&index| {
                self.constraints[index]
                    .references()
                    .iter()
                    .any(|&(entity, _)| entity == id)
            })
            .collect();
        for &index in &dependents {
            self.set_state(index, ConstraintState::Removed);
        }
        dependents
    }

    /// Whether entity `id` has been removed.
    pub fn is_entity_removed(&self, id: EntityId) -> bool {
        self.removed.get(id).copied().unwrap_or(false)
    }

    /// A counter that changes whenever the structure of the system does:
    /// an entity or constraint added or removed, a constraint enabled,
    /// suppressed or replaced, a parameter fixed. Changing a dimension value
    /// or moving geometry leaves it alone, so anything derived from the
    /// structure alone (the components, for one) can be kept while it holds.
    pub fn topology_revision(&self) -> u64 {
        self.revision
    }

    // ── Model check ──────────────────────────────────────────────────────────

    /// Check the model for defects that would make a solve meaningless or
    /// panic: constraints referring to entities that do not exist, were
    /// removed, or are of the wrong kind, and values or parameters that are
    /// not finite. Removed constraints are not checked.
    pub fn validate(&self) -> Result<(), Vec<SketchError>> {
        let mut errors = Vec::new();
        for (index, value) in self.params.iter().enumerate() {
            if !value.is_finite() && !self.is_param_of_removed_entity(index) {
                errors.push(SketchError::NonFiniteParameter { index });
            }
        }
        for (constraint, c) in self.constraints.iter().enumerate() {
            if self.constraint_state(constraint) == ConstraintState::Removed {
                continue;
            }
            for (entity, expected) in c.references() {
                match self.entities.get(entity) {
                    None => errors.push(SketchError::EntityOutOfRange { constraint, entity }),
                    Some(_) if self.is_entity_removed(entity) => {
                        errors.push(SketchError::EntityRemoved { constraint, entity })
                    }
                    Some(found) if !expected.contains(&found.kind) => {
                        errors.push(SketchError::WrongEntityKind {
                            constraint,
                            entity,
                            found: found.kind,
                            expected,
                        })
                    }
                    Some(_) => {}
                }
            }
            if !c.has_finite_values() {
                errors.push(SketchError::NonFiniteValue { constraint });
            }
        }
        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }

    fn is_param_of_removed_entity(&self, index: usize) -> bool {
        self.entities.iter().enumerate().any(|(id, e)| {
            self.is_entity_removed(id)
                && (e.param_start..e.param_start + e.kind.param_count()).contains(&index)
        })
    }

    /// The enforced constraints of this sketch, as the pipeline sees them.
    fn equations(&self) -> SketchEquations<'_> {
        SketchEquations::with_states(&self.entities, &self.constraints, &self.states)
    }

    // ── DOF analysis ─────────────────────────────────────────────────────────

    /// True degrees of freedom: free parameters − rank of the constraint
    /// Jacobian at the current parameters (the same account as
    /// [`Sketch::diagnose`]).
    ///
    /// - `dof > 0`: under-constrained (the sketch can still move).
    /// - `dof == 0`: no freedom left.
    ///
    /// Redundant equations consume no freedom, so the value is never
    /// negative and does not signal over-constraint by itself: read the
    /// redundant / conflicting sets from [`Sketch::diagnose`] for that.
    /// This runs a rank-revealing factorisation per connected component;
    /// [`Sketch::structural_dof_estimate`] is the cheap count.
    pub fn dof(&self) -> i64 {
        self.diagnose().dof
    }

    /// Count-based estimate: (free parameters) − (equations of the active
    /// constraints).
    ///
    /// No factorisation, but only an estimate: it is exact when every
    /// equation is independent, too low when some are redundant, and it
    /// counts across the whole sketch (a surplus in one place hides a
    /// freedom in another). A negative value does mean the system has
    /// redundant or conflicting equations.
    pub fn structural_dof_estimate(&self) -> i64 {
        let free = self.fixed.iter().filter(|&&f| !f).count() as i64;
        free - self.active_equation_count() as i64
    }

    fn active_equation_count(&self) -> usize {
        let eqs = self.equations();
        system::equation_total(&eqs, &system::all_indices(&eqs))
    }

    /// Detailed DOF analysis: free parameters, equations, Jacobian rank,
    /// the true DOF, the count-based estimate, and the free parameter
    /// indices.
    pub fn dof_analysis(&self) -> DofReport {
        let free_param_indices: Vec<usize> =
            (0..self.params.len()).filter(|&i| !self.fixed[i]).collect();
        let equations = self.active_equation_count();
        let diagnosis = self.diagnose();
        DofReport {
            free_params: free_param_indices.len(),
            equations,
            rank: diagnosis.rank,
            dof: diagnosis.dof,
            structural_dof: free_param_indices.len() as i64 - equations as i64,
            free_param_indices,
        }
    }

    // ── Solver ────────────────────────────────────────────────────────────────

    /// Run the Newton-Raphson GCS solver.
    ///
    /// Modifies entity parameters in-place.  Returns a [`SolveResult`]
    /// describing convergence. A connected component that does not converge
    /// gets its starting parameters back; use [`Sketch::solve_with`] to keep
    /// the last iterate instead.
    pub fn solve(&mut self) -> SolveResult {
        self.solve_with(&SolveOptions::default())
    }

    /// [`Sketch::solve`] with explicit [`SolveOptions`] (tolerance,
    /// iteration limit, what a failed solve leaves behind).
    pub fn solve_with(&mut self, options: &SolveOptions) -> SolveResult {
        self.solve_selected_with(None, options)
    }

    /// Rank-revealing diagnosis of the constraint system at the current
    /// parameters (true DOF, redundant / conflicting constraint indices).
    /// Suppressed and removed constraints take no part. See
    /// [`crate::diagnosis`].
    pub fn diagnose(&self) -> crate::diagnosis::Diagnosis {
        system::diagnose(&self.equations(), &self.params, &self.fixed, &NoGauge)
    }

    /// The state of every entity and constraint at the current parameters
    /// (see [`crate::status`]): one [`EntityStatus`](crate::EntityStatus)
    /// per entity — the most serious of conflicting, unsolved,
    /// over-defined, under-defined and well-defined that applies — and one
    /// [`ConstraintStatus`](crate::ConstraintStatus) per constraint. It
    /// runs [`Sketch::diagnose`] and [`Sketch::mobility`] (a factorisation
    /// each per component); call it after a solve, not every frame.
    pub fn status(&self) -> crate::status::Status {
        let movable = self.mobility();
        let movable_entity = |e: usize| self.entity_param_range(e).any(|p| movable[p]);
        let removed = |e: usize| self.is_entity_removed(e);
        let of_constraint = |c: usize| {
            self.constraints[c]
                .references()
                .into_iter()
                .map(|(e, _)| e)
                .collect()
        };
        crate::status::status_of(
            &self.equations(),
            &self.params,
            self.diagnose(),
            crate::status::Entities {
                count: self.entities.len(),
                removed: &removed,
                movable: &movable_entity,
                of_constraint: &of_constraint,
            },
        )
    }

    /// Connected components of the constraint graph: groups of constraint
    /// indices that share parameters (see [`crate::graph`]). The solver
    /// runs per component; localized re-solve = solving only the components
    /// an edit landed in. Suppressed and removed constraints belong to no
    /// component.
    pub fn components(&self) -> Vec<Vec<usize>> {
        system::components(&self.equations(), &self.params)
    }

    /// Localized solve: run Newton only on the selected components
    /// (ordinals into [`Sketch::components`]); `None` selects all.
    /// Unselected components are not iterated and their parameters stay
    /// untouched (solver roadmap M3.2).
    pub fn solve_selected(&mut self, only: Option<&[usize]>) -> SolveResult {
        self.solve_selected_with(only, &SolveOptions::default())
    }

    /// [`Sketch::solve_selected`] with explicit [`SolveOptions`].
    pub fn solve_selected_with(
        &mut self,
        only: Option<&[usize]>,
        options: &SolveOptions,
    ) -> SolveResult {
        let key = self.plan_key();
        let eqs = SketchEquations::with_states(&self.entities, &self.constraints, &self.states);
        let params = &self.params;
        let plans = self.plans.get(key, || system::Plans::new(&eqs, params));
        system::solve_planned(
            &eqs,
            &mut self.params,
            &self.fixed,
            &plans,
            only,
            &CoreOptions::from(options),
            |_, _| {},
        )
    }

    // ── Convenience point-ref helpers ─────────────────────────────────────────

    /// Return a [`PointRef`] for a `Point` entity.
    pub fn point(&self, id: EntityId) -> PointRef {
        debug_assert_eq!(self.entities[id].kind, EntityKind::Point);
        PointRef::Point(id)
    }

    /// Return a [`PointRef`] for the start of a `Line` entity.
    pub fn line_start(&self, id: EntityId) -> PointRef {
        debug_assert_eq!(self.entities[id].kind, EntityKind::Line);
        PointRef::LineStart(id)
    }

    /// Return a [`PointRef`] for the end of a `Line` entity.
    pub fn line_end(&self, id: EntityId) -> PointRef {
        debug_assert_eq!(self.entities[id].kind, EntityKind::Line);
        PointRef::LineEnd(id)
    }

    /// Return a [`PointRef`] for the center of a `Circle` or `Arc` entity.
    pub fn center(&self, id: EntityId) -> PointRef {
        debug_assert!(matches!(
            self.entities[id].kind,
            EntityKind::Circle | EntityKind::Arc
        ));
        PointRef::Center(id)
    }
}

impl Default for Sketch {
    fn default() -> Self {
        Self::new()
    }
}
