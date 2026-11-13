//! 3D parametric sketch — container for space entities and constraints.
//!
//! Analogous to the 2D [`crate::Sketch`] but operates in 3D space with the
//! entity kinds of [`SpaceEntityKind`]: points, lines, planes, spheres,
//! cylinders, circles, arcs and Bezier splines.

use crate::edit::{ConstraintState, SketchError};
use crate::numerics::{SolveOptions, SolveResult};
use crate::space3d::constraint::SpaceConstraint;
use crate::space3d::entity::{SpaceEntity, SpaceEntityId, SpaceEntityKind, SpaceFrame3};
use crate::space3d::solver::SpaceEquations;
use crate::system;

/// A defect found by [`SpaceSketch::validate`].
pub type SpaceSketchError = SketchError<SpaceEntityKind>;

/// Detailed DOF analysis result for 3D sketches.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub struct SpaceDofReport {
    /// Number of free (non-fixed) parameters.
    pub free_params: usize,
    /// Total number of constraint equations.
    pub equations: usize,
    /// Rank of the constraint Jacobian over the free parameters at the
    /// current configuration.
    #[cfg_attr(feature = "serde", serde(default))]
    pub rank: usize,
    /// True degrees of freedom = effective free parameters − rank (the
    /// representation gauge of planes and cylinders is discounted).
    pub dof: i64,
    /// Count-based estimate = `free_params − equations`.
    #[cfg_attr(feature = "serde", serde(default))]
    pub structural_dof: i64,
    /// Indices of free parameters in the sketch's parameter vector.
    pub free_param_indices: Vec<usize>,
}

/// A 3D geometric constraint sketch.
///
/// Holds a flat parameter vector and a list of 3D entities and constraints.
/// Call [`SpaceSketch::solve`] to run the Newton-Raphson solver.
///
/// Entity ids and constraint indices are stable: removing something leaves
/// its slot behind instead of shifting what comes after it. Edit through
/// the methods ([`SpaceSketch::add_constraint`],
/// [`SpaceSketch::remove_constraint`], [`SpaceSketch::remove_entity`], …);
/// pushing into or deleting from the public lists directly goes unnoticed
/// by [`SpaceSketch::topology_revision`], and deleting from them shifts
/// indices.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub struct SpaceSketch {
    pub params: Vec<f64>,
    pub entities: Vec<SpaceEntity>,
    /// Constraint slots, including suppressed and removed ones (see
    /// [`SpaceSketch::constraint_state`]).
    pub constraints: Vec<SpaceConstraint>,
    fixed: Vec<bool>,
    /// The state of each constraint slot (a missing entry means active).
    #[cfg_attr(feature = "serde", serde(default))]
    states: Vec<ConstraintState>,
    /// Per entity: removed (a missing entry means present).
    #[cfg_attr(feature = "serde", serde(default))]
    removed: Vec<bool>,
    /// See [`SpaceSketch::topology_revision`].
    #[cfg_attr(feature = "serde", serde(default))]
    revision: u64,
}

impl Default for SpaceSketch {
    fn default() -> Self {
        Self::new()
    }
}

impl SpaceSketch {
    /// Create a new empty 3D sketch.
    pub fn new() -> Self {
        Self {
            params: Vec::new(),
            entities: Vec::new(),
            constraints: Vec::new(),
            fixed: Vec::new(),
            states: Vec::new(),
            removed: Vec::new(),
            revision: 0,
        }
    }

    /// Add a 3D point entity at the given coordinates.
    pub fn add_point(&mut self, x: f64, y: f64, z: f64) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        self.revision += 1;
        self.params.push(x);
        self.params.push(y);
        self.params.push(z);
        self.fixed.push(false);
        self.fixed.push(false);
        self.fixed.push(false);
        self.entities
            .push(SpaceEntity::new(SpaceEntityKind::SpacePoint, param_start));
        id
    }

    /// Add a 3D line entity from (x1,y1,z1) to (x2,y2,z2).
    pub fn add_line(
        &mut self,
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
    ) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        self.revision += 1;
        for &v in &[x1, y1, z1, x2, y2, z2] {
            self.params.push(v);
            self.fixed.push(false);
        }
        self.entities
            .push(SpaceEntity::new(SpaceEntityKind::SpaceLine, param_start));
        id
    }

    /// Add a plane entity with normal (nx,ny,nz) and distance d from origin.
    ///
    /// The `(normal, offset)` pair is stored in the fixed unit-normal
    /// representation: the normal is normalized on insertion (and the offset
    /// scaled with it), removing the redundant scale coordinate from the
    /// parameter vector. Uniform scaling leaves every constraint residual
    /// unchanged, so this does not alter the geometry.
    pub fn add_plane(&mut self, nx: f64, ny: f64, nz: f64, d: f64) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        self.revision += 1;
        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        let inv = if len > crate::numerics::SOLVER_NORM_TOL {
            1.0 / len
        } else {
            1.0
        };
        for &v in &[nx * inv, ny * inv, nz * inv, d * inv] {
            self.params.push(v);
            self.fixed.push(false);
        }
        self.entities
            .push(SpaceEntity::new(SpaceEntityKind::Plane, param_start));
        id
    }

    /// Add a sphere entity at (cx,cy,cz) with radius r.
    pub fn add_sphere(&mut self, cx: f64, cy: f64, cz: f64, r: f64) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        self.revision += 1;
        for &v in &[cx, cy, cz, r] {
            self.params.push(v);
            self.fixed.push(false);
        }
        self.entities
            .push(SpaceEntity::new(SpaceEntityKind::Sphere, param_start));
        id
    }

    /// Add a cylinder entity with axis point (cx,cy,cz), axis direction (ax,ay,az), and radius r.
    #[allow(clippy::too_many_arguments)]
    pub fn add_cylinder(
        &mut self,
        cx: f64,
        cy: f64,
        cz: f64,
        ax: f64,
        ay: f64,
        az: f64,
        r: f64,
    ) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        self.revision += 1;
        for &v in &[cx, cy, cz, ax, ay, az, r] {
            self.params.push(v);
            self.fixed.push(false);
        }
        self.entities
            .push(SpaceEntity::new(SpaceEntityKind::Cylinder, param_start));
        id
    }

    /// Add a circle entity (center, radius) over an explicit in-plane frame.
    ///
    /// The frame is fixed metadata (datum orientation), so the circle has
    /// four free parameters — center (3) plus radius (1) — and never carries
    /// a degenerate normal in the parameter vector. Degenerate input is
    /// sanitized deterministically by [`SpaceFrame3::sanitized`]; a
    /// non-positive or non-finite radius is pinned to the positive solver
    /// floor so every downstream residual stays finite.
    pub fn add_circle(
        &mut self,
        cx: f64,
        cy: f64,
        cz: f64,
        r: f64,
        frame: SpaceFrame3,
    ) -> SpaceEntityId {
        self.add_curve(SpaceEntityKind::SpaceCircle, cx, cy, cz, r, frame)
    }

    /// Add a circular arc entity (center, radius) over an explicit in-plane
    /// frame; the angles are carried by the frame's `start_angle` and
    /// `end_angle` (the end is re-expressed as `start + span` on every
    /// evaluation, so a wrapped angle pair is interpreted as the
    /// counter-clockwise span).
    pub fn add_arc(
        &mut self,
        cx: f64,
        cy: f64,
        cz: f64,
        r: f64,
        frame: SpaceFrame3,
    ) -> SpaceEntityId {
        self.add_curve(SpaceEntityKind::SpaceArc, cx, cy, cz, r, frame)
    }

    fn add_curve(
        &mut self,
        kind: SpaceEntityKind,
        cx: f64,
        cy: f64,
        cz: f64,
        r: f64,
        frame: SpaceFrame3,
    ) -> SpaceEntityId {
        // A degenerate frame can never enter the entity metadata: the
        // sanitizer normalizes/orthogonalizes deterministically, so every
        // residual over this curve sees a unit orthogonal basis.
        let frame = SpaceFrame3::sanitized(
            frame.normal,
            frame.x_dir,
            frame.start_angle,
            frame.end_angle,
        );
        let id = self.entities.len();
        let param_start = self.params.len();
        self.revision += 1;
        let center = if cx.is_finite() && cy.is_finite() && cz.is_finite() {
            [cx, cy, cz]
        } else {
            [0.0, 0.0, 0.0]
        };
        let radius = if r.is_finite() && r > crate::numerics::SOLVER_NORM_TOL {
            r
        } else {
            crate::numerics::SOLVER_NORM_TOL
        };
        for &v in &[center[0], center[1], center[2], radius] {
            self.params.push(v);
            self.fixed.push(false);
        }
        self.entities
            .push(SpaceEntity::with_frame(kind, param_start, frame));
        id
    }

    /// Center and radius of a circle/arc entity.
    pub fn curve_geometry(&self, id: SpaceEntityId) -> Option<(f64, f64, f64, f64)> {
        let e = self.entities.get(id)?;
        if !e.kind.is_curve() {
            return None;
        }
        Some((
            self.params[e.param(0)],
            self.params[e.param(1)],
            self.params[e.param(2)],
            self.params[e.param(3)],
        ))
    }

    /// Fixed frame of a circle/arc entity.
    pub fn curve_frame(&self, id: SpaceEntityId) -> Option<SpaceFrame3> {
        let e = self.entities.get(id)?;
        if !e.kind.is_curve() {
            return None;
        }
        Some(e.curve_frame())
    }

    /// Add a Bezier spline entity from control points.
    ///
    /// Each control point contributes 3 parameters (x, y, z). Minimum 4
    /// control points (cubic Bezier). Returns the entity id.
    pub fn add_spline(&mut self, control_points: &[(f64, f64, f64)]) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        self.revision += 1;
        let n = control_points.len();
        for &(x, y, z) in control_points {
            self.params.push(x);
            self.params.push(y);
            self.params.push(z);
            self.fixed.push(false);
            self.fixed.push(false);
            self.fixed.push(false);
        }
        self.entities.push(SpaceEntity::with_control_points(
            SpaceEntityKind::SpaceSpline,
            param_start,
            n,
        ));
        id
    }

    /// Control points of a spline entity.
    pub fn spline_control_points(&self, id: SpaceEntityId) -> Option<Vec<(f64, f64, f64)>> {
        let e = self.entities.get(id)?;
        if e.kind != SpaceEntityKind::SpaceSpline {
            return None;
        }
        let n = e.control_point_count;
        Some(
            (0..n)
                .map(|i| {
                    (
                        self.params[e.param(i * 3)],
                        self.params[e.param(i * 3 + 1)],
                        self.params[e.param(i * 3 + 2)],
                    )
                })
                .collect(),
        )
    }

    // ── Constraint management ─────────────────────────────────────────────────

    /// Add a constraint. Returns its index, which stays valid for the life
    /// of the sketch: removing other constraints does not shift it.
    pub fn add_constraint(&mut self, constraint: SpaceConstraint) -> usize {
        self.constraints.push(constraint);
        self.states
            .resize(self.constraints.len(), ConstraintState::Active);
        self.revision += 1;
        self.constraints.len() - 1
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
    /// drives it, and [`SpaceSketch::measure`] reads what it currently is.
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
    pub fn replace_constraint(&mut self, index: usize, constraint: SpaceConstraint) {
        self.constraints[index] = constraint;
        self.states
            .resize(self.constraints.len(), ConstraintState::Active);
        self.states[index] = ConstraintState::Active;
        self.revision += 1;
    }

    /// The driving value of dimensional constraint `index` (see
    /// [`SpaceConstraint::value`]).
    pub fn constraint_value(&self, index: usize) -> Option<f64> {
        self.constraints[index].value()
    }

    /// Change the driving value of dimensional constraint `index` in place.
    /// The geometry is not touched until the next solve, which starts from
    /// where the sketch stands. Returns `false` for a constraint without a
    /// scalar value.
    ///
    /// This does not change which parameters any constraint touches, so
    /// [`SpaceSketch::topology_revision`] stays the same.
    pub fn set_constraint_value(&mut self, index: usize, value: f64) -> bool {
        self.constraints[index].set_value(value)
    }

    /// What dimensional constraint `index` measures in the current
    /// geometry, whatever its state (see [`SpaceConstraint::measure`]).
    pub fn measure(&self, index: usize) -> Option<f64> {
        self.constraints[index].measure(&self.params, &self.entities)
    }

    fn set_state(&mut self, index: usize, state: ConstraintState) -> bool {
        assert!(
            index < self.constraints.len(),
            "constraint index {index} out of range"
        );
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
    pub fn remove_entity(&mut self, id: SpaceEntityId) -> Vec<usize> {
        self.removed.resize(self.entities.len(), false);
        if self.removed[id] {
            return Vec::new();
        }
        self.removed[id] = true;
        let entity = &self.entities[id];
        for i in entity.param_start..entity.param_start + entity.param_count() {
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
    pub fn is_entity_removed(&self, id: SpaceEntityId) -> bool {
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
    /// removed, or are of the wrong kind; spline control points that do not
    /// exist; and values or parameters that are not finite. Removed
    /// constraints are not checked.
    pub fn validate(&self) -> Result<(), Vec<SpaceSketchError>> {
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
            if let Some((entity, indices)) = c.control_point_indices() {
                // A reference that is not a spline was reported above.
                let spline = self
                    .entities
                    .get(entity)
                    .filter(|e| e.kind == SpaceEntityKind::SpaceSpline);
                if let Some(spline) = spline {
                    let count = spline.control_point_count;
                    for index in indices.into_iter().filter(|&index| index >= count) {
                        errors.push(SketchError::ControlPointOutOfRange {
                            constraint,
                            entity,
                            index,
                            count,
                        });
                    }
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
                && (e.param_start..e.param_start + e.param_count()).contains(&index)
        })
    }

    /// The enforced constraints of this sketch, as the pipeline sees them.
    fn equations(&self) -> SpaceEquations<'_> {
        SpaceEquations::with_states(&self.entities, &self.constraints, &self.states)
    }

    fn active_equation_count(&self) -> usize {
        let eqs = self.equations();
        system::equation_total(&eqs, &system::all_indices(&eqs))
    }

    // ── Fixing ────────────────────────────────────────────────────────────────

    /// Fix a specific parameter so the solver won't change it.
    pub fn fix_param(&mut self, param_idx: usize) {
        if param_idx < self.fixed.len() {
            self.fixed[param_idx] = true;
            self.revision += 1;
        }
    }

    /// Fix all parameters of an entity.
    pub fn fix_entity(&mut self, entity_id: SpaceEntityId) {
        if let Some(entity) = self.entities.get(entity_id) {
            // The per-entity count: a spline's parameter count depends on
            // its control point count (the kind-level count reports 0).
            let count = entity.param_count();
            for i in 0..count {
                let idx = entity.param_start + i;
                if idx < self.fixed.len() {
                    self.fixed[idx] = true;
                }
            }
            self.revision += 1;
        }
    }

    /// Get the coordinates of a SpacePoint entity.
    pub fn point_coords(&self, id: SpaceEntityId) -> Option<(f64, f64, f64)> {
        let e = self.entities.get(id)?;
        if e.kind != SpaceEntityKind::SpacePoint {
            return None;
        }
        Some((
            self.params[e.param(0)],
            self.params[e.param(1)],
            self.params[e.param(2)],
        ))
    }

    /// Get the endpoints of a SpaceLine entity.
    pub fn line_endpoints(&self, id: SpaceEntityId) -> Option<(f64, f64, f64, f64, f64, f64)> {
        let e = self.entities.get(id)?;
        if e.kind != SpaceEntityKind::SpaceLine {
            return None;
        }
        Some((
            self.params[e.param(0)],
            self.params[e.param(1)],
            self.params[e.param(2)],
            self.params[e.param(3)],
            self.params[e.param(4)],
            self.params[e.param(5)],
        ))
    }

    // ── Parameter access ──────────────────────────────────────────────────────

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

    /// The solver's fixed mask (one flag per parameter).
    pub fn fixed_mask(&self) -> &[bool] {
        &self.fixed
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

    /// Re-pin every fully free plane to the unit-normal representation.
    /// Residual-invariant gauge fix (see [`crate::space3d::solver`]).
    pub fn normalize_planes(&mut self) {
        crate::space3d::solver::normalize_planes(&mut self.params, &self.entities, &self.fixed);
    }

    // ── Structural analysis ───────────────────────────────────────────────────

    /// Connected components of this sketch's constraint graph: groups of
    /// constraint indices that share parameters (see
    /// [`crate::space3d::graph`]). The solver runs per component; localized
    /// re-solve = solving only the components an edit landed in. Suppressed
    /// and removed constraints belong to no component.
    pub fn components(&self) -> Vec<Vec<usize>> {
        system::components(&self.equations(), &self.params)
    }

    /// Rank-revealing diagnosis of the constraint system at the current
    /// parameters (true DOF, redundant / conflicting constraint indices,
    /// dependency groups). Suppressed and removed constraints take no
    /// part. See [`crate::space3d::diagnosis`].
    pub fn diagnose(&self) -> crate::space3d::diagnosis::SpaceDiagnosis {
        crate::space3d::diagnosis::diagnose_equations(&self.equations(), &self.params, &self.fixed)
    }

    // ── Solver ────────────────────────────────────────────────────────────────

    /// Run the Newton-Raphson solver on this sketch (per connected
    /// component; the current parameters are the warm start).
    ///
    /// A component that does not converge gets its starting parameters
    /// back; use [`SpaceSketch::solve_with`] to keep the last iterate
    /// instead.
    pub fn solve(&mut self) -> SolveResult {
        self.solve_with(&SolveOptions::default())
    }

    /// [`SpaceSketch::solve`] with explicit [`SolveOptions`] (tolerance,
    /// iteration limit, what a failed solve leaves behind).
    pub fn solve_with(&mut self, options: &SolveOptions) -> SolveResult {
        // Gauge fix before iterating: the solver's variable space is the
        // (unit normal, offset) representation.
        self.normalize_planes();
        self.solve_selected_with(None, options)
    }

    /// Localized solve: run Newton only on the selected components
    /// (ordinals into [`SpaceSketch::components`]); `None` selects all.
    /// Unselected components are not iterated and their parameters stay
    /// untouched.
    pub fn solve_selected(&mut self, only: Option<&[usize]>) -> SolveResult {
        self.solve_selected_with(only, &SolveOptions::default())
    }

    /// [`SpaceSketch::solve_selected`] with explicit [`SolveOptions`].
    pub fn solve_selected_with(
        &mut self,
        only: Option<&[usize]>,
        options: &SolveOptions,
    ) -> SolveResult {
        let eqs = SpaceEquations::with_states(&self.entities, &self.constraints, &self.states);
        let components = system::components(&eqs, &self.params);
        crate::space3d::solver::solve_equations(
            &eqs,
            &mut self.params,
            &self.fixed,
            &components,
            only,
            options,
        )
    }

    /// Incremental re-solve: iterate only the components that contain one of
    /// the `changed_constraints` (indices into [`SpaceSketch::constraints`]),
    /// leaving every other component's parameters untouched. An empty (or
    /// unmatched) selection does not iterate and reports the current
    /// residual state instead — the warm start is the parameter vector as it
    /// stands.
    pub fn solve_incremental(&mut self, changed_constraints: &[usize]) -> SolveResult {
        self.solve_incremental_with(changed_constraints, &SolveOptions::default())
    }

    /// [`SpaceSketch::solve_incremental`] with explicit [`SolveOptions`].
    pub fn solve_incremental_with(
        &mut self,
        changed_constraints: &[usize],
        options: &SolveOptions,
    ) -> SolveResult {
        let eqs = SpaceEquations::with_states(&self.entities, &self.constraints, &self.states);
        let components = system::components(&eqs, &self.params);
        let selected: Vec<usize> = components
            .iter()
            .enumerate()
            .filter(|(_, group)| group.iter().any(|c| changed_constraints.contains(c)))
            .map(|(ordinal, _)| ordinal)
            .collect();
        if selected.is_empty() {
            return crate::space3d::solver::residual_report(&eqs, &self.params, options.tolerance);
        }
        crate::space3d::solver::solve_equations(
            &eqs,
            &mut self.params,
            &self.fixed,
            &components,
            Some(&selected),
            options,
        )
    }

    /// True degrees of freedom: effective free parameters − rank of the
    /// constraint Jacobian at the current parameters (the same account as
    /// [`SpaceSketch::diagnose`], which discounts the representation gauge
    /// of planes and cylinders).
    ///
    /// Redundant equations consume no freedom, so the value does not signal
    /// over-constraint by itself: read the redundant / conflicting sets
    /// from [`SpaceSketch::diagnose`] for that.
    /// [`SpaceSketch::structural_dof_estimate`] is the cheap count.
    pub fn dof(&self) -> i64 {
        self.diagnose().dof
    }

    /// Count-based estimate: (free parameters) − (equations of the active
    /// constraints).
    ///
    /// No factorisation, but only an estimate: it counts raw storage
    /// parameters (a plane has four for three physical freedoms, a cylinder
    /// seven for five) and treats every equation as independent.
    pub fn structural_dof_estimate(&self) -> i64 {
        let free = self.fixed.iter().filter(|&&f| !f).count() as i64;
        free - self.active_equation_count() as i64
    }

    /// Detailed DOF analysis: free parameters, equations, Jacobian rank,
    /// the true DOF, the count-based estimate, and the free parameter
    /// indices.
    pub fn dof_analysis(&self) -> SpaceDofReport {
        let free_param_indices: Vec<usize> =
            (0..self.params.len()).filter(|&i| !self.fixed[i]).collect();
        let equations = self.active_equation_count();
        let diagnosis = self.diagnose();
        SpaceDofReport {
            free_params: free_param_indices.len(),
            equations,
            rank: diagnosis.rank,
            dof: diagnosis.dof,
            structural_dof: free_param_indices.len() as i64 - equations as i64,
            free_param_indices,
        }
    }
}
