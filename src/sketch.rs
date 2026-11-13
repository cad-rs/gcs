//! The main [`Sketch`] API.

use glam::DVec2;

use crate::constraint::Constraint;
use crate::entity::{Entity, EntityId, EntityKind, PointRef};
use crate::solver::{self, SolveResult};

/// Detailed DOF analysis result.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub struct DofReport {
    /// Number of free (non-fixed) parameters.
    pub free_params: usize,
    /// Total number of constraint equations.
    pub equations: usize,
    /// Net degrees of freedom = free_params - equations.
    pub dof: i64,
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
    /// Constraint list.
    pub(crate) constraints: Vec<Constraint>,
}

impl Sketch {
    /// Create an empty sketch.
    pub fn new() -> Self {
        Self {
            params: Vec::new(),
            fixed: Vec::new(),
            entities: Vec::new(),
            constraints: Vec::new(),
        }
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
        self.fixed.extend(std::iter::repeat(false).take(init.len()));
        self.entities.push(Entity::new(kind, param_start));
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
    /// [`crate::mobility::mobility`]).
    pub fn mobility(&self) -> Vec<bool> {
        crate::mobility::mobility(&self.params, &self.fixed, &self.entities, &self.constraints)
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
    }

    /// Fix all parameters of an entity (equivalent to a rigid body constraint).
    pub fn fix_entity(&mut self, id: EntityId) {
        let e = &self.entities[id];
        for i in e.param_start..e.param_start + e.kind.param_count() {
            self.fixed[i] = true;
        }
    }

    // ── Constraint management ─────────────────────────────────────────────────

    /// Add a constraint.
    pub fn add_constraint(&mut self, c: Constraint) {
        self.constraints.push(c);
    }

    /// The constraint list.
    pub fn constraints(&self) -> &[Constraint] {
        &self.constraints
    }

    /// Number of constraints (their gcs indices are the positions here —
    /// the mapping target for callers that track edited constraints).
    pub fn constraints_len(&self) -> usize {
        self.constraints.len()
    }

    // ── DOF analysis ─────────────────────────────────────────────────────────

    /// Degrees of freedom = (free parameters) − (constraint equations).
    ///
    /// - `dof > 0`: under-constrained (sketch can still move).
    /// - `dof == 0`: fully constrained.
    /// - `dof < 0`: over-constrained (conflicting constraints).
    pub fn dof(&self) -> i64 {
        let free = self.fixed.iter().filter(|&&f| !f).count() as i64;
        let eqs: i64 = self
            .constraints
            .iter()
            .map(|c| c.equation_count() as i64)
            .sum();
        free - eqs
    }

    /// Detailed DOF analysis.
    ///
    /// Returns the number of free parameters, constraint equations,
    /// the net DOF, and the indices of free parameters.
    pub fn dof_analysis(&self) -> DofReport {
        let free_param_indices: Vec<usize> =
            (0..self.params.len()).filter(|&i| !self.fixed[i]).collect();
        let free = free_param_indices.len() as i64;
        let eqs: i64 = self
            .constraints
            .iter()
            .map(|c| c.equation_count() as i64)
            .sum();
        DofReport {
            free_params: free_param_indices.len(),
            equations: eqs as usize,
            dof: free - eqs,
            free_param_indices,
        }
    }

    // ── Solver ────────────────────────────────────────────────────────────────

    /// Run the Newton-Raphson GCS solver.
    ///
    /// Modifies entity parameters in-place.  Returns a [`SolveResult`]
    /// describing convergence.
    pub fn solve(&mut self) -> SolveResult {
        solver::solve(
            &mut self.params,
            &self.fixed,
            &self.entities,
            &self.constraints,
        )
    }

    /// Rank-revealing diagnosis of the constraint system at the current
    /// parameters (true DOF, redundant / conflicting constraint indices).
    /// See [`crate::diagnosis`].
    pub fn diagnose(&self) -> crate::diagnosis::Diagnosis {
        crate::diagnosis::diagnose(&self.params, &self.fixed, &self.entities, &self.constraints)
    }

    /// Connected components of this sketch's constraint graph: groups of
    /// constraint indices that share parameters (see [`crate::graph`]).
    /// The solver runs per component; localized re-solve = solving only
    /// the components an edit landed in.
    pub fn components(&self) -> Vec<Vec<usize>> {
        crate::graph::constraint_components(&self.params, &self.entities, &self.constraints)
    }

    /// Localized solve: run Newton only on the selected components
    /// (ordinals into [`Sketch::components`]); `None` selects all.
    /// Unselected components are not iterated and their parameters stay
    /// untouched (solver roadmap M3.2).
    pub fn solve_selected(&mut self, only: Option<&[usize]>) -> SolveResult {
        let components = self.components();
        crate::solver::solve_components(
            &mut self.params,
            &self.fixed,
            &self.entities,
            &self.constraints,
            &components,
            only,
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
