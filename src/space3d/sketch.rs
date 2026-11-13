//! 3D parametric sketch — container for space entities and constraints.
//!
//! Analogous to the 2D [`crate::Sketch`] but operates in 3D space with the
//! entity kinds of [`SpaceEntityKind`]: points, lines, planes, spheres,
//! cylinders, circles, arcs and Bezier splines.

use crate::numerics::{SolveOptions, SolveResult};
use crate::space3d::constraint::SpaceConstraint;
use crate::space3d::entity::{SpaceEntity, SpaceEntityId, SpaceEntityKind, SpaceFrame3};
use crate::space3d::solver::solve_space_with;

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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub struct SpaceSketch {
    pub params: Vec<f64>,
    pub entities: Vec<SpaceEntity>,
    pub constraints: Vec<SpaceConstraint>,
    fixed: Vec<bool>,
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
        }
    }

    /// Add a 3D point entity at the given coordinates.
    pub fn add_point(&mut self, x: f64, y: f64, z: f64) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
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

    /// Add a constraint.
    pub fn add_constraint(&mut self, constraint: SpaceConstraint) {
        self.constraints.push(constraint);
    }

    /// Fix a specific parameter so the solver won't change it.
    pub fn fix_param(&mut self, param_idx: usize) {
        if param_idx < self.fixed.len() {
            self.fixed[param_idx] = true;
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
    /// re-solve = solving only the components an edit landed in.
    pub fn components(&self) -> Vec<Vec<usize>> {
        crate::space3d::graph::constraint_components(
            &self.params,
            &self.entities,
            &self.constraints,
        )
    }

    /// Rank-revealing diagnosis of the constraint system at the current
    /// parameters (true DOF, redundant / conflicting constraint indices,
    /// dependency groups). See [`crate::space3d::diagnosis`].
    pub fn diagnose(&self) -> crate::space3d::diagnosis::SpaceDiagnosis {
        crate::space3d::diagnosis::diagnose(
            &self.params,
            &self.fixed,
            &self.entities,
            &self.constraints,
        )
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
        solve_space_with(
            &mut self.params,
            &self.fixed,
            &self.entities,
            &self.constraints,
            options,
        )
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
        let components = self.components();
        crate::space3d::solver::solve_components_with(
            &mut self.params,
            &self.fixed,
            &self.entities,
            &self.constraints,
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
        let components = self.components();
        let selected: Vec<usize> = components
            .iter()
            .enumerate()
            .filter(|(_, group)| group.iter().any(|c| changed_constraints.contains(c)))
            .map(|(ordinal, _)| ordinal)
            .collect();
        if selected.is_empty() {
            return crate::space3d::solver::residual_report(
                &self.params,
                &self.entities,
                &self.constraints,
                options.tolerance,
            );
        }
        crate::space3d::solver::solve_components_with(
            &mut self.params,
            &self.fixed,
            &self.entities,
            &self.constraints,
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

    /// Count-based estimate: (free parameters) − (constraint equations).
    ///
    /// No factorisation, but only an estimate: it counts raw storage
    /// parameters (a plane has four for three physical freedoms, a cylinder
    /// seven for five) and treats every equation as independent.
    pub fn structural_dof_estimate(&self) -> i64 {
        let free = self.fixed.iter().filter(|&&f| !f).count() as i64;
        let eqs: i64 = self
            .constraints
            .iter()
            .map(|c| c.equation_count() as i64)
            .sum();
        free - eqs
    }

    /// Detailed DOF analysis: free parameters, equations, Jacobian rank,
    /// the true DOF, the count-based estimate, and the free parameter
    /// indices.
    pub fn dof_analysis(&self) -> SpaceDofReport {
        let free_param_indices: Vec<usize> =
            (0..self.params.len()).filter(|&i| !self.fixed[i]).collect();
        let equations: usize = self.constraints.iter().map(|c| c.equation_count()).sum();
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
