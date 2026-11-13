//! Rigid-pose solve: the assembly side of the shared pipeline.
//!
//! Per component: translation `(tx, ty, tz)` plus incremental orientation via
//! quaternion integration (`q ← normalize(q · dq)`). Mate offsets `point1` /
//! `point2` are in component space; world points use `t + R·local`. Direction
//! mates (`Parallel`, `Angle`, `Perpendicular`) use unit vectors
//! `normalize(R·local_dir)`. Jacobian rows are analytic (translation ±I,
//! rotation via the `ω × r` linearization, direction rows via `δd = ω × d`)
//! except the orientation rows of `Fixed` and `Lock`, which use central
//! finite differences.
//!
//! This module supplies the equations (residuals and Jacobian rows per
//! constraint) and the solver front end. The iteration itself — the
//! normal-equation step, the line search, the convergence polish — is the
//! Newton core the other families use, driven per connected cluster. A body
//! that a `Fixed` constraint pins in full is grounded: placed at its pose and
//! left out of the unknowns, so the clusters hanging off it solve
//! independently.

use glam::{DQuat, DVec3, EulerRot};

use super::constraint::{Constraint, ConstraintKind, ConstraintSystem};
use super::equations::{pack, unpack, RigidEquations, COLUMNS_PER_BODY};
use crate::numerics::{CoreOptions, FailurePolicy, ResidualNorm};
use crate::system;

/// Finite-difference step for direction-mate Jacobian rows.
const FD_EPS: f64 = 1e-7;

/// Result of solving.
#[derive(Debug, Clone)]
pub struct PoseSolveResult {
    /// Whether solving succeeded
    pub success: bool,
    /// Final error value
    pub error: f64,
    /// Number of iterations
    pub iterations: usize,
    /// Solved positions per component (translation)
    pub positions: Vec<DVec3>,
    /// Solved orientations as normalized quaternions (same order as [`Self::positions`])
    pub rotations: Vec<DQuat>,
    /// Solved orientations per component (Euler XYZ in radians, from the internal quaternion)
    pub orientations: Vec<DVec3>,
}

/// Solver options.
#[derive(Debug, Clone)]
pub struct SolveOptions {
    /// Maximum iterations
    pub max_iterations: usize,
    /// Convergence tolerance (L2 norm of residual vector)
    pub tolerance: f64,
    /// Scale factor applied to each Gauss–Newton step.
    pub damping: f64,
    /// Ridge added to the diagonal of `JᵀJ` for Levenberg–Marquardt-style stability.
    pub lm_lambda: f64,
}

impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance: 1e-6,
            damping: 1.0,
            lm_lambda: 1e-10,
        }
    }
}

/// Constraint solver.
#[derive(Debug, Clone)]
pub struct ConstraintSolver {
    options: SolveOptions,
}

/// Dense residual/Jacobian containers: `j[row][col]` with one column per
/// rigid-DoF (`6 · num_components` columns).
pub(crate) struct ResidualJacobian {
    pub(crate) r: Vec<f64>,
    pub(crate) j: Vec<Vec<f64>>,
}

impl ConstraintSolver {
    /// Create a new solver with default options.
    pub fn new() -> Self {
        Self {
            options: SolveOptions::default(),
        }
    }

    /// Create solver with custom options.
    pub fn with_options(options: SolveOptions) -> Self {
        Self { options }
    }

    /// Solve the constraint system using Gauss–Newton (`6 · num_components` unknowns).
    /// Fixed bodies seed the start; everything else begins at zero/identity.
    pub fn solve(&self, system: &ConstraintSystem) -> PoseSolveResult {
        let (positions, rots) = Self::seed_fixed(system);
        self.solve_from(system, positions, rots)
    }

    /// Seed fixed bodies from their constraints.
    pub fn seed_fixed(system: &ConstraintSystem) -> (Vec<DVec3>, Vec<DQuat>) {
        let n = system.num_components;
        let mut positions: Vec<DVec3> = (0..n).map(|_| DVec3::ZERO).collect();
        let mut rots: Vec<DQuat> = (0..n).map(|_| DQuat::IDENTITY).collect();
        for constraint in &system.constraints {
            if let ConstraintKind::Fixed {
                body,
                position,
                orientation,
            } = constraint.kind
            {
                positions[body] = position;
                if let Some(q) = orientation {
                    rots[body] = q.normalize();
                }
            }
        }
        (positions, rots)
    }

    /// Connected components over the BODY graph: constraints join the
    /// bodies they reference (solver roadmap M3, assembly edition).
    /// Groups of constraint indices, ordered by first appearance.
    pub fn components(system: &ConstraintSystem) -> Vec<Vec<usize>> {
        let n = system.num_components;
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                x = parent[x];
            }
            x
        }
        // Pass 1: union the bodies of every constraint.
        for c in &system.constraints {
            let (first, second) = c.kind.bodies();
            let a = find(&mut parent, first);
            if let Some(c2) = second {
                let b = find(&mut parent, c2);
                if a != b {
                    parent[a] = b;
                }
            }
        }
        // Pass 2: group constraints by their FINAL root. Grouping during
        // pass 1 would key earlier constraints by stale roots when a later
        // constraint merges two established groups (regression-tested).
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut first_root: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        for (idx, c) in system.constraints.iter().enumerate() {
            let (first, second) = c.kind.bodies();
            let mut root = find(&mut parent, first);
            if let Some(c2) = second {
                root = find(&mut parent, c2);
            }
            let slot = first_root.entry(root).or_insert_with(|| {
                groups.push(Vec::new());
                groups.len() - 1
            });
            groups[*slot].push(idx);
        }
        groups
    }

    /// Localized assembly solve: run Newton only on the selected body
    /// components (ordinals into [`Self::components`]); untouched clusters
    /// are not iterated and their poses stay untouched (M3, assembly
    /// edition).
    pub fn solve_components(
        &self,
        system: &ConstraintSystem,
        components: &[Vec<usize>],
        only: Option<&[usize]>,
    ) -> PoseSolveResult {
        let (positions, rots) = Self::seed_fixed(system);
        self.solve_components_from(system, components, only, positions, rots)
    }

    /// [`Self::solve_components`] from an explicit warm start — the
    /// chaining form for piecewise multi-component solves.
    ///
    /// A body that a `Fixed` constraint pins in full (position and
    /// orientation) is grounded: it is placed at that pose and is not an
    /// unknown. Mates to the ground therefore couple nothing, and a
    /// selected component is solved as the independent sub-clusters it
    /// falls into once the ground is taken out — a hub with forty parts
    /// pinned to it is forty small solves, not one large one.
    ///
    /// Each sub-cluster is one Newton solve of the shared pipeline over its
    /// own bodies (six columns each): the normal-equation step, the line
    /// search and the convergence polish are the ones the 2D and
    /// 3D-geometric families use. A sub-cluster converges when the L2 norm
    /// of its own residuals is below the tolerance; `success` says that
    /// every one did, and `error` is the L2 norm of the whole system's
    /// residuals at the returned poses.
    pub fn solve_components_from(
        &self,
        system: &ConstraintSystem,
        components: &[Vec<usize>],
        only: Option<&[usize]>,
        positions: Vec<DVec3>,
        rots: Vec<DQuat>,
    ) -> PoseSolveResult {
        // Grounded bodies are known, not unknowns: place them, and let the
        // clusters that hang off them come apart.
        let eqs = RigidEquations::grounded(system);
        let mut params = pack(&positions, &rots);
        let mut clusters: Vec<Vec<usize>> = Vec::new();
        for (ordinal, group) in components.iter().enumerate() {
            if only.is_some_and(|selected| !selected.contains(&ordinal)) {
                continue;
            }
            eqs.snap_grounded(group, &mut params);
            clusters.extend(system::components_of(&eqs, group, &params));
        }

        let fixed = vec![false; system.num_components * COLUMNS_PER_BODY];
        let result = system::solve_components(
            &eqs,
            &mut params,
            &fixed,
            &clusters,
            None,
            &self.core_options(),
            |_, _| {},
        );

        // The reported error counts every constraint as an equation again,
        // the grounding ones included.
        let error = eqs.residual_norm_of_every_constraint(&params);
        let (positions, rots) = unpack(&params);
        PoseSolveResult {
            success: result.converged,
            error,
            iterations: result.iterations,
            orientations: quats_to_euler_xyz(&rots),
            positions,
            rotations: rots,
        }
    }

    /// Solve from an explicit warm start (poses are NOT re-seeded) — the
    /// M4-style round-trip and drift-probe entry point.
    ///
    /// Solves every connected cluster on its own (a full solve is the
    /// composition of the per-cluster solves, as in the other families). A
    /// body no constraint references keeps its pose untouched.
    pub fn solve_from(
        &self,
        system: &ConstraintSystem,
        positions: Vec<DVec3>,
        rots: Vec<DQuat>,
    ) -> PoseSolveResult {
        let components = Self::components(system);
        self.solve_components_from(system, &components, None, positions, rots)
    }

    /// This solver's options in the Newton core's terms. A solve that does
    /// not converge returns the poses it reached, as it always has.
    fn core_options(&self) -> CoreOptions {
        CoreOptions {
            tolerance: self.options.tolerance,
            norm: ResidualNorm::L2,
            max_iterations: self.options.max_iterations,
            on_failure: FailurePolicy::KeepLast,
            lambda: self.options.lm_lambda.max(1e-16),
            initial_step: self.options.damping,
            // One polish step, as this family has always done.
            polish_steps: 1,
            polish_singular: true,
        }
    }

    /// Equation count of a single constraint (the stacked-row width the
    /// diagnosis uses to map rows back to mate ids).
    pub(crate) fn equation_count_of(c: &Constraint) -> usize {
        match c.kind {
            ConstraintKind::Fixed { orientation, .. } => {
                3 + if orientation.is_some() { 3 } else { 0 }
            }
            ConstraintKind::Coincident { .. } => 3,
            ConstraintKind::Distance { .. } => 1,
            ConstraintKind::Parallel { .. } => 2,
            ConstraintKind::Angle { .. } => 1,
            ConstraintKind::Perpendicular { .. } => 1,
            ConstraintKind::Concentric { .. } => 4,
            ConstraintKind::PointOnPlane { .. } => 1,
            ConstraintKind::Lock { .. } => 6,
        }
    }

    pub(crate) fn equation_count(system: &ConstraintSystem) -> usize {
        system.constraints.iter().map(Self::equation_count_of).sum()
    }

    #[cfg(test)]
    pub(crate) fn build_residual_jacobian(
        system: &ConstraintSystem,
        positions: &[DVec3],
        rots: &[DQuat],
    ) -> ResidualJacobian {
        let mut rj = ResidualJacobian {
            r: Vec::new(),
            j: Vec::new(),
        };
        Self::fill_residual_jacobian(system, positions, rots, &mut rj);
        rj
    }

    /// [`Self::build_residual_jacobian`] into the buffers of `rj`, which are
    /// resized and overwritten.
    pub(crate) fn fill_residual_jacobian(
        system: &ConstraintSystem,
        positions: &[DVec3],
        rots: &[DQuat],
        rj: &mut ResidualJacobian,
    ) {
        let n = system.num_components;
        let cols = n * 6;
        let m = Self::equation_count(system);

        rj.r.clear();
        rj.r.resize(m, 0.0_f64);
        rj.j.resize_with(m, Vec::new);
        for row in rj.j.iter_mut() {
            row.clear();
            row.resize(cols, 0.0_f64);
        }

        let mut row = 0usize;
        for constraint in &system.constraints {
            match constraint.kind {
                ConstraintKind::Fixed {
                    body: c,
                    position: t,
                    orientation,
                } => {
                    let p = positions[c];
                    for k in 0..3 {
                        rj.r[row + k] = [p.x - t.x, p.y - t.y, p.z - t.z][k];
                        rj.j[row + k][c * 6 + k] = 1.0;
                    }
                    row += 3;
                    if let Some(q_tgt) = orientation {
                        let res = fixed_orientation_residual(rots[c], q_tgt);
                        rj.r[row] = res.x;
                        rj.r[row + 1] = res.y;
                        rj.r[row + 2] = res.z;
                        Self::fill_vec3_residual_fd(
                            |_, rp| fixed_orientation_residual(rp[c], q_tgt),
                            positions,
                            rots,
                            row,
                            cols,
                            &mut rj.j,
                        );
                        row += 3;
                    }
                }
                ConstraintKind::Coincident {
                    body1: c1,
                    point1,
                    body2: c2,
                    point2,
                } => {
                    let rw1 = rots[c1].mul_vec3(point1);
                    let rw2 = rots[c2].mul_vec3(point2);
                    let w1 = positions[c1] + rw1;
                    let w2 = positions[c2] + rw2;
                    let diff = w1 - w2;
                    rj.r[row] = diff.x;
                    rj.r[row + 1] = diff.y;
                    rj.r[row + 2] = diff.z;
                    for k in 0..3 {
                        // Translation
                        for jt in 0..3 {
                            rj.j[row + k][c1 * 6 + jt] = if jt == k { 1.0 } else { 0.0 };
                            rj.j[row + k][c2 * 6 + jt] = if jt == k { -1.0 } else { 0.0 };
                        }
                        // Rotation columns: the solver's increment is a
                        // BODY-frame omega (quaternion right-multiply), so
                        // the perturbed world offset is
                        // rw + delta*(R*e_j)xrw — the axis must be the body
                        // axis carried to world, not the world axis (the
                        // analytic-vs-FD gate catches this).
                        for jw in 0..3 {
                            let axis1 = rots[c1].mul_vec3(e_axis(jw));
                            let axis2 = rots[c2].mul_vec3(e_axis(jw));
                            let om1 = axis1.cross(rw1);
                            let om2 = axis2.cross(rw2);
                            rj.j[row + k][c1 * 6 + 3 + jw] = om1[k];
                            rj.j[row + k][c2 * 6 + 3 + jw] = -om2[k];
                        }
                    }
                    row += 3;
                }
                ConstraintKind::Distance {
                    body1: c1,
                    point1,
                    body2: c2,
                    point2,
                    distance,
                } => {
                    let rw1 = rots[c1].mul_vec3(point1);
                    let rw2 = rots[c2].mul_vec3(point2);
                    let w1 = positions[c1] + rw1;
                    let w2 = positions[c2] + rw2;
                    let delta = w1 - w2;
                    let len = delta.length();
                    let u = if len > 1e-14 { delta / len } else { DVec3::X };
                    rj.r[row] = len - distance;
                    for jt in 0..3 {
                        rj.j[row][c1 * 6 + jt] = u[jt];
                        rj.j[row][c2 * 6 + jt] = -u[jt];
                    }
                    for jw in 0..3 {
                        // Body-frame axis carried to world (see the
                        // Coincident rotation-column note).
                        let axis1 = rots[c1].mul_vec3(e_axis(jw));
                        let axis2 = rots[c2].mul_vec3(e_axis(jw));
                        let om1 = axis1.cross(rw1);
                        let om2 = axis2.cross(rw2);
                        rj.j[row][c1 * 6 + 3 + jw] = u.dot(om1);
                        rj.j[row][c2 * 6 + 3 + jw] = -u.dot(om2);
                    }
                    row += 1;
                }
                ConstraintKind::Parallel {
                    body1: c1,
                    dir1,
                    body2: c2,
                    dir2,
                } => {
                    // Two rows for two freedoms: see
                    // `fill_direction_tangent_rows`.
                    let basis = world_tangent_basis(rots[c1], dir1);
                    let d2w = world_unit_dir(rots[c2], dir2);
                    Self::fill_direction_tangent_rows((c1, basis), (c2, d2w), rots, row, rj);
                    row += 2;
                }
                ConstraintKind::Angle {
                    body1: c1,
                    dir1,
                    body2: c2,
                    dir2,
                    angle,
                } => {
                    let d1w = world_unit_dir(rots[c1], dir1);
                    let d2w = world_unit_dir(rots[c2], dir2);
                    rj.r[row] = d1w.dot(d2w) - angle.cos();
                    Self::fill_direction_dot_jacobian((c1, d1w), (c2, d2w), rots, row, &mut rj.j);
                    row += 1;
                }
                ConstraintKind::Perpendicular {
                    body1: c1,
                    dir1,
                    body2: c2,
                    dir2,
                } => {
                    let d1w = world_unit_dir(rots[c1], dir1);
                    let d2w = world_unit_dir(rots[c2], dir2);
                    rj.r[row] = d1w.dot(d2w);
                    Self::fill_direction_dot_jacobian((c1, d1w), (c2, d2w), rots, row, &mut rj.j);
                    row += 1;
                }
                ConstraintKind::Concentric {
                    body1: c1,
                    point1,
                    dir1,
                    body2: c2,
                    point2,
                    dir2,
                } => {
                    // Four rows for four freedoms. Rows 0-1: axis 2 is
                    // parallel to axis 1, as in `Parallel`. Rows 2-3: the
                    // offset between the two axis points has no component
                    // across axis 1, r_k = e_k · (p2w − p1w), over the same
                    // basis e_1, e_2 of the plane perpendicular to axis 1.
                    let basis = world_tangent_basis(rots[c1], dir1);
                    let d2w = world_unit_dir(rots[c2], dir2);
                    Self::fill_direction_tangent_rows((c1, basis), (c2, d2w), rots, row, rj);

                    // Rotation columns differentiate rw = R·local (rotation
                    // only — translation enters separately), NOT the
                    // translated p_w.
                    let rw1 = rots[c1].mul_vec3(point1);
                    let rw2 = rots[c2].mul_vec3(point2);
                    let off = (positions[c2] + rw2) - (positions[c1] + rw1);
                    for (k, e) in basis.into_iter().enumerate() {
                        let r = row + 2 + k;
                        rj.r[r] = e.dot(off);
                        for jt in 0..3 {
                            rj.j[r][c1 * 6 + jt] += -e[jt];
                            rj.j[r][c2 * 6 + jt] += e[jt];
                        }
                        for jw in 0..3 {
                            // Body 1 turns the basis vector and moves p1:
                            // δr = δe·off − e·δp1.
                            let axis1 = rots[c1].mul_vec3(e_axis(jw));
                            rj.j[r][c1 * 6 + 3 + jw] +=
                                axis1.cross(e).dot(off) - e.dot(axis1.cross(rw1));
                            // Body 2 moves p2: δr = e·δp2.
                            let axis2 = rots[c2].mul_vec3(e_axis(jw));
                            rj.j[r][c2 * 6 + 3 + jw] += e.dot(axis2.cross(rw2));
                        }
                    }
                    row += 4;
                }
                ConstraintKind::PointOnPlane {
                    body1: c1,
                    plane_point,
                    plane_normal,
                    body2: c2,
                    point,
                } => {
                    // The plane rides body 1; the point is on body 2. One
                    // row: the signed perpendicular distance (w2 - w1) . n1w.
                    let rw1 = rots[c1].mul_vec3(plane_point);
                    let w1 = positions[c1] + rw1;
                    let n1w = world_unit_dir(rots[c1], plane_normal);
                    let rw2 = rots[c2].mul_vec3(point);
                    let w2 = positions[c2] + rw2;
                    rj.r[row] = (w2 - w1).dot(n1w);
                    for jt in 0..3 {
                        rj.j[row][c1 * 6 + jt] = -n1w[jt];
                        rj.j[row][c2 * 6 + jt] = n1w[jt];
                    }
                    // Rotation columns: dr = -dw1.n1w + diff.dn1w with
                    // dw1 = axis1 x rw1, dn1w = axis1 x n1w, which reduces
                    // to axis1 . (n1w x off), off = w2 - t1.
                    let off = w2 - positions[c1];
                    for jw in 0..3 {
                        let axis1 = rots[c1].mul_vec3(e_axis(jw));
                        rj.j[row][c1 * 6 + 3 + jw] = axis1.dot(n1w.cross(off));
                        let axis2 = rots[c2].mul_vec3(e_axis(jw));
                        rj.j[row][c2 * 6 + 3 + jw] = axis2.dot(rw2.cross(n1w));
                    }
                    row += 1;
                }
                ConstraintKind::Lock {
                    body1: c1,
                    body2: c2,
                    translation: rel_t,
                    rotation: rel_q,
                } => {
                    // Relative-pose lock. Translation rows are analytic;
                    // the rotation rows reuse the Fixed-orientation FD
                    // helper (the rotation-vector Jacobian is not worth an
                    // analytic form here).
                    let r1i = rots[c1].inverse();
                    let t_rel = r1i.mul_vec3(positions[c2] - positions[c1]);
                    rj.r[row] = t_rel.x - rel_t.x;
                    rj.r[row + 1] = t_rel.y - rel_t.y;
                    rj.r[row + 2] = t_rel.z - rel_t.z;
                    let q_rel = r1i * rots[c2];
                    let res = (q_rel * rel_q.inverse()).normalize().to_scaled_axis();
                    rj.r[row + 3] = res.x;
                    rj.r[row + 4] = res.y;
                    rj.r[row + 5] = res.z;
                    for k in 0..3 {
                        for jt in 0..3 {
                            let col = r1i.mul_vec3(e_axis(jt));
                            rj.j[row + k][c1 * 6 + jt] = -col[k];
                            rj.j[row + k][c2 * 6 + jt] = col[k];
                        }
                    }
                    for jw in 0..3 {
                        let col = t_rel.cross(e_axis(jw));
                        for k in 0..3 {
                            rj.j[row + k][c1 * 6 + 3 + jw] = col[k];
                        }
                    }
                    Self::fill_vec3_residual_fd(
                        |_, rp| {
                            let q_err = (rp[c1].inverse() * rp[c2] * rel_q.inverse()).normalize();
                            q_err.to_scaled_axis()
                        },
                        positions,
                        rots,
                        row + 3,
                        cols,
                        &mut rj.j,
                    );
                    row += 6;
                }
            }
        }

        debug_assert_eq!(row, m);
    }

    fn fill_vec3_residual_fd<F>(
        residual_fn: F,
        positions: &[DVec3],
        rots: &[DQuat],
        row: usize,
        cols: usize,
        j: &mut [Vec<f64>],
    ) where
        F: Fn(&[DVec3], &[DQuat]) -> DVec3 + Copy,
    {
        // One working copy, perturbed and restored a single body at a time
        // (each column moves exactly one body).
        let mut pp = positions.to_vec();
        let mut rp = rots.to_vec();
        for col in 0..cols {
            let body = col / 6;

            apply_parameter_delta(&mut pp, &mut rp, col, FD_EPS);
            let vp = residual_fn(&pp, &rp);
            pp[body] = positions[body];
            rp[body] = rots[body];

            apply_parameter_delta(&mut pp, &mut rp, col, -FD_EPS);
            let vm = residual_fn(&pp, &rp);
            pp[body] = positions[body];
            rp[body] = rots[body];

            for k in 0..3 {
                j[row + k][col] = (vp[k] - vm[k]) / (2.0 * FD_EPS);
            }
        }
    }

    /// Two rows saying that the direction `d2w` (on body `c2`) is parallel
    /// to a direction on body `c1`, given by a `basis` e_1, e_2 of the plane
    /// perpendicular to it: `r_k = e_k · d2`.
    ///
    /// Both rows vanish exactly when the directions are parallel, either
    /// way round, and their gradients are independent there, so the
    /// condition costs two equations for the two freedoms it removes. (The
    /// cross product `d1 × d2` says the same in three rows of rank two; the
    /// spare row then reads as a redundancy in the diagnosis.) The two forms
    /// have the same sum of squares, `sin²` of the angle between the
    /// directions.
    ///
    /// δr_k = δe_k·d2 + e_k·δd2, where a vector carried by body i turns by
    /// (R_i·e_j)×v for the body-frame omega column j. The translation
    /// columns are zero. Entries are added, so that a mate between a body
    /// and itself sums its two contributions.
    fn fill_direction_tangent_rows(
        (c1, basis): (usize, [DVec3; 2]),
        (c2, d2w): (usize, DVec3),
        rots: &[DQuat],
        row: usize,
        rj: &mut ResidualJacobian,
    ) {
        for (k, e) in basis.into_iter().enumerate() {
            rj.r[row + k] = e.dot(d2w);
            for jw in 0..3 {
                let axis1 = rots[c1].mul_vec3(e_axis(jw));
                rj.j[row + k][c1 * 6 + 3 + jw] += axis1.cross(e).dot(d2w);
                let axis2 = rots[c2].mul_vec3(e_axis(jw));
                rj.j[row + k][c2 * 6 + 3 + jw] += e.dot(axis2.cross(d2w));
            }
        }
    }

    /// Analytic Jacobian of the row `r = d1·d2 (− const)` of two world
    /// directions: δr = δd1·d2 + d1·δd2. Serves `Angle` and
    /// `Perpendicular`.
    fn fill_direction_dot_jacobian(
        (c1, d1w): (usize, DVec3),
        (c2, d2w): (usize, DVec3),
        rots: &[DQuat],
        row: usize,
        j: &mut [Vec<f64>],
    ) {
        for jw in 0..3 {
            let axis1 = rots[c1].mul_vec3(e_axis(jw));
            let a1 = axis1.cross(d1w);
            j[row][c1 * 6 + 3 + jw] = a1.dot(d2w);
            let axis2 = rots[c2].mul_vec3(e_axis(jw));
            let a2 = axis2.cross(d2w);
            j[row][c2 * 6 + 3 + jw] = d1w.dot(a2);
        }
    }
}

#[inline]
fn e_axis(j: usize) -> DVec3 {
    match j {
        0 => DVec3::X,
        1 => DVec3::Y,
        _ => DVec3::Z,
    }
}

#[cfg(test)]
fn l2_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn apply_parameter_delta(pos: &mut [DVec3], rot: &mut [DQuat], col: usize, delta: f64) {
    let comp = col / 6;
    let dim = col % 6;
    match dim {
        0 => pos[comp].x += delta,
        1 => pos[comp].y += delta,
        2 => pos[comp].z += delta,
        3 => {
            rot[comp] = (rot[comp] * DQuat::from_scaled_axis(DVec3::X * delta)).normalize();
        }
        4 => {
            rot[comp] = (rot[comp] * DQuat::from_scaled_axis(DVec3::Y * delta)).normalize();
        }
        5 => {
            rot[comp] = (rot[comp] * DQuat::from_scaled_axis(DVec3::Z * delta)).normalize();
        }
        _ => {}
    }
}

fn quats_to_euler_xyz(quats: &[DQuat]) -> Vec<DVec3> {
    quats
        .iter()
        .map(|q| {
            let (xr, yr, zr) = q.to_euler(EulerRot::XYZ);
            DVec3::new(xr, yr, zr)
        })
        .collect()
}

fn world_unit_dir(rot: DQuat, local: DVec3) -> DVec3 {
    let v = rot.mul_vec3(local);
    let len = v.length();
    if len < 1e-14 {
        DVec3::Z
    } else {
        v / len
    }
}

/// A basis of the plane perpendicular to a body's local direction, carried
/// to world: two orthonormal vectors that turn with the body.
///
/// The basis is chosen in the body frame, from the coordinate axis the
/// direction is least aligned with. It therefore depends only on the
/// constraint's own data, not on the pose, and the world vectors vary
/// smoothly as the body moves.
fn world_tangent_basis(rot: DQuat, local: DVec3) -> [DVec3; 2] {
    let unit = local.try_normalize().unwrap_or(DVec3::Z);
    let (ax, ay, az) = (unit.x.abs(), unit.y.abs(), unit.z.abs());
    let helper = if ax <= ay && ax <= az {
        DVec3::X
    } else if ay <= az {
        DVec3::Y
    } else {
        DVec3::Z
    };
    let e1 = unit.cross(helper).normalize();
    let e2 = unit.cross(e1);
    [rot.mul_vec3(e1), rot.mul_vec3(e2)]
}

/// Rotation vector residual locking current orientation `q` to target `q_tgt`.
fn fixed_orientation_residual(q: DQuat, q_tgt: DQuat) -> DVec3 {
    let q_err = (q * q_tgt.inverse()).normalize();
    q_err.to_scaled_axis()
}

impl Default for ConstraintSolver {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rigid3d::constraint::Constraint;

    #[test]
    fn test_solve_fixed_constraint() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::new(1.0, 2.0, 3.0)));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        assert!(result.success);
        assert!((result.positions[0].x - 1.0).abs() < 0.01);
        assert!((result.positions[0].y - 2.0).abs() < 0.01);
        assert!((result.positions[0].z - 3.0).abs() < 0.01);
    }

    #[test]
    fn test_solve_fixed_with_orientation() {
        let tgt = DQuat::from_rotation_z(std::f64::consts::FRAC_PI_4);
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed_with_orientation(
            0,
            0,
            DVec3::new(0.0, 1.0, 2.0),
            tgt,
        ));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        assert!(result.success, "error={}", result.error);
        assert!((result.positions[0] - DVec3::new(0.0, 1.0, 2.0)).length() < 1e-5);
        let (ex, ey, ez) = tgt.to_euler(EulerRot::XYZ);
        assert!((result.orientations[0].x - ex).abs() < 1e-3);
        assert!((result.orientations[0].y - ey).abs() < 1e-3);
        assert!((result.orientations[0].z - ez).abs() < 1e-3);
    }

    #[test]
    fn test_solve_coincident_constraint() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
        system.add_constraint(Constraint::coincident(
            1,
            0,
            DVec3::ZERO,
            1,
            DVec3::new(1.0, 0.0, 0.0),
        ));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        assert!(result.success);
    }

    #[test]
    fn polish_recovers_fixed_pose_bitwise_from_warm_start() {
        // The M4-style round trip: warm start slightly off the authored
        // pose. The loop converges under the tolerance, and the
        // convergence polish must pin the exactly-representable fixed
        // pose BITWISE — not merely within tolerance.
        let target = DVec3::new(1.0, 2.0, 3.0);
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, target));

        let solver = ConstraintSolver::new();
        let start = target + DVec3::new(1.0e-7, -2.0e-7, 5.0e-8);
        let raw_error = l2_norm(
            &ConstraintSolver::build_residual_jacobian(&system, &[start], &[DQuat::IDENTITY]).r,
        );
        assert!(
            raw_error < solver.options.tolerance,
            "warm start must sit inside the convergence tolerance for the polish path: {raw_error:.3e}",
        );

        let result = solver.solve_from(&system, vec![start], vec![DQuat::IDENTITY]);

        assert!(result.success, "error={}", result.error);
        assert!(result.error < raw_error, "polish must reduce the error");
        assert_eq!(
            result.positions[0], target,
            "polish must land bitwise on the fixed pose"
        );
    }

    #[test]
    fn test_solve_distance_constraint() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
        system.add_constraint(Constraint::distance(
            1,
            0,
            DVec3::ZERO,
            1,
            DVec3::new(1.0, 0.0, 0.0),
            5.0,
        ));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        let p0 = result.positions[0];
        let p1 = result.positions[1];
        let world_p0 = p0 + DVec3::ZERO;
        let world_p1 = p1 + DVec3::new(1.0, 0.0, 0.0);
        let dist = (world_p1 - world_p0).length();
        assert!(
            (dist - 5.0).abs() < 0.1,
            "Distance is {} but expected 5.0",
            dist
        );
    }

    #[test]
    fn test_three_component_chain() {
        let mut system = ConstraintSystem::with_components(3);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
        system.add_constraint(Constraint::coincident(
            1,
            0,
            DVec3::new(1.0, 0.0, 0.0),
            1,
            DVec3::ZERO,
        ));
        system.add_constraint(Constraint::coincident(
            2,
            1,
            DVec3::new(1.0, 0.0, 0.0),
            2,
            DVec3::ZERO,
        ));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        assert!(result.success);
    }

    #[test]
    fn test_four_component_square() {
        let mut system = ConstraintSystem::with_components(4);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));

        system.add_constraint(Constraint::distance(
            1,
            0,
            DVec3::ZERO,
            1,
            DVec3::ZERO,
            10.0,
        ));

        system.add_constraint(Constraint::distance(
            2,
            1,
            DVec3::ZERO,
            2,
            DVec3::ZERO,
            10.0,
        ));

        system.add_constraint(Constraint::distance(
            3,
            0,
            DVec3::ZERO,
            3,
            DVec3::ZERO,
            10.0,
        ));
        system.add_constraint(Constraint::distance(
            4,
            2,
            DVec3::ZERO,
            3,
            DVec3::ZERO,
            10.0,
        ));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        assert!(result.iterations > 0);
    }

    #[test]
    fn test_parallel_constraints() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
        system.add_constraint(Constraint::parallel(
            1,
            0,
            DVec3::new(0.0, 0.0, 1.0),
            1,
            DVec3::new(0.0, 0.0, 1.0),
        ));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        assert!(
            result.success,
            "parallel Z mates with identity rotation should vanish"
        );
        assert!(result.error < 1e-5);
    }

    #[test]
    fn test_perpendicular_constraint() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
        system.add_constraint(Constraint::perpendicular(
            1,
            0,
            DVec3::new(0.0, 0.0, 1.0),
            1,
            DVec3::new(1.0, 0.0, 0.0),
        ));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);
        assert!(result.success);
        assert!(result.error < 1e-5);
    }

    #[test]
    fn test_angle_constraint_aligned() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
        system.add_constraint(Constraint::angle(
            1,
            0,
            DVec3::new(0.0, 0.0, 1.0),
            1,
            DVec3::new(0.0, 0.0, 1.0),
            0.0,
        ));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);
        assert!(result.success);
        assert!(result.error < 1e-5);
    }

    #[test]
    fn test_overconstrained_system() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
        system.add_constraint(Constraint::fixed(1, 0, DVec3::new(1.0, 0.0, 0.0)));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        assert!((result.positions[0].x - 0.5).abs() < 0.01);
        assert!(!result.success);
    }

    #[test]
    fn test_many_components() {
        let n = 10;
        let mut system = ConstraintSystem::with_components(n);

        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));

        for i in 1..n {
            system.add_constraint(Constraint::distance(
                i,
                i - 1,
                DVec3::ZERO,
                i,
                DVec3::ZERO,
                1.0,
            ));
        }

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        assert!(result.iterations <= 100);
    }

    #[test]
    fn test_tree_structure() {
        let mut system = ConstraintSystem::with_components(5);

        system.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));

        system.add_constraint(Constraint::distance(1, 0, DVec3::ZERO, 1, DVec3::ZERO, 5.0));
        system.add_constraint(Constraint::distance(2, 1, DVec3::ZERO, 2, DVec3::ZERO, 3.0));

        system.add_constraint(Constraint::distance(3, 0, DVec3::ZERO, 3, DVec3::ZERO, 7.0));
        system.add_constraint(Constraint::distance(4, 3, DVec3::ZERO, 4, DVec3::ZERO, 4.0));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&system);

        assert!(result.iterations > 0);
    }

    #[test]
    fn test_solver_options() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 0, DVec3::new(5.0, 5.0, 5.0)));
        system.add_constraint(Constraint::coincident(1, 0, DVec3::ZERO, 1, DVec3::ZERO));

        let options = SolveOptions {
            max_iterations: 50,
            tolerance: 1e-8,
            damping: 0.8,
            ..SolveOptions::default()
        };

        let solver = ConstraintSolver::with_options(options);
        let result = solver.solve(&system);

        assert!(result.success);
        assert!(result.iterations <= 50);
    }
}

#[cfg(test)]
mod gate_tests {
    use super::*;
    use crate::rigid3d::constraint::ConstraintType;

    /// The 3D analytic-vs-FD gate (solver roadmap M6.1): at NON-identity
    /// body orientations, every analytic Jacobian column must match the
    /// central finite difference of the residual map under the solver's own
    /// body-frame quaternion parameterization.
    /// The two-row form of "parallel" and the four-row form of
    /// "concentric" measure the same misalignment as the cross products
    /// they replace: |d1 × d2| for the directions, and |(p2 − p1) × d1| for
    /// the offset between the axes.
    #[test]
    fn tangent_rows_have_the_norm_of_the_cross_products() {
        let positions = [DVec3::new(0.4, -1.0, 0.3), DVec3::new(1.2, -0.4, 0.8)];
        let rots = [
            DQuat::from_euler(EulerRot::XYZ, -0.3, 0.6, 1.1),
            DQuat::from_euler(EulerRot::XYZ, 0.5, -0.9, 0.3),
        ];
        let (point1, dir1) = (DVec3::new(0.4, -0.1, 0.2), DVec3::new(0.3, -1.2, 0.5));
        let (point2, dir2) = (DVec3::new(-0.2, 0.3, 0.1), DVec3::new(2.0, 0.1, -0.7));
        let d1w = world_unit_dir(rots[0], dir1);
        let d2w = world_unit_dir(rots[1], dir2);
        let off =
            (positions[1] + rots[1].mul_vec3(point2)) - (positions[0] + rots[0].mul_vec3(point1));

        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::parallel(0, 0, dir1, 1, dir2));
        system.add_constraint(Constraint::concentric(1, 0, point1, dir1, 1, point2, dir2));
        let r = ConstraintSolver::build_residual_jacobian(&system, &positions, &rots).r;
        assert_eq!(r.len(), 2 + 4);

        let directions = d1w.cross(d2w).length();
        let offset = off.cross(d1w).length();
        assert!(directions > 0.1 && offset > 0.1, "a generic pose");
        assert!((l2_norm(&r[0..2]) - directions).abs() < 1e-14);
        assert!((l2_norm(&r[2..4]) - directions).abs() < 1e-14);
        assert!((l2_norm(&r[4..6]) - offset).abs() < 1e-14);
    }

    #[test]
    fn analytic_jacobian_matches_fd_at_rotated_pose() {
        let mut system = ConstraintSystem::with_components(3);
        system.add_constraint(Constraint::fixed_with_orientation(
            0,
            0,
            DVec3::ZERO,
            DQuat::from_rotation_y(0.7),
        ));
        system.add_constraint(Constraint::coincident(
            1,
            0,
            DVec3::new(0.5, -0.2, 0.3),
            1,
            DVec3::new(1.0, 0.5, -0.5),
        ));
        system.add_constraint(Constraint::distance(
            2,
            0,
            DVec3::new(0.2, 0.4, 0.1),
            1,
            DVec3::new(-0.3, 0.2, 0.6),
            2.5,
        ));
        // Direction-type mates: analytic Jacobians under test. The refs
        // are UNIT local directions (the analytic form omits the
        // normalization derivative by authoring contract).
        system.add_constraint(Constraint::parallel(
            3,
            1,
            DVec3::new(0.0, 1.0, 0.0),
            2,
            DVec3::new(0.6, 0.8, 0.0),
        ));
        system.add_constraint(Constraint::angle(
            4,
            0,
            DVec3::new(0.0, 0.0, 1.0),
            1,
            DVec3::new(1.0, 0.0, 0.0),
            1.1,
        ));
        system.add_constraint(Constraint::perpendicular(
            5,
            1,
            DVec3::new(0.0, 0.0, 1.0),
            2,
            DVec3::new(0.0, 1.0, 0.0),
        ));
        system.add_constraint(Constraint::concentric(
            6,
            1,
            DVec3::new(0.4, -0.1, 0.2),
            DVec3::new(0.0, 0.0, 1.0),
            2,
            DVec3::new(-0.2, 0.3, 0.1),
            DVec3::new(0.0, 1.0, 0.0),
        ));
        // Point-on-plane: the plane rides body 1, the point is on body 2.
        system.add_constraint(Constraint::point_on_plane(
            7,
            1,
            DVec3::new(0.3, 0.2, -0.1),
            DVec3::new(0.0, 0.0, 1.0),
            2,
            DVec3::new(0.1, -0.4, 0.6),
        ));
        // Relative-pose lock between bodies 1 and 2.
        system.add_constraint(Constraint::lock(
            8,
            1,
            2,
            DVec3::new(0.7, -0.2, 0.4),
            DQuat::from_euler(EulerRot::XYZ, 0.2, 0.3, -0.1),
        ));

        // The same two mates again, on directions that are neither unit
        // length nor along a coordinate axis: the rows are built on a basis
        // of the plane across the first direction, and its choice must not
        // show in the derivative.
        system.add_constraint(Constraint::parallel(
            9,
            2,
            DVec3::new(0.3, -1.2, 0.5),
            1,
            DVec3::new(2.0, 0.1, -0.7),
        ));
        system.add_constraint(Constraint::concentric(
            10,
            2,
            DVec3::new(0.6, 0.1, -0.3),
            DVec3::new(-0.4, 0.2, 1.5),
            0,
            DVec3::new(0.2, -0.5, 0.4),
            DVec3::new(1.0, 1.0, 0.2),
        ));
        // Both directions on one body: the residual cannot change, and the
        // two contributions to each column have to cancel.
        system.add_constraint(Constraint::parallel(
            11,
            1,
            DVec3::new(0.2, 0.9, -0.4),
            1,
            DVec3::new(1.0, 0.3, 0.5),
        ));

        // Body 1 sits at a non-identity orientation — the regime where a
        // world-axis vs body-axis basis mixup in the analytic columns shows.
        let positions = vec![
            DVec3::ZERO,
            DVec3::new(1.2, -0.4, 0.8),
            DVec3::new(0.3, 2.0, -1.0),
        ];
        let rots = vec![
            DQuat::from_rotation_y(0.7),
            DQuat::from_euler(EulerRot::XYZ, 0.5, -0.9, 0.3),
            DQuat::from_rotation_z(-1.1),
        ];

        let rj = ConstraintSolver::build_residual_jacobian(&system, &positions, &rots);

        // FD reference helper: a full rebuild per perturbed pose, residuals
        // extracted from the joint residual+Jacobian pass.
        let residuals_at = |pp: &[DVec3], rp: &[DQuat]| -> Vec<f64> {
            ConstraintSolver::build_residual_jacobian(&system, pp, rp).r
        };

        // Full finite-difference reference over the same parameterization.
        let m = rj.r.len();
        let cols = positions.len() * 6;
        let mut jfd = vec![vec![0.0_f64; cols]; m];
        for col in 0..cols {
            let mut pp = positions.to_vec();
            let mut rp = rots.to_vec();
            apply_parameter_delta(&mut pp, &mut rp, col, FD_EPS);
            let mut pm = positions.to_vec();
            let mut rm = rots.to_vec();
            apply_parameter_delta(&mut pm, &mut rm, col, -FD_EPS);
            let fp = residuals_at(&pp, &rp);
            let fm = residuals_at(&pm, &rm);
            for row in 0..m {
                jfd[row][col] = (fp[row] - fm[row]) / (2.0 * FD_EPS);
            }
        }

        let mut bad = 0usize;
        for row in 0..m {
            for col in 0..cols {
                let (va, vf) = (rj.j[row][col], jfd[row][col]);
                if (va - vf).abs() > 1.0e-5 * (1.0 + vf.abs()) {
                    println!("DRIFT row {row} col {col}: analytic {va} vs fd {vf}");
                    bad += 1;
                }
            }
        }
        assert_eq!(bad, 0, "analytic Jacobian drifts at a rotated pose");
        let _ = ConstraintType::Fixed;
    }
}

#[cfg(test)]
mod assembly_acceptance {
    use super::*;
    use crate::rigid3d::constraint::Constraint;

    /// Ground + one body pinned by THREE coincident pairs (9 consistent
    /// equations over 6 DOF): a fully determined assembly whose solution
    /// is the authored pose. M6.1 acceptance: converges; pose drifts
    /// re-solve back within 1e-9; re-solve from the same start is bitwise
    /// deterministic.
    fn build() -> ConstraintSystem {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed_with_orientation(
            0,
            0,
            DVec3::ZERO,
            DQuat::IDENTITY,
        ));
        // Body 1 authored at t=(1,0,1), R=I: its locals (0,0,0)/(2,0,0)/
        // (0,3,1) land on the ground targets below.
        system.add_constraint(Constraint::coincident(
            1,
            0,
            DVec3::new(1.0, 0.0, 1.0),
            1,
            DVec3::new(0.0, 0.0, 0.0),
        ));
        system.add_constraint(Constraint::coincident(
            2,
            0,
            DVec3::new(3.0, 0.0, 1.0),
            1,
            DVec3::new(2.0, 0.0, 0.0),
        ));
        system.add_constraint(Constraint::coincident(
            3,
            0,
            DVec3::new(1.0, 3.0, 2.0),
            1,
            DVec3::new(0.0, 3.0, 1.0),
        ));
        system
    }

    const TOL: f64 = 1e-9;

    #[test]
    fn two_body_assembly_solves_and_round_trips() {
        let system = build();
        let solver = ConstraintSolver::new();
        let r0 = solver.solve(&system);
        assert!(r0.success, "error {}", r0.error);
        let t0 = r0.positions[1];
        assert!(
            (t0 - DVec3::new(1.0, 0.0, 1.0)).length() < 1e-7,
            "t0 {t0:?}"
        );

        // Pose excursion: shove the free body well off its solution,
        // re-solve — it must return within the contract tolerance.
        let mut pos = r0.positions.clone();
        let mut rot = r0.rotations.clone();
        pos[1] += DVec3::new(0.4, -0.3, 0.6);
        rot[1] = (rot[1] * DQuat::from_scaled_axis(DVec3::new(0.05, -0.08, 0.03))).normalize();
        let r1 = solver.solve_from(&system, pos, rot);
        assert!(r1.success, "re-solve error {}", r1.error);
        assert!(
            (r1.positions[1] - t0).length() < TOL,
            "round-trip translation drift {:?}",
            r1.positions[1] - t0
        );
        let qerr = (r1.rotations[1] * r0.rotations[1].inverse()).normalize();
        assert!(
            qerr.to_scaled_axis().length() < TOL,
            "round-trip rotation drift {:?}",
            qerr.to_scaled_axis()
        );

        // Determinism: the same excursion re-solved is bitwise identical.
        let mut pos2 = r0.positions.clone();
        let mut rot2 = r0.rotations.clone();
        pos2[1] += DVec3::new(0.4, -0.3, 0.6);
        rot2[1] = (rot2[1] * DQuat::from_scaled_axis(DVec3::new(0.05, -0.08, 0.03))).normalize();
        let r2 = solver.solve_from(&system, pos2, rot2);
        assert_eq!(
            r1.positions[1], r2.positions[1],
            "re-solve must be deterministic"
        );
    }

    /// Localized assembly solve (M3, assembly edition): a lone fixed body
    /// forms its own component; re-solving only the main cluster leaves
    /// the lone body bitwise untouched, and the piecewise solve lands on
    /// the same poses as the all-at-once solve.
    #[test]
    fn per_component_solve_matches_and_skips_untouched() {
        // build() + a lone fixed body far away.
        let mut system = build();
        system.num_components = 3;
        system.add_constraint(Constraint::fixed_with_orientation(
            4,
            2,
            DVec3::new(50.0, 0.0, 0.0),
            DQuat::from_rotation_z(0.3),
        ));

        let solver = ConstraintSolver::new();
        let comps = ConstraintSolver::components(&system);
        assert_eq!(comps.len(), 2, "main cluster + lone body: {comps:?}");

        let all = system.clone();
        let r_all = solver.solve(&all);
        assert!(r_all.success);

        // Piecewise: solve component 0, then component 1, chaining the
        // warm start through.
        let (seed_pos, seed_rot) = ConstraintSolver::seed_fixed(&system);
        let r0 = solver.solve_components_from(
            &system,
            &comps,
            Some(&[0]),
            seed_pos.clone(),
            seed_rot.clone(),
        );
        assert!(r0.success);
        let r1 =
            solver.solve_components_from(&system, &comps, Some(&[1]), r0.positions, r0.rotations);
        assert!(r1.success);
        for (a, b) in r_all.positions.iter().zip(r1.positions.iter()) {
            assert!(
                (a - b).length() < 1e-9,
                "per-component divergence {a:?} vs {b:?}"
            );
        }

        // Untouched-cluster skip: excursions to the main cluster leave the
        // lone body bitwise at its seeded pose.
        let lone_before = seed_pos[2];
        let mut pos = r_all.positions.clone();
        let mut rot = r_all.rotations.clone();
        pos[1] += DVec3::new(0.7, 0.2, -0.4);
        rot[1] = (rot[1] * DQuat::from_scaled_axis(DVec3::new(0.1, 0.0, -0.05))).normalize();
        let r = solver.solve_components_from(&system, &comps, Some(&[0]), pos, rot);
        assert!(r.success);
        assert_eq!(r.positions[2], lone_before, "untouched cluster moved");
    }
}
