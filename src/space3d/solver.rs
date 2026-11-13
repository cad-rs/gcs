//! Newton-Raphson GCS solver for 3D geometric constraints.
//!
//! A thin family adapter, exactly like the 2D solver: the 3D constraint set
//! contributes the residuals and the analytic Jacobian rows; the
//! decomposition and the per-component solve are the shared pipeline of
//! the crate-internal `system` module over the Newton core of [`crate::numerics`]. What is
//! specific to this family is the plane representation, which is pinned to
//! a unit normal around every solve.
//!
//! The solve runs PER CONNECTED COMPONENT of the constraint graph
//! ([`crate::space3d::graph`]). A full solve is therefore the composition
//! of independent component solves, and a localized (incremental) re-solve
//! iterates only the components an edit landed in — untouched components
//! keep their parameters bitwise.

use crate::edit::ConstraintState;
use crate::numerics::{self, CoreOptions, SolveOptions, SolveResult, SolveStatus};
use crate::space3d::constraint::SpaceConstraint;
use crate::space3d::entity::{SpaceEntity, SpaceEntityKind};
use crate::system::{self, Equations, SparseRow};

pub use crate::numerics::{FD_H, LAMBDA, MAX_ITER, PIVOT_TOL, RESIDUAL_TOL, SOLVER_NORM_TOL};

/// The 3D constraint list as the pipeline's [`Equations`].
pub(crate) struct SpaceEquations<'a> {
    pub(crate) entities: &'a [SpaceEntity],
    constraints: &'a [SpaceConstraint],
    /// The state of each constraint slot; a slot beyond the end of this
    /// list, or no list at all, counts as active.
    states: Option<&'a [ConstraintState]>,
}

impl<'a> SpaceEquations<'a> {
    /// Every constraint in the list is enforced.
    pub(crate) fn new(entities: &'a [SpaceEntity], constraints: &'a [SpaceConstraint]) -> Self {
        Self {
            entities,
            constraints,
            states: None,
        }
    }

    /// Only the constraints whose slot is [`ConstraintState::Active`] are
    /// enforced. A suppressed or removed one contributes no equation and
    /// keeps its index.
    pub(crate) fn with_states(
        entities: &'a [SpaceEntity],
        constraints: &'a [SpaceConstraint],
        states: &'a [ConstraintState],
    ) -> Self {
        Self {
            entities,
            constraints,
            states: Some(states),
        }
    }

    fn is_active(&self, index: usize) -> bool {
        self.states
            .and_then(|states| states.get(index))
            .is_none_or(|state| *state == ConstraintState::Active)
    }
}

impl Equations for SpaceEquations<'_> {
    fn constraint_count(&self) -> usize {
        self.constraints.len()
    }

    fn equation_count(&self, index: usize) -> usize {
        if self.is_active(index) {
            self.constraints[index].equation_count()
        } else {
            0
        }
    }

    fn residuals(&self, index: usize, params: &[f64], out: &mut [f64]) {
        if self.is_active(index) {
            self.constraints[index].residuals(params, self.entities, out);
        }
    }

    fn jacobian(&self, index: usize, params: &[f64], out: &mut [SparseRow]) {
        if self.is_active(index) {
            self.constraints[index].jacobian(params, self.entities, out);
        }
    }
}

/// Evaluate the stacked residual vector of `constraints` at `params`.
#[cfg(test)]
pub(crate) fn eval_residuals(
    params: &[f64],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
    n_eq: usize,
) -> Vec<f64> {
    let eqs = SpaceEquations::new(entities, constraints);
    let f = system::eval_residuals(&eqs, &system::all_indices(&eqs), params);
    debug_assert_eq!(f.len(), n_eq);
    f
}

/// Assemble the analytic Jacobian (m equations × the columns of `col_of`).
/// `out` arrives zeroed.
#[cfg(test)]
pub(crate) fn analytic_jacobian(
    params: &[f64],
    col_of: &std::collections::HashMap<usize, usize>,
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
    out: &mut [Vec<f64>],
) {
    let eqs = SpaceEquations::new(entities, constraints);
    system::analytic_jacobian(&eqs, &system::all_indices(&eqs), params, col_of, out);
}

/// The plane gauge fix: rescale `(nx, ny, nz, d)` to a unit normal (scaling
/// the offset with it).
///
/// Uniform scaling of `(n, d)` leaves every constraint residual unchanged
/// (plane equations are homogeneous in the pair; oriented quantities divide
/// by `|n|`), so the 3D plane has three physical degrees of freedom — two
/// orientation plus one offset — while its storage has four parameters. The
/// extra direction is pure gauge. Pinning the representation to `|n| = 1`
/// keeps iterates and diagnosis well-conditioned; constraints that read the
/// plane through normalized quantities are unaffected.
///
/// Planes with any fixed parameter are left untouched: normalization would
/// move a masked value, breaking the fix mask contract.
pub fn normalize_planes(params: &mut [f64], entities: &[SpaceEntity], fixed: &[bool]) {
    for e in entities {
        if e.kind != SpaceEntityKind::Plane {
            continue;
        }
        let idx = [e.param(0), e.param(1), e.param(2), e.param(3)];
        if idx.iter().any(|&i| fixed.get(i).copied().unwrap_or(true)) {
            continue;
        }
        normalize_plane_at(params, idx);
    }
}

/// [`normalize_planes`] restricted to the planes that own at least one
/// parameter in `touched` — the incremental-solve form that leaves
/// unselected components bitwise untouched.
pub fn normalize_planes_touched(
    params: &mut [f64],
    entities: &[SpaceEntity],
    fixed: &[bool],
    touched: &[usize],
) {
    let mut is_touched = vec![false; params.len()];
    for &p in touched {
        if let Some(flag) = is_touched.get_mut(p) {
            *flag = true;
        }
    }
    for e in entities {
        if e.kind != SpaceEntityKind::Plane {
            continue;
        }
        let idx = [e.param(0), e.param(1), e.param(2), e.param(3)];
        if !idx.iter().any(|&i| is_touched[i]) {
            continue;
        }
        if idx.iter().any(|&i| fixed.get(i).copied().unwrap_or(true)) {
            continue;
        }
        normalize_plane_at(params, idx);
    }
}

fn normalize_plane_at(params: &mut [f64], idx: [usize; 4]) {
    let (nx, ny, nz, d) = (
        params[idx[0]],
        params[idx[1]],
        params[idx[2]],
        params[idx[3]],
    );
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    // Skip the division when the representation is already unit to
    // roundoff: re-normalizing an already normalized plane would perturb
    // the last bits and break the bitwise warm-start idempotence.
    if len > SOLVER_NORM_TOL && (len - 1.0).abs() > 1.0e-12 {
        let inv = 1.0 / len;
        params[idx[0]] = nx * inv;
        params[idx[1]] = ny * inv;
        params[idx[2]] = nz * inv;
        params[idx[3]] = d * inv;
    }
}

/// Run the 3D GCS solver over every connected component of the constraint
/// graph.
///
/// - `params`: full parameter vector (modified in-place for free params).
/// - `fixed`: mask — `fixed[i] == true` means `params[i]` is held constant.
/// - `entities`: entity metadata (for constraint evaluation).
/// - `constraints`: list of 3D constraints.
///
/// Uses the default [`SolveOptions`]: a component that does not converge
/// gets its starting parameters back (see
/// [`FailurePolicy`][crate::numerics::FailurePolicy]). Planes are still
/// left in their normalized `(unit normal, offset)` representation, which
/// describes the same plane.
pub fn solve_space(
    params: &mut [f64],
    fixed: &[bool],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
) -> SolveResult {
    solve_space_with(
        params,
        fixed,
        entities,
        constraints,
        &SolveOptions::default(),
    )
}

/// [`solve_space`] with explicit [`SolveOptions`].
pub fn solve_space_with(
    params: &mut [f64],
    fixed: &[bool],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
    options: &SolveOptions,
) -> SolveResult {
    // Gauge fix before iterating: the solver's variable space is the
    // (unit normal, offset) representation.
    normalize_planes(params, entities, fixed);
    let components = crate::space3d::graph::constraint_components(params, entities, constraints);
    solve_components_with(
        params,
        fixed,
        entities,
        constraints,
        &components,
        None,
        options,
    )
}

/// Solve only the selected connected components (`only` holds ordinals
/// into `components`; `None` selects all). Unselected components are not
/// iterated and their parameters are not touched — the localized-solving
/// entry point behind incremental re-solve (3D-0).
pub fn solve_components(
    params: &mut [f64],
    fixed: &[bool],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
    components: &[Vec<usize>],
    only: Option<&[usize]>,
) -> SolveResult {
    solve_components_with(
        params,
        fixed,
        entities,
        constraints,
        components,
        only,
        &SolveOptions::default(),
    )
}

/// [`solve_components`] with explicit [`SolveOptions`].
///
/// The reported status is the most severe one among the solved components;
/// the failure policy applies to each component separately.
pub fn solve_components_with(
    params: &mut [f64],
    fixed: &[bool],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
    components: &[Vec<usize>],
    only: Option<&[usize]>,
    options: &SolveOptions,
) -> SolveResult {
    solve_equations(
        &SpaceEquations::new(entities, constraints),
        params,
        fixed,
        components,
        only,
        options,
    )
}

/// [`solve_components_with`] over an equation view, which is how a sketch
/// leaves its suppressed and removed constraints out.
pub(crate) fn solve_equations(
    eqs: &SpaceEquations<'_>,
    params: &mut [f64],
    fixed: &[bool],
    components: &[Vec<usize>],
    only: Option<&[usize]>,
    options: &SolveOptions,
) -> SolveResult {
    let entities = eqs.entities;
    // Plane representation is pinned only for the planes the solved
    // components touched (bitwise guarantee for untouched components).
    system::solve_components(
        eqs,
        params,
        fixed,
        components,
        only,
        &CoreOptions::from(options),
        |params, touched| normalize_planes_touched(params, entities, fixed, touched),
    )
}

/// Residual report for the current parameters without iterating — the
/// return value of an incremental solve that selected no component (a warm
/// start already at a fixed point, or an edit that touches no constraint).
/// With nothing iterated, an unsatisfied state reports
/// [`SolveStatus::MaxIterations`] (a budget of zero iterations).
pub(crate) fn residual_report(
    eqs: &SpaceEquations<'_>,
    params: &[f64],
    tolerance: f64,
) -> SolveResult {
    let f = system::eval_residuals(eqs, &system::all_indices(eqs), params);
    let residual = numerics::rms(&f);
    let converged = residual < tolerance;
    SolveResult {
        converged,
        residual,
        iterations: 0,
        status: if converged {
            SolveStatus::Converged
        } else {
            SolveStatus::MaxIterations
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space3d::constraint::SpaceConstraint;
    use crate::space3d::entity::{SpaceEntity, SpaceEntityKind, SpacePointRef};
    use std::collections::HashMap;

    #[test]
    fn fixed_point_3d() {
        let entities = vec![SpaceEntity::new(SpaceEntityKind::SpacePoint, 0)];
        let constraints = vec![SpaceConstraint::Fixed {
            point: SpacePointRef::Point(0),
            x: 1.0,
            y: 2.0,
            z: 3.0,
        }];
        let mut params = vec![10.0_f64, 20.0, 30.0];
        let fixed = vec![false, false, false];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        assert!((params[0] - 1.0).abs() < 1e-7);
        assert!((params[1] - 2.0).abs() < 1e-7);
        assert!((params[2] - 3.0).abs() < 1e-7);
    }

    #[test]
    fn coincident_3d() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 3),
        ];
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::Point(0),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::Coincident(SpacePointRef::Point(0), SpacePointRef::Point(1)),
        ];
        let mut params = vec![0.0, 0.0, 0.0, 5.0, 5.0, 5.0];
        let fixed = vec![true, true, true, false, false, false];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        assert!((params[3] - 0.0).abs() < 1e-7);
        assert!((params[4] - 0.0).abs() < 1e-7);
        assert!((params[5] - 0.0).abs() < 1e-7);
    }

    #[test]
    fn point_distance_3d() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 3),
        ];
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::Point(0),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::PointDistance {
                p1: SpacePointRef::Point(0),
                p2: SpacePointRef::Point(1),
                distance: 5.0,
            },
        ];
        let mut params = vec![0.0, 0.0, 0.0, 3.0, 4.0, 0.0];
        let fixed = vec![true, true, true, false, false, false];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        let dx = params[3];
        let dy = params[4];
        let dz = params[5];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!((dist - 5.0).abs() < 1e-7, "distance = {dist}");
    }

    #[test]
    fn point_on_plane() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
            SpaceEntity::new(SpaceEntityKind::Plane, 3),
        ];
        // Plane: z = 10 (normal = (0,0,1), d = 10)
        // Fix the plane params and point x,y, leave z free for PointOnPlane to solve.
        let constraints = vec![SpaceConstraint::PointOnPlane {
            point: SpacePointRef::Point(0),
            plane: 1,
        }];
        let mut params = vec![1.0, 2.0, 0.0, 0.0, 0.0, 1.0, 10.0];
        // Fix plane params (indices 3-6) and point x,y (indices 0,1); z (index 2) is free.
        let fixed = vec![false, false, false, true, true, true, true];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        assert!(
            (params[2] - 10.0).abs() < 1e-7,
            "z should be 10, got {}",
            params[2]
        );
    }

    #[test]
    fn line_length_3d() {
        let entities = vec![SpaceEntity::new(SpaceEntityKind::SpaceLine, 0)];
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::LineStart(0),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::LineLength {
                line: 0,
                length: 3.0,
            },
        ];
        let mut params = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let fixed = vec![true, true, true, false, false, false];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        let dx = params[3];
        let dy = params[4];
        let dz = params[5];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!((len - 3.0).abs() < 1e-7, "length = {len}");
    }

    #[test]
    fn sphere_tangent_external() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::Sphere, 0),
            SpaceEntity::new(SpaceEntityKind::Sphere, 4),
        ];
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::SphereCenter(0),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::SphereRadius {
                sphere: 0,
                radius: 2.0,
            },
            SpaceConstraint::SphereRadius {
                sphere: 1,
                radius: 1.0,
            },
            SpaceConstraint::Fixed {
                point: SpacePointRef::SphereCenter(1),
                x: 5.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::SphereTangent {
                s1: 0,
                s2: 1,
                external: true,
            },
        ];
        let mut params = vec![0.0, 0.0, 0.0, 2.0, 5.0, 0.0, 0.0, 1.0];
        // Fix sphere 0 center + radius, sphere 1 center + radius
        let fixed = vec![true, true, true, true, true, true, true, true];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        // All params fixed → nothing to solve, and it must not panic. The
        // spheres sit 5 apart with radii 2 and 1, so the tangency does not
        // hold; 0.1.0 reported this as converged, now it says what it is.
        assert!(!result.converged);
        assert_eq!(result.status, SolveStatus::Stalled);
        assert_eq!(result.iterations, 0);
        assert_eq!(params, vec![0.0, 0.0, 0.0, 2.0, 5.0, 0.0, 0.0, 1.0]);
    }

    #[test]
    fn line_perpendicular_line() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpaceLine, 0),
            SpaceEntity::new(SpaceEntityKind::SpaceLine, 6),
        ];
        // Line 0: along X axis from (0,0,0) to (1,0,0)
        // Line 1: along Y axis from (0,0,0) to (0,1,0)
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::LineStart(0),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::Fixed {
                point: SpacePointRef::LineEnd(0),
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::Fixed {
                point: SpacePointRef::LineStart(1),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::LinePerpendicularLine(0, 1),
        ];
        let mut params = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let fixed = vec![
            true, true, true, true, true, true, true, true, true, false, false, false,
        ];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        // Line 1 should end up perpendicular to line 0 (along Y or Z)
        let dx2 = params[9];
        let dy2 = params[10];
        let dz2 = params[11];
        let dot = 1.0 * dx2 + 0.0 * dy2 + 0.0 * dz2;
        assert!(dot.abs() < 1e-7, "dot product should be ~0, got {dot}");
    }

    /// The 2-row constraints now report two equations; the solver must
    /// consume both rows (regression for the pre-3D-0 equation_count gap).
    #[test]
    fn parallel_direction_constraint_uses_two_rows() {
        let entities = vec![SpaceEntity::new(SpaceEntityKind::SpaceLine, 0)];
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::LineStart(0),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::LineParallelToDir {
                line: 0,
                dx: 1.0,
                dy: 0.0,
                dz: 0.0,
            },
        ];
        assert_eq!(constraints[1].equation_count(), 2);
        let mut params = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let fixed = vec![true, true, true, false, false, false];
        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "residual {}", result.residual);
        assert!(
            params[5].abs() < 1e-7,
            "z direction should be zero, got {}",
            params[5]
        );
    }

    /// 3D-0 acceptance measurement (ignored by default; run explicitly):
    /// `cargo test -p rcad-gcs --lib -- --ignored --nocapture analytic_jacobian_assembly_bench`
    /// Analytic assembly vs a finite-difference build at ~1500 parameters.
    #[test]
    #[ignore]
    fn analytic_jacobian_assembly_bench() {
        use std::time::Instant;

        let n = 500usize;
        let entities: Vec<SpaceEntity> = (0..n)
            .map(|i| SpaceEntity::new(SpaceEntityKind::SpacePoint, i * 3))
            .collect();
        let mut constraints = vec![SpaceConstraint::fix_point(0usize, 0.0, 0.0, 0.0)];
        for i in 0..n - 1 {
            constraints.push(SpaceConstraint::point_distance(i, i + 1, 1.0));
        }
        let params: Vec<f64> = (0..n)
            .flat_map(|i| [i as f64, (i % 7) as f64 * 0.1, 0.0_f64])
            .collect();
        let n_eq: usize = constraints.iter().map(|c| c.equation_count()).sum();
        let free: Vec<usize> = (0..params.len()).collect();
        let col_of: HashMap<usize, usize> = free.iter().enumerate().map(|(c, &p)| (p, c)).collect();

        let mut ja = vec![vec![0.0_f64; free.len()]; n_eq];
        let t0 = Instant::now();
        analytic_jacobian(&params, &col_of, &entities, &constraints, &mut ja);
        let analytic_us = t0.elapsed().as_secs_f64() * 1e6;

        let mut eval = |p: &[f64], out: &mut [f64]| {
            let v = eval_residuals(p, &entities, &constraints, n_eq);
            out.copy_from_slice(&v);
        };
        let mut jfd = vec![vec![0.0_f64; free.len()]; n_eq];
        let t1 = Instant::now();
        numerics::numerical_jacobian_into(&params, &free, n_eq, &mut eval, &mut jfd);
        let fd_us = t1.elapsed().as_secs_f64() * 1e6;

        let max_err = ja
            .iter()
            .zip(jfd.iter())
            .flat_map(|(a, b)| a.iter().zip(b.iter()))
            .map(|(x, y)| (x - y).abs())
            .fold(0.0_f64, f64::max);
        println!("params: {n_eq} equations x {} free", free.len());
        println!("analytic assembly: {analytic_us:.0} us");
        println!("finite differences: {fd_us:.0} us");
        println!(
            "speedup: {:.1}x, max |analytic - fd| = {max_err:.2e}",
            fd_us / analytic_us
        );
    }

    #[test]
    fn planes_are_normalized_by_the_gauge_fix() {
        // The gauge fix scales the offset with the normal: (6,-3,0,2)
        // becomes ((6,-3,0)/|n|, 2/|n|).
        let mut params = vec![6.0, -3.0, 0.0, 2.0];
        let fixed = vec![false; 4];
        let entities = vec![SpaceEntity::new(SpaceEntityKind::Plane, 0)];
        normalize_planes(&mut params, &entities, &fixed);
        let len = (params[0] * params[0] + params[1] * params[1] + params[2] * params[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "plane normal not unit: {len}");
        assert!(
            (params[0] * (2.0 / 6.0) - params[3]).abs() < 1e-12,
            "offset not scaled with normal: n={:?} d={}",
            &params[0..3],
            params[3]
        );

        // A solve leaves the representation unit-normal (the point and the
        // plane move, but the gauge stays pinned).
        let mut params = vec![0.0, 0.0, 0.0, 6.0, -3.0, 0.0, 2.0];
        let fixed = vec![false; 7];
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
            SpaceEntity::new(SpaceEntityKind::Plane, 3),
        ];
        let constraints = vec![SpaceConstraint::PointOnPlane {
            point: SpacePointRef::Point(0),
            plane: 1,
        }];
        let _ = solve_space(&mut params, &fixed, &entities, &constraints);
        let len = (params[3] * params[3] + params[4] * params[4] + params[5] * params[5]).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-12,
            "plane normal not unit after solve: {len}"
        );
    }
}
