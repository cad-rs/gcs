//! Newton-Raphson GCS solver for 2D sketch constraints.
//!
//! A thin family adapter: the 2D constraint set contributes the residuals
//! and the analytic Jacobian rows (each [`Constraint`] derivative stays in
//! lockstep with its residual via the generic finite-difference CI test);
//! the decomposition and the per-component solve are the shared pipeline of
//! the crate-internal `system` module over the Newton core of [`crate::numerics`].

use crate::constraint::Constraint;
use crate::edit::ConstraintState;
use crate::entity::{Entity, EntityKind};
use crate::numerics::CoreOptions;
use crate::system::{self, Equations, SparseRow};

pub use crate::numerics::{
    SolveOptions, SolveResult, SolveStatus, FD_H, LAMBDA, MAX_ITER, PIVOT_TOL, RESIDUAL_TOL,
    SOLVER_NORM_TOL,
};

/// The 2D constraint list as the pipeline's [`Equations`].
pub(crate) struct SketchEquations<'a> {
    entities: &'a [Entity],
    constraints: &'a [Constraint],
    /// The state of each constraint slot; a slot beyond the end of this
    /// list, or no list at all, counts as active.
    states: Option<&'a [ConstraintState]>,
}

impl<'a> SketchEquations<'a> {
    /// Every constraint in the list is enforced.
    pub(crate) fn new(entities: &'a [Entity], constraints: &'a [Constraint]) -> Self {
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
        entities: &'a [Entity],
        constraints: &'a [Constraint],
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

impl Equations for SketchEquations<'_> {
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

    /// Every 2D equation is a length (the angular relations are scaled by
    /// the lengths of their lines). The size of a component is the larger
    /// side of the box around the entities its active constraints refer
    /// to, circles and arcs taken whole.
    fn length_scale(&self, group: &[usize], params: &[f64]) -> Option<f64> {
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        let mut add = |x: f64, y: f64, r: f64| {
            lo[0] = lo[0].min(x - r);
            lo[1] = lo[1].min(y - r);
            hi[0] = hi[0].max(x + r);
            hi[1] = hi[1].max(y + r);
        };
        for &index in group {
            if !self.is_active(index) {
                continue;
            }
            for (id, _) in self.constraints[index].references() {
                let e = &self.entities[id];
                let p = |k: usize| params[e.param(k)];
                match e.kind {
                    EntityKind::Point => add(p(0), p(1), 0.0),
                    EntityKind::Line => {
                        add(p(0), p(1), 0.0);
                        add(p(2), p(3), 0.0);
                    }
                    EntityKind::Circle | EntityKind::Arc => add(p(0), p(1), p(2).abs()),
                }
            }
        }
        let size = (hi[0] - lo[0]).max(hi[1] - lo[1]);
        size.is_finite().then_some(size)
    }
}

/// Evaluate the stacked residual vector of `constraints` at `params`.
#[cfg(test)]
pub(crate) fn eval_residuals(
    params: &[f64],
    entities: &[Entity],
    constraints: &[Constraint],
    n_eq: usize,
) -> Vec<f64> {
    let eqs = SketchEquations::new(entities, constraints);
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
    entities: &[Entity],
    constraints: &[Constraint],
    out: &mut [Vec<f64>],
) {
    let eqs = SketchEquations::new(entities, constraints);
    system::analytic_jacobian(&eqs, &system::all_indices(&eqs), params, col_of, out);
}

/// Run the Newton-Raphson solver over every connected component of the
/// constraint graph (solver roadmap M3: full solve = per-component solve,
/// which is what makes incremental re-solve result-identical by
/// construction).
///
/// - `params`: full parameter vector (modified in-place for free params).
/// - `fixed`: mask — `fixed[i] == true` means `params[i]` is held constant.
/// - `entities`: entity metadata (for constraint evaluation).
/// - `constraints`: list of constraints.
///
/// Uses the default [`SolveOptions`]: a component that does not converge
/// gets its starting parameters back (see
/// [`FailurePolicy`][crate::numerics::FailurePolicy]).
pub fn solve(
    params: &mut [f64],
    fixed: &[bool],
    entities: &[Entity],
    constraints: &[Constraint],
) -> SolveResult {
    solve_with(
        params,
        fixed,
        entities,
        constraints,
        &SolveOptions::default(),
    )
}

/// [`solve`] with explicit [`SolveOptions`].
pub fn solve_with(
    params: &mut [f64],
    fixed: &[bool],
    entities: &[Entity],
    constraints: &[Constraint],
    options: &SolveOptions,
) -> SolveResult {
    let components = crate::graph::constraint_components(params, entities, constraints);
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
/// entry point behind the sketcher's incremental rebuild (M3.2).
pub fn solve_components(
    params: &mut [f64],
    fixed: &[bool],
    entities: &[Entity],
    constraints: &[Constraint],
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
    entities: &[Entity],
    constraints: &[Constraint],
    components: &[Vec<usize>],
    only: Option<&[usize]>,
    options: &SolveOptions,
) -> SolveResult {
    let eqs = SketchEquations::new(entities, constraints);
    system::solve_components(
        &eqs,
        params,
        fixed,
        components,
        only,
        &CoreOptions::from(options),
        |_, _| {},
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numerics::gaussian_elimination;

    #[test]
    fn debug_fixed_point_inline() {
        use crate::constraint::Constraint;
        use crate::entity::{Entity, EntityKind, PointRef};

        let entities = vec![Entity::new(EntityKind::Point, 0)];
        let constraints = vec![Constraint::Fixed {
            point: PointRef::Point(0),
            x: 0.0,
            y: 0.0,
        }];
        let mut params = vec![3.0_f64, 7.0];
        let fixed = vec![false, false];

        let n_eq = 2usize;
        let f = eval_residuals(&params, &entities, &constraints, n_eq);
        println!("initial f = {:?}", f);

        let result = solve(&mut params, &fixed, &entities, &constraints);
        println!("params after = {:?}", params);
        println!(
            "converged={} residual={}",
            result.converged, result.residual
        );
        assert!(result.converged, "not converged: {}", result.residual);
        assert!((params[0] - 0.0).abs() < 1e-7);
        assert!((params[1] - 0.0).abs() < 1e-7);
        let _ = (f, fixed);
    }

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

    // ─── Constraint solver edge case tests ───────────────────────────────────

    #[test]
    fn under_constrained_system_converges() {
        // One point with only Fixed X — under-constrained (Y is free).
        // The solver should still converge (Tikhonov regularization handles this).
        use crate::constraint::Constraint;
        use crate::entity::{Entity, EntityKind, PointRef};

        let entities = vec![Entity::new(EntityKind::Point, 0)];
        let constraints = vec![Constraint::Fixed {
            point: PointRef::Point(0),
            x: 5.0,
            y: 0.0,
        }];
        let mut params = vec![0.0, 3.0]; // x=0 (wrong), y=3 (free)
        let fixed = vec![false, false];

        let result = solve(&mut params, &fixed, &entities, &constraints);
        assert!(
            result.converged,
            "under-constrained should converge: residual={}",
            result.residual
        );
        assert!(
            (params[0] - 5.0).abs() < 1e-7,
            "x should be 5.0, got {}",
            params[0]
        );
    }

    #[test]
    fn over_constrained_consistent_system() {
        // A point fixed at (3,4) by two Fixed constraints (redundant but consistent).
        use crate::constraint::Constraint;
        use crate::entity::{Entity, EntityKind, PointRef};

        let entities = vec![Entity::new(EntityKind::Point, 0)];
        let constraints = vec![
            Constraint::Fixed {
                point: PointRef::Point(0),
                x: 3.0,
                y: 4.0,
            },
            Constraint::Fixed {
                point: PointRef::Point(0),
                x: 3.0,
                y: 4.0,
            }, // duplicate
        ];
        let mut params = vec![0.0, 0.0];
        let fixed = vec![false, false];

        let result = solve(&mut params, &fixed, &entities, &constraints);
        assert!(
            result.converged,
            "over-constrained consistent should converge: residual={}, iterations={}",
            result.residual, result.iterations
        );
        assert!(
            (params[0] - 3.0).abs() < 1e-6,
            "x should be 3.0, got {}",
            params[0]
        );
        assert!(
            (params[1] - 4.0).abs() < 1e-6,
            "y should be 4.0, got {}",
            params[1]
        );
    }

    #[test]
    fn over_constrained_inconsistent_system() {
        // Two points with conflicting constraints:
        // p0 fixed at (0,0), p1 fixed at (1,1), but PointDistance=5 (impossible, max dist=sqrt(2)).
        use crate::constraint::Constraint;
        use crate::entity::{Entity, EntityKind, PointRef};

        let entities = vec![
            Entity::new(EntityKind::Point, 0),
            Entity::new(EntityKind::Point, 1),
        ];
        let constraints = vec![
            Constraint::Fixed {
                point: PointRef::Point(0),
                x: 0.0,
                y: 0.0,
            },
            Constraint::Fixed {
                point: PointRef::Point(1),
                x: 1.0,
                y: 1.0,
            },
            Constraint::PointDistance {
                p1: PointRef::Point(0),
                p2: PointRef::Point(1),
                distance: 5.0,
            },
        ];
        let mut params = vec![0.0, 0.0, 1.0, 1.0];
        let fixed = vec![false, false, false, false];

        let result = solve(&mut params, &fixed, &entities, &constraints);
        // This system is inconsistent — it may not fully converge but should not panic.
        // The solver should report non-convergence.
        assert!(!result.converged, "inconsistent system should not converge");
    }

    #[test]
    fn grid_system_converges() {
        // A 4x4 grid of points with distance constraints in both directions:
        // 32 parameters, 24 distance rows plus three anchored points.
        //
        // This test used to be `large_system_makes_progress` and only asked
        // for a residual below 0.5. It never converged because the entities
        // were laid out with `param_start = i` instead of `2 * i`, so
        // neighbouring points shared a coordinate and the system was
        // inconsistent. With the layout fixed it converges outright.
        use crate::constraint::Constraint;
        use crate::entity::{Entity, EntityKind, PointRef};

        let n = 4; // 4x4 grid = 16 points
        let spacing = 1.0;
        let n_points = n * n;

        let entities: Vec<Entity> = (0..n_points)
            .map(|i| Entity::new(EntityKind::Point, 2 * i))
            .collect();

        let mut constraints = Vec::new();

        // Horizontal distance constraints
        for r in 0..n {
            for c in 0..n - 1 {
                let i = r * n + c;
                let j = r * n + c + 1;
                constraints.push(Constraint::PointDistance {
                    p1: PointRef::Point(i),
                    p2: PointRef::Point(j),
                    distance: spacing,
                });
            }
        }

        // Vertical distance constraints
        for r in 0..n - 1 {
            for c in 0..n {
                let i = r * n + c;
                let j = (r + 1) * n + c;
                constraints.push(Constraint::PointDistance {
                    p1: PointRef::Point(i),
                    p2: PointRef::Point(j),
                    distance: spacing,
                });
            }
        }

        // Fix corner point at origin (removes translation DOF)
        constraints.push(Constraint::Fixed {
            point: PointRef::Point(0),
            x: 0.0,
            y: 0.0,
        });

        // Fix point (0,1) on Y axis (removes rotation DOF)
        constraints.push(Constraint::Fixed {
            point: PointRef::Point(n),
            x: 0.0,
            y: spacing,
        });

        // Fix point (1,0) on X axis (removes reflection DOF)
        constraints.push(Constraint::Fixed {
            point: PointRef::Point(1),
            x: spacing,
            y: 0.0,
        });

        // Initialize params: place all points on a regular grid with perturbation
        let mut params: Vec<f64> = vec![0.0; n_points * 2];
        for r in 0..n {
            for c in 0..n {
                let i = r * n + c;
                params[2 * i] = c as f64 * spacing + 0.01;
                params[2 * i + 1] = r as f64 * spacing + 0.01;
            }
        }

        let fixed = vec![false; n_points * 2];

        let result = solve(&mut params, &fixed, &entities, &constraints);
        assert!(
            result.converged,
            "grid should converge: status={:?}, residual={}, iterations={}",
            result.status, result.residual, result.iterations
        );
        assert_eq!(result.status, SolveStatus::Converged);
        // The three anchors landed where they were pinned.
        for (point, x, y) in [(0, 0.0, 0.0), (n, 0.0, spacing), (1, spacing, 0.0)] {
            assert!((params[2 * point] - x).abs() < 1e-7, "point {point} x");
            assert!((params[2 * point + 1] - y).abs() < 1e-7, "point {point} y");
        }
    }

    #[test]
    fn empty_constraints_converge_immediately() {
        use crate::entity::{Entity, EntityKind};

        let entities = vec![Entity::new(EntityKind::Point, 0)];
        let constraints: Vec<Constraint> = vec![];
        let mut params = vec![42.0, 99.0];
        let fixed = vec![false, false];

        let result = solve(&mut params, &fixed, &entities, &constraints);
        assert!(
            result.converged,
            "empty constraints should trivially converge"
        );
        assert_eq!(result.residual, 0.0);
        assert_eq!(result.iterations, 0);
    }

    #[test]
    fn all_fixed_params_converge_immediately() {
        use crate::entity::{Entity, EntityKind, PointRef};

        let entities = vec![Entity::new(EntityKind::Point, 0)];
        let constraints = vec![Constraint::Coincident(
            PointRef::Point(0),
            PointRef::Point(0),
        )];
        let mut params = vec![1.0, 2.0];
        let fixed = vec![true, true]; // all params fixed

        let result = solve(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "all-fixed should trivially converge");
        assert_eq!(params, vec![1.0, 2.0], "fixed params should not change");
    }
}
