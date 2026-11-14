//! What a failed solve leaves behind, and what it reports.
//!
//! The default contract: a connected component that does not converge gets
//! its starting parameters back, so a host never has to repair geometry
//! after a failed solve. `FailurePolicy::KeepLast` opts out per call.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::entity::PointRef;
use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::SpaceSketch;
use rcad_gcs::{FailurePolicy, Sketch, SolveOptions, SolveStatus};

/// Two points pinned 1.41 apart and asked to be 5 apart: no solution.
fn conflicting_sketch() -> Sketch {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(p2, 1.0, 1.0));
    sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
    sk
}

fn keep_last() -> SolveOptions {
    SolveOptions {
        on_failure: FailurePolicy::KeepLast,
        ..SolveOptions::default()
    }
}

#[test]
fn failed_solve_restores_the_start_by_default() {
    let mut sk = conflicting_sketch();
    let before = sk.params_snapshot();

    let result = sk.solve();

    assert!(!result.converged);
    // A conflict settles at its least-squares compromise, where no step of
    // any length lowers the residual further. The solve says so, well
    // before the iteration limit.
    assert_eq!(result.status, SolveStatus::Stalled);
    assert!(
        result.iterations < 20,
        "{} iterations to find out",
        result.iterations
    );
    assert_eq!(
        sk.params_snapshot(),
        before,
        "a failed solve must leave the geometry untouched"
    );
    // The reported residual is the one of the state left behind: only the
    // distance row is unsatisfied, by 5 − √2, over five equations.
    let expected = (5.0 - 2.0_f64.sqrt()) / 5.0_f64.sqrt();
    assert!(
        (result.residual - expected).abs() < 1e-12,
        "residual {} vs {expected}",
        result.residual
    );
}

#[test]
fn keep_last_leaves_the_last_iterate() {
    let mut restored = conflicting_sketch();
    let start_residual = restored.solve().residual;

    let mut sk = conflicting_sketch();
    let before = sk.params_snapshot();
    let result = sk.solve_with(&keep_last());

    assert!(!result.converged);
    assert_ne!(
        sk.params_snapshot(),
        before,
        "KeepLast keeps whatever the iteration reached"
    );
    assert!(
        result.residual < start_residual,
        "the least-squares compromise is closer than the start: {} vs {start_residual}",
        result.residual
    );
}

#[test]
fn restore_is_per_component() {
    // Component A is solvable, component B is the conflicting pair.
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 0.0);
    let a2 = sk.add_point(2.0, 0.5);
    sk.add_constraint(Constraint::fix_point(a1, 0.0, 0.0));
    sk.add_constraint(Constraint::point_distance(a1, a2, 3.0));

    let b1 = sk.add_point(10.0, 0.0);
    let b2 = sk.add_point(11.0, 1.0);
    sk.add_constraint(Constraint::fix_point(b1, 10.0, 0.0));
    sk.add_constraint(Constraint::fix_point(b2, 11.0, 1.0));
    sk.add_constraint(Constraint::point_distance(b1, b2, 5.0));

    assert_eq!(sk.components().len(), 2);
    let before = sk.params_snapshot();
    let result = sk.solve();

    assert!(
        !result.converged,
        "the conflicting component fails the solve"
    );
    assert_ne!(result.status, SolveStatus::Converged);

    // Component A was solved: a2 now sits 3 away from the origin.
    let a2_now = sk.point_coords(PointRef::Point(a2));
    assert!(
        (a2_now.length() - 3.0).abs() < 1e-9,
        "a2 at distance {}",
        a2_now.length()
    );
    // Component B was put back exactly.
    let after = sk.params_snapshot();
    assert_eq!(&after[4..8], &before[4..8], "component B must be untouched");
}

#[test]
fn converged_solve_reports_converged_status() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.3, -0.2);
    let p2 = sk.add_point(2.0, 3.0);
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
    sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));

    let result = sk.solve();
    assert!(result.converged);
    assert_eq!(result.status, SolveStatus::Converged);
}

/// The apex of a 3-4-5 triangle on a held base, started off its place:
/// after the fixed points are eliminated two distance rows remain, which
/// take several Newton iterations. (A single distance from one fixed point
/// does not: once that point is in place the row is linear along the step.)
fn triangle_apex() -> Sketch {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.3, -0.2);
    let b = sk.add_point(3.7, 0.4);
    let apex = sk.add_point(0.8, 2.0);
    sk.add_constraint(Constraint::fix_point(a, 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(b, 4.0, 0.0));
    sk.add_constraint(Constraint::point_distance(a, apex, 3.0));
    sk.add_constraint(Constraint::point_distance(b, apex, 5.0));
    sk
}

#[test]
fn iteration_limit_is_reported_and_restores() {
    let mut sk = triangle_apex();
    let before = sk.params_snapshot();

    // One iteration is not enough for the nonlinear distance rows.
    let result = sk.solve_with(&SolveOptions {
        max_iterations: 1,
        ..SolveOptions::default()
    });

    assert!(!result.converged);
    assert_eq!(result.status, SolveStatus::MaxIterations);
    assert_eq!(result.iterations, 1);
    assert_eq!(sk.params_snapshot(), before);
}

#[test]
fn looser_tolerance_is_honoured() {
    let build = triangle_apex;
    let strict = build().solve();
    let loose = build().solve_with(&SolveOptions {
        tolerance: 1e-2,
        ..SolveOptions::default()
    });
    assert!(strict.converged && loose.converged);
    assert!(
        loose.iterations < strict.iterations,
        "a looser tolerance stops earlier: {} vs {}",
        loose.iterations,
        strict.iterations
    );
}

#[test]
fn space_sketch_restores_on_failure_too() {
    let mut sk = SpaceSketch::new();
    let p1 = sk.add_point(0.0, 0.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0, 1.0);
    sk.add_constraint(SpaceConstraint::fix_point(p1, 0.0, 0.0, 0.0));
    sk.add_constraint(SpaceConstraint::fix_point(p2, 1.0, 1.0, 1.0));
    sk.add_constraint(SpaceConstraint::point_distance(p1, p2, 9.0));
    let before = sk.params_snapshot();

    let result = sk.solve();
    assert!(!result.converged);
    assert_ne!(result.status, SolveStatus::Converged);
    assert_eq!(sk.params_snapshot(), before);

    let kept = sk.solve_with(&keep_last());
    assert!(!kept.converged);
    assert_ne!(sk.params_snapshot(), before);
}

#[test]
fn unsatisfied_equations_with_nothing_free_are_not_converged() {
    // A fully fixed line of length 3 asked to be 7 long: nothing can move.
    let mut sk = Sketch::new();
    let line = sk.add_line(0.0, 0.0, 3.0, 0.0);
    sk.fix_entity(line);
    sk.add_constraint(Constraint::LineLength { line, length: 7.0 });

    let result = sk.solve();
    assert!(!result.converged, "0.1.0 reported this as converged");
    assert_eq!(result.status, SolveStatus::Stalled);
    assert_eq!(result.iterations, 0);
    assert!((result.residual - 4.0).abs() < 1e-12);
    assert_eq!(sk.diagnose().conflicting, vec![0]);

    // The same fixed line asked for the length it has is satisfied.
    let mut sk = Sketch::new();
    let line = sk.add_line(0.0, 0.0, 3.0, 0.0);
    sk.fix_entity(line);
    sk.add_constraint(Constraint::LineLength { line, length: 3.0 });
    let result = sk.solve();
    assert!(result.converged);
    assert_eq!(result.status, SolveStatus::Converged);
}

// ── Rigid-pose assemblies ───────────────────────────────────────────────────

mod assembly {
    use glam::{DQuat, DVec3};
    use rcad_gcs::rigid3d::{Constraint, ConstraintSolver, ConstraintSystem, PoseSolveOptions};
    use rcad_gcs::{FailurePolicy, SolveStatus};

    /// A ground (body 0), a part that can be mated to it (body 1), and a
    /// part pinned at two places at once (body 2), which cannot be solved.
    /// The two parts share no mate, so they are separate clusters.
    fn half_solvable() -> ConstraintSystem {
        let mut system = ConstraintSystem::with_components(3);
        system.add_constraint(Constraint::fixed_with_orientation(
            0,
            0,
            DVec3::ZERO,
            DQuat::IDENTITY,
        ));
        system.add_constraint(Constraint::coincident(
            1,
            0,
            DVec3::new(2.0, 0.0, 0.0),
            1,
            DVec3::new(0.0, 1.0, 0.0),
        ));
        system.add_constraint(Constraint::fixed(2, 2, DVec3::new(5.0, 0.0, 0.0)));
        system.add_constraint(Constraint::fixed(3, 2, DVec3::new(5.0, 4.0, 0.0)));
        system
    }

    /// Where the parts are before the solve: nowhere in particular.
    fn start() -> (Vec<DVec3>, Vec<DQuat>) {
        (
            vec![
                DVec3::ZERO,
                DVec3::new(0.3, -0.7, 1.1),
                DVec3::new(4.0, 1.0, -2.0),
            ],
            vec![
                DQuat::IDENTITY,
                DQuat::from_rotation_x(0.4),
                DQuat::from_rotation_z(-0.9),
            ],
        )
    }

    #[test]
    fn a_failed_cluster_is_put_back_and_a_solved_one_kept() {
        let system = half_solvable();
        let (positions, rotations) = start();
        let result =
            ConstraintSolver::new().solve_from(&system, positions.clone(), rotations.clone());

        assert!(!result.success);
        assert_ne!(result.status, SolveStatus::Converged);

        // The part that cannot be solved is exactly where it was.
        assert_eq!(result.positions[2], positions[2]);
        assert_eq!(result.rotations[2], rotations[2]);

        // The part that can be solved was: its mate point is on the ground's.
        let mated = result.positions[1] + result.rotations[1] * DVec3::new(0.0, 1.0, 0.0);
        assert!(
            (mated - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-6,
            "{mated}"
        );

        // The error is that of the poses handed back: body 2 misses its two
        // pins by what it missed them by at the start.
        let miss = |pin: DVec3| (positions[2] - pin).length_squared();
        let expected = (miss(DVec3::new(5.0, 0.0, 0.0)) + miss(DVec3::new(5.0, 4.0, 0.0))).sqrt();
        assert!(
            (result.error - expected).abs() < 1e-6,
            "{} vs {expected}",
            result.error
        );
    }

    #[test]
    fn keep_last_leaves_the_compromise() {
        let system = half_solvable();
        let (positions, rotations) = start();
        let solver = ConstraintSolver::with_options(PoseSolveOptions {
            on_failure: FailurePolicy::KeepLast,
            ..PoseSolveOptions::default()
        });
        let result = solver.solve_from(&system, positions, rotations);

        assert!(!result.success);
        // Halfway between the two pins is the best that can be done.
        let midway = DVec3::new(5.0, 2.0, 0.0);
        assert!(
            (result.positions[2] - midway).length() < 1e-3,
            "{}",
            result.positions[2]
        );
    }

    #[test]
    fn a_solved_assembly_reports_converged() {
        let mut system = half_solvable();
        system.constraints.pop();
        let (positions, rotations) = start();
        let result = ConstraintSolver::new().solve_from(&system, positions, rotations);
        assert!(result.success, "error {:.3e}", result.error);
        assert_eq!(result.status, SolveStatus::Converged);
        assert!((result.positions[2] - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-9);
    }
}
