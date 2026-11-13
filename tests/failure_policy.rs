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
    // A conflict settles toward its least-squares compromise and runs out
    // of iterations there.
    assert_eq!(result.status, SolveStatus::MaxIterations);
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

#[test]
fn iteration_limit_is_reported_and_restores() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.3, -0.2);
    let p2 = sk.add_point(2.0, 3.0);
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
    sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
    let before = sk.params_snapshot();

    // One iteration is not enough for the nonlinear distance row.
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
    let build = || {
        let mut sk = Sketch::new();
        let p1 = sk.add_point(0.3, -0.2);
        let p2 = sk.add_point(2.0, 3.0);
        sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
        sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
        sk
    };
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
