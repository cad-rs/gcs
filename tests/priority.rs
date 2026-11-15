//! Constraint priorities (P5.2): a soft constraint gives way.
//!
//! Level 0 is a constraint that has to hold; a higher level holds only if
//! it can, and gives way to every lower level. A solve takes the levels in
//! turn and stops at the first it cannot meet, leaving the geometry where
//! the level below put it and reporting that level.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::{FailurePolicy, Sketch, SolveOptions};

/// `a` pinned at the origin and `b` free on the x axis: one hard
/// constraint each, and room for a soft distance to decide where `b` is.
/// Returns the sketch and the two points.
fn axis() -> (Sketch, usize, usize) {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(3.0, 0.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Point(a), 0.0, 0.0));
    sk.add_constraint(Constraint::distance_y(
        PointRef::Point(a),
        PointRef::Point(b),
        0.0,
    ));
    (sk, a, b)
}

/// A soft constraint that can be met is met, and the level says so.
#[test]
fn a_soft_constraint_that_can_be_met_is() {
    let (mut sk, a, b) = axis();
    let span = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        5.0,
    ));
    assert!(sk.set_constraint_priority(span, 1));
    assert_eq!(sk.constraint_priority(span), 1);
    let result = sk.solve();
    assert!(result.converged, "residual {:.3e}", result.residual);
    assert_eq!(result.satisfied_level, 1, "the soft level was within reach");
    assert!(
        (sk.point_coords(sk.point(b)).x.abs() - 5.0).abs() < 1e-9,
        "b is at {:?}",
        sk.point_coords(sk.point(b))
    );
}

/// A soft constraint that cannot be met is left out, the solve converges
/// all the same, and the geometry is the one the hard constraints give.
#[test]
fn a_soft_constraint_that_cannot_be_met_gives_way() {
    let (mut sk, a, b) = axis();
    // b is pinned by a hard constraint, so the soft distance is hopeless.
    // `distance_x` is signed (`x1 − x2`), so b is three to the right of a.
    sk.add_constraint(Constraint::distance_x(
        PointRef::Point(b),
        PointRef::Point(a),
        3.0,
    ));
    let span = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        5.0,
    ));
    assert!(sk.set_constraint_priority(span, 1));
    let result = sk.solve();
    assert!(
        result.converged,
        "a hopeless wish must not fail the solve: {:?} {:.3e}",
        result.status, result.residual
    );
    assert_eq!(result.satisfied_level, 0, "only the hard level holds");
    let at = sk.point_coords(sk.point(b));
    assert!(
        (at.x - 3.0).abs() < 1e-9 && at.y.abs() < 1e-9,
        "b should be where the hard constraints put it, not {at:?}"
    );
}

/// A level that cannot be met leaves the geometry where the level below
/// put it, whatever a failed solve is told to leave behind: not meeting a
/// wish is not a failure, so `FailurePolicy::KeepLast` has no say in it.
#[test]
fn a_level_that_cannot_be_met_keeps_the_level_below() {
    let (mut sk, a, b) = axis();
    let five = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        5.0,
    ));
    let eight = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        8.0,
    ));
    assert!(sk.set_constraint_priority(five, 1));
    assert!(sk.set_constraint_priority(eight, 2));
    let options = SolveOptions {
        on_failure: FailurePolicy::KeepLast,
        ..SolveOptions::default()
    };
    let result = sk.solve_with(&options);
    assert!(result.converged, "residual {:.3e}", result.residual);
    assert_eq!(result.satisfied_level, 1);
    assert!(
        (sk.point_coords(sk.point(b)).x.abs() - 5.0).abs() < 1e-9,
        "the level-2 iterate was kept: b is at {:?}",
        sk.point_coords(sk.point(b))
    );
}

/// The levels are taken in turn and the walk stops at the first one that
/// cannot be met: the level below is what the geometry is left at.
#[test]
fn the_levels_are_taken_in_turn() {
    let (mut sk, a, b) = axis();
    let five = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        5.0,
    ));
    let eight = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        8.0,
    ));
    // Levels need not be consecutive; only their order matters.
    assert!(sk.set_constraint_priority(five, 1));
    assert!(sk.set_constraint_priority(eight, 7));
    let result = sk.solve();
    assert!(result.converged, "residual {:.3e}", result.residual);
    assert_eq!(result.satisfied_level, 1, "level 7 contradicts level 1");
    assert!(
        (sk.point_coords(sk.point(b)).x.abs() - 5.0).abs() < 1e-9,
        "b is at {:?}",
        sk.point_coords(sk.point(b))
    );
    // The other way round: the weaker wish is the one that can be met.
    assert!(sk.set_constraint_priority(five, 7));
    assert!(sk.set_constraint_priority(eight, 1));
    let result = sk.solve();
    assert!(result.converged, "residual {:.3e}", result.residual);
    assert_eq!(result.satisfied_level, 1);
    assert!(
        (sk.point_coords(sk.point(b)).x.abs() - 8.0).abs() < 1e-9,
        "b is at {:?}",
        sk.point_coords(sk.point(b))
    );
}

/// A suppressed soft constraint takes no part, in either direction.
#[test]
fn a_suppressed_soft_constraint_takes_no_part() {
    let (mut sk, a, b) = axis();
    let span = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        5.0,
    ));
    assert!(sk.set_constraint_priority(span, 1));
    assert!(sk.set_constraint_enabled(span, false));
    let result = sk.solve();
    assert!(result.converged);
    assert_eq!(
        result.satisfied_level, 0,
        "a suppressed level is not a level"
    );
    assert!(
        (sk.point_coords(sk.point(b)).x - 3.0).abs() < 1e-9,
        "b should not have moved: {:?}",
        sk.point_coords(sk.point(b))
    );
    // It still remembers its level, and counts again when enabled.
    assert_eq!(sk.constraint_priority(span), 1);
    assert!(sk.set_constraint_enabled(span, true));
    assert_eq!(sk.solve().satisfied_level, 1);
}

/// A sketch with no priorities solves exactly as it did, and says level 0.
#[test]
fn a_sketch_without_priorities_is_unchanged() {
    let (mut sk, a, b) = axis();
    sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        5.0,
    ));
    let mut same = sk.clone();
    let plain = sk.solve();
    assert_eq!(plain.satisfied_level, 0);
    // Setting level 0 explicitly is no priority at all.
    for index in 0..same.constraints_len() {
        assert!(same.set_constraint_priority(index, 0));
    }
    let explicit = same.solve();
    assert_eq!(explicit.satisfied_level, 0);
    assert_eq!(
        sk.params_snapshot(),
        same.params_snapshot(),
        "level 0 everywhere changed the solution"
    );
    assert_eq!(plain.iterations, explicit.iterations);
}

/// An index that is not a slot is refused.
#[test]
fn an_index_that_is_not_a_slot_is_refused() {
    let (mut sk, _, _) = axis();
    let slots = sk.constraints_len();
    assert!(!sk.set_constraint_priority(slots, 1));
    assert_eq!(sk.constraint_priority(slots), 0);
    assert!(sk.set_constraint_priority(slots - 1, 3));
    assert_eq!(sk.constraint_priority(slots - 1), 3);
}

/// A level that holds only because a lower one was met: the soft levels
/// are added to the levels below, not solved on their own.
#[test]
fn a_level_is_added_to_the_ones_below() {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(4.0, 0.0);
    let c = sk.add_point(4.0, 3.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Point(a), 0.0, 0.0));
    // Hard: the three points make a right angle at b.
    sk.add_constraint(Constraint::distance_y(
        PointRef::Point(a),
        PointRef::Point(b),
        0.0,
    ));
    sk.add_constraint(Constraint::distance_x(
        PointRef::Point(b),
        PointRef::Point(c),
        0.0,
    ));
    // Level 1 fixes the leg a–b, level 2 asks for b–c to match it.
    let leg = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        5.0,
    ));
    let other = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(b),
        PointRef::Point(c),
        5.0,
    ));
    assert!(sk.set_constraint_priority(leg, 1));
    assert!(sk.set_constraint_priority(other, 2));
    let result = sk.solve();
    assert!(result.converged, "residual {:.3e}", result.residual);
    assert_eq!(result.satisfied_level, 2, "both wishes fit");
    let (at_b, at_c) = (sk.point_coords(sk.point(b)), sk.point_coords(sk.point(c)));
    assert!((at_b.x.abs() - 5.0).abs() < 1e-9, "b is at {at_b:?}");
    assert!(
        ((at_c.y - at_b.y).abs() - 5.0).abs() < 1e-9,
        "c is at {at_c:?}"
    );
    // Both are still what they were as hard constraints would have them.
    assert!(sk.diagnose().conflicting.is_empty());
    let _ = Horizontal(a);
}
