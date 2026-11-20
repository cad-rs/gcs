//! Dimension jumps (P5.4): a value a solve cannot go straight to is
//! walked there a step at a time.
//!
//! The walk only happens after a solve that did not arrive, so a sketch
//! that solves is solved exactly as it was: `a_jump_a_solve_makes_is_left
//! _alone` is the guard for that. What changes is what a host gets when
//! the jump is too big — a value reached continuously instead of a
//! mirrored sketch, or an honest "this is how far it got".

use rcad_gcs::constraint::Constraint;
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::{EntityId, PointRef};
use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::entity::SpacePointRef;
use rcad_gcs::space3d::SpaceSketch;
use rcad_gcs::{Sketch, SolveStatus};

/// A rectangle with horizontal and vertical sides, its width driven by a
/// signed x distance and its height by a signed y distance, so that
/// asking for the other sign asks for the rectangle inside out. Returns
/// the sketch, the bottom line, and the width and height dimensions.
fn driven_rectangle(width: f64, height: f64) -> (Sketch, EntityId, usize, usize) {
    let mut sk = Sketch::new();
    let origin = sk.add_point(0.0, 0.0);
    let bottom = sk.add_line(0.0, 0.0, width, 0.0);
    let right = sk.add_line(width, 0.0, width, height);
    let top = sk.add_line(width, height, 0.0, height);
    let left = sk.add_line(0.0, height, 0.0, 0.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Point(origin), 0.0, 0.0));
    sk.add_constraint(Coincident(
        PointRef::LineStart(bottom),
        PointRef::Point(origin),
    ));
    for (end, start) in [(bottom, right), (right, top), (top, left), (left, bottom)] {
        sk.add_constraint(Coincident(
            PointRef::LineEnd(end),
            PointRef::LineStart(start),
        ));
    }
    for line in [bottom, top] {
        sk.add_constraint(Horizontal(line));
    }
    for line in [right, left] {
        sk.add_constraint(Vertical(line));
    }
    let height_dim = sk.add_constraint(Constraint::distance_y(
        PointRef::LineEnd(right),
        PointRef::LineStart(right),
        height,
    ));
    let width_dim = sk.add_constraint(Constraint::distance_x(
        PointRef::LineEnd(bottom),
        PointRef::Point(origin),
        width,
    ));
    (sk, bottom, width_dim, height_dim)
}

/// A four-bar linkage: `a` and `d` pinned, three driven link lengths.
/// Returns the sketch and the three dimension slots.
fn four_bar(span: f64, crank: f64, coupler: f64, rocker: f64) -> (Sketch, [usize; 3]) {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(0.0, crank);
    let c = sk.add_point(span - 1.0, rocker * 0.9);
    let d = sk.add_point(span, 0.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Point(a), 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(PointRef::Point(d), span, 0.0));
    let mut dims = [0usize; 3];
    for (i, (p, q, length)) in [(a, b, crank), (b, c, coupler), (c, d, rocker)]
        .into_iter()
        .enumerate()
    {
        dims[i] = sk.add_constraint(Constraint::point_distance(
            PointRef::Point(p),
            PointRef::Point(q),
            length,
        ));
    }
    (sk, dims)
}

/// Two spheres tangent on the inside, the smaller inside the larger, on
/// a line through the larger one's centre. Returns the sketch and the
/// big sphere's radius dimension.
fn nested_spheres() -> (SpaceSketch, usize) {
    let mut sk = SpaceSketch::new();
    let big = sk.add_sphere(0.0, 0.0, 0.0, 5.0);
    let small = sk.add_sphere(3.0, 0.0, 0.0, 2.0);
    sk.add_constraint(SpaceConstraint::Fixed {
        point: SpacePointRef::SphereCenter(big),
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    for (dy, dz) in [(1.0, 0.0), (0.0, 1.0)] {
        sk.add_constraint(SpaceConstraint::PointDistanceAlongDir {
            p1: SpacePointRef::SphereCenter(big),
            p2: SpacePointRef::SphereCenter(small),
            dx: 0.0,
            dy,
            dz,
            distance: 0.0,
        });
    }
    sk.add_constraint(SpaceConstraint::SphereTangent {
        s1: big,
        s2: small,
        external: false,
    });
    let big_radius = sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere: big,
        radius: 5.0,
    });
    sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere: small,
        radius: 2.0,
    });
    (sk, big_radius)
}

/// A width driven to the other side of zero is walked through zero —
/// where a rectangle has no side to be on — and arrives, instead of
/// jumping straight to the mirror image and reporting the flip.
#[test]
fn a_jump_through_zero_is_walked_there() {
    let (mut sk, bottom, width, _) = driven_rectangle(4.0, 3.0);
    assert_eq!(sk.solve().status, SolveStatus::Converged, "the fixture");

    assert!(sk.set_constraint_value(width, -4.0));
    let result = sk.solve();
    assert_eq!(
        result.status,
        SolveStatus::Converged,
        "the walk went through zero width, so nothing turned over"
    );
    assert!(result.converged);
    assert!(
        (sk.point_coords(PointRef::LineEnd(bottom)).x + 4.0).abs() < 1e-9,
        "the width the host asked for is the one it is at: {}",
        sk.point_coords(PointRef::LineEnd(bottom)).x
    );
}

/// The same jump without the walk — a localized solve — still reports
/// the flip: the walk is what [`Sketch::solve`] adds, not a change to
/// what a branch is.
#[test]
fn a_localized_solve_does_not_walk() {
    let (mut sk, bottom, width, _) = driven_rectangle(4.0, 3.0);
    assert!(sk.solve().converged, "the fixture");
    assert!(sk.set_constraint_value(width, -4.0));
    let result = sk.solve_selected(None);
    assert_eq!(result.status, SolveStatus::BranchFlip);
    assert!(sk.point_coords(PointRef::LineEnd(bottom)).x < 0.0);
}

/// Two dimensions changed at once are walked together, on one fraction
/// of the whole jump: both values are the host's at the end.
#[test]
fn two_dimensions_walk_together() {
    let (mut sk, bottom, width, height) = driven_rectangle(4.0, 3.0);
    assert!(sk.solve().converged, "the fixture");

    assert!(sk.set_constraint_value(width, -4.0));
    assert!(sk.set_constraint_value(height, 11.0));
    let result = sk.solve();
    assert!(result.converged, "{:?}", result.status);
    assert_eq!(result.status, SolveStatus::Converged);
    let end = sk.point_coords(PointRef::LineEnd(bottom));
    assert!((end.x + 4.0).abs() < 1e-9, "width: {}", end.x);
    assert_eq!(sk.constraint_value(height), Some(11.0));
    assert!(
        (sk.measure(height).unwrap() - 11.0).abs() < 1e-9,
        "height: {:?}",
        sk.measure(height)
    );
}

/// A jump a solve makes in one go is solved exactly as it was before
/// there was any walking: the same geometry, bit for bit, and the same
/// number of iterations as the same change made by replacing the
/// constraint (which no jump is remembered for).
#[test]
fn a_jump_a_solve_makes_is_left_alone() {
    let (mut sk, _, width, _) = driven_rectangle(4.0, 3.0);
    assert!(sk.solve().converged, "the fixture");

    let mut walked = sk.clone();
    assert!(walked.set_constraint_value(width, 6.0));
    let a = walked.solve();

    let mut replaced = sk.clone();
    let mut swapped = replaced.constraints()[width].clone();
    assert!(swapped.set_value(6.0));
    replaced.replace_constraint(width, swapped);
    let b = replaced.solve();

    assert_eq!(a.iterations, b.iterations, "the solve took another path");
    assert_eq!(a.status, b.status);
    assert_eq!(
        walked.params_snapshot(),
        replaced.params_snapshot(),
        "a value set and a value replaced solve to the same geometry"
    );
}

/// A value the host takes back before solving is no jump: the solve is
/// the ordinary one.
#[test]
fn a_value_taken_back_is_no_jump() {
    let (mut sk, _, width, _) = driven_rectangle(4.0, 3.0);
    assert!(sk.solve().converged, "the fixture");
    let before = sk.params_snapshot();

    assert!(sk.set_constraint_value(width, -40.0));
    assert!(sk.set_constraint_value(width, 4.0));
    let result = sk.solve();
    assert_eq!(result.status, SolveStatus::Converged);
    assert_eq!(sk.params_snapshot(), before, "nothing had to move");
}

/// A value nothing can reach: the dimension stays the one the host
/// asked for, the geometry stays at the furthest the walk got to, and
/// the status is the one the direct solve ended with.
#[test]
fn a_jump_out_of_reach_says_how_far_it_got() {
    let (mut sk, dims) = four_bar(6.0, 3.0, 5.0, 3.0);
    assert!(sk.solve().converged, "the fixture");
    let coupler = dims[1];
    let reach = 6.0 + 3.0 + 3.0;

    assert!(sk.set_constraint_value(coupler, 40.0));
    let result = sk.solve();
    assert!(
        !result.converged,
        "40 is past the reach of {reach}: {:?}",
        result.status
    );
    assert_eq!(
        sk.constraint_value(coupler),
        Some(40.0),
        "the library does not put a dimension back"
    );
    let got = sk.measure(coupler).unwrap();
    assert!(
        got > 5.0 + 1e-6 && got < reach + 1e-6,
        "the walk should have got partway: {got}"
    );
    assert!(
        result.residual > 1e-3,
        "the residual is the one of the state left behind, which does not meet 40: {:.3e} at {got}",
        result.residual
    );
}

/// A later solve of a jump that was given up on starts from where the
/// walk stopped, and a value back within reach is met.
#[test]
fn a_jump_given_up_on_can_be_asked_for_again() {
    let (mut sk, dims) = four_bar(6.0, 3.0, 5.0, 3.0);
    assert!(sk.solve().converged, "the fixture");
    let coupler = dims[1];

    assert!(sk.set_constraint_value(coupler, 40.0));
    assert!(!sk.solve().converged);

    assert!(sk.set_constraint_value(coupler, 7.0));
    let result = sk.solve();
    assert!(result.converged, "{:?}", result.status);
    assert!((sk.measure(coupler).unwrap() - 7.0).abs() < 1e-6);
}

/// In 3D: the bigger sphere driven smaller than the one tangent inside
/// it cannot be walked there. The radii can only pass each other through
/// the configuration where the centres coincide, and Newton has nowhere
/// to step from there — so the walk gives up, and what the host gets is
/// the solution the direct solve had found, flip and all, rather than a
/// geometry stuck halfway.
#[test]
fn a_three_d_jump_that_cannot_be_walked_keeps_the_solution_it_found() {
    let (mut sk, big_radius) = nested_spheres();
    assert_eq!(sk.solve().status, SolveStatus::Converged, "the fixture");

    assert!(sk.set_constraint_value(big_radius, 1.0));
    let result = sk.solve();
    assert_eq!(result.status, SolveStatus::BranchFlip);
    assert!(
        result.converged,
        "a flip is still a solution: {:.3e}",
        result.residual
    );
    assert!(
        result.residual < 1e-9,
        "the geometry meets the values: {:.3e}",
        result.residual
    );
    assert_eq!(sk.constraint_value(big_radius), Some(1.0));
    assert!(
        (sk.measure(big_radius).unwrap() - 1.0).abs() < 1e-9,
        "the radius asked for is the one it is at: {:?}",
        sk.measure(big_radius)
    );
}

/// In 3D, a localized solve does not walk either.
#[test]
fn a_three_d_localized_solve_does_not_walk() {
    let (mut sk, big_radius) = nested_spheres();
    assert!(sk.solve().converged, "the fixture");
    assert!(sk.set_constraint_value(big_radius, 1.0));
    assert_eq!(sk.solve_selected(None).status, SolveStatus::BranchFlip);
}
