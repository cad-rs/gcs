//! Which branch a solve came out on (P5.3): `SolveStatus::BranchFlip`.
//!
//! A drag refuses to cross over (see `tests/drag.rs`); a solve keeps the
//! solution it found and says where it ended up, which is what a host
//! needs to offer an undo or to subdivide the change into smaller ones.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::entity::SpacePointRef;
use rcad_gcs::space3d::SpaceSketch;
use rcad_gcs::{Sketch, SolveStatus};

/// A rectangle with horizontal and vertical sides, its width driven by a
/// signed x distance, so that asking for the other sign asks for the
/// rectangle inside out. Returns the sketch, the bottom line and the
/// index of the width dimension.
fn driven_rectangle(width: f64) -> (Sketch, usize, usize) {
    let mut sk = Sketch::new();
    let origin = sk.add_point(0.0, 0.0);
    let bottom = sk.add_line(0.0, 0.0, width, 0.0);
    let right = sk.add_line(width, 0.0, width, 3.0);
    let top = sk.add_line(width, 3.0, 0.0, 3.0);
    let left = sk.add_line(0.0, 3.0, 0.0, 0.0);
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
    sk.add_constraint(Constraint::distance_y(
        PointRef::LineStart(right),
        PointRef::LineEnd(right),
        -3.0,
    ));
    let dimension = sk.add_constraint(Constraint::distance_x(
        PointRef::LineEnd(bottom),
        PointRef::Point(origin),
        width,
    ));
    (sk, bottom, dimension)
}

/// A dimension changed to the other side of zero mirrors the sketch, and
/// the solve says so while keeping what it found.
///
/// A localized solve is the one that goes straight there: a full
/// [`Sketch::solve`] walks a jump it cannot make in one go, and this one
/// it can walk — through zero width, where a rectangle has no side to be
/// on (see `tests/jump.rs`).
#[test]
fn a_dimension_change_that_mirrors_says_so() {
    let (mut sk, bottom, dimension) = driven_rectangle(4.0);
    let first = sk.solve();
    assert_eq!(first.status, SolveStatus::Converged, "the fixture");
    assert!(sk.point_coords(PointRef::LineEnd(bottom)).x > 0.0);

    assert!(sk.set_constraint_value(dimension, -4.0));
    let result = sk.solve_selected(None);
    assert!(
        result.converged,
        "the mirrored rectangle is a solution: {:?} {:.3e}",
        result.status, result.residual
    );
    assert_eq!(
        result.status,
        SolveStatus::BranchFlip,
        "the solve came out mirrored and did not say so"
    );
    assert!(
        sk.point_coords(PointRef::LineEnd(bottom)).x < 0.0,
        "the solution it found is the one it kept"
    );
}

/// A dimension change that keeps the sketch the way round it was says
/// nothing.
#[test]
fn a_dimension_change_that_keeps_its_side_is_silent() {
    let (mut sk, bottom, dimension) = driven_rectangle(4.0);
    assert!(sk.solve().converged, "the fixture");
    for width in [6.0, 40.0, 0.5] {
        assert!(sk.set_constraint_value(dimension, width));
        let result = sk.solve();
        assert_eq!(
            result.status,
            SolveStatus::Converged,
            "width {width} is the same way round"
        );
        assert!((sk.point_coords(PointRef::LineEnd(bottom)).x - width).abs() < 1e-9);
    }
}

/// In 3D: two spheres tangent on the inside, the smaller one inside the
/// larger. Driving the radii past each other turns that round — the
/// residual `|r1 − r2|` holds either way — and the solve says so.
#[test]
fn spheres_that_swap_which_is_inside_say_so() {
    let mut sk = SpaceSketch::new();
    let big = sk.add_sphere(0.0, 0.0, 0.0, 5.0);
    let small = sk.add_sphere(3.0, 0.0, 0.0, 2.0);
    sk.add_constraint(SpaceConstraint::Fixed {
        point: SpacePointRef::SphereCenter(big),
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    sk.add_constraint(SpaceConstraint::PointDistanceAlongDir {
        p1: SpacePointRef::SphereCenter(big),
        p2: SpacePointRef::SphereCenter(small),
        dx: 0.0,
        dy: 1.0,
        dz: 0.0,
        distance: 0.0,
    });
    sk.add_constraint(SpaceConstraint::PointDistanceAlongDir {
        p1: SpacePointRef::SphereCenter(big),
        p2: SpacePointRef::SphereCenter(small),
        dx: 0.0,
        dy: 0.0,
        dz: 1.0,
        distance: 0.0,
    });
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
    let first = sk.solve();
    assert!(first.converged, "the fixture: {:?}", first.status);
    assert_eq!(first.status, SolveStatus::Converged);

    // The big sphere becomes the small one.
    assert!(sk.set_constraint_value(big_radius, 1.0));
    let result = sk.solve();
    assert!(
        result.converged,
        "it is still a solution: {:?} {:.3e}",
        result.status, result.residual
    );
    assert_eq!(
        result.status,
        SolveStatus::BranchFlip,
        "which sphere is inside turned round without a word"
    );
}

/// In an assembly: a `Parallel` mate holds either way round, so a body
/// that the rest of its mates turn half a turn comes out on the other
/// side of it. The solve keeps the poses it found and says so.
#[test]
fn a_body_assembled_the_other_way_round_says_so() {
    use glam::{DQuat, DVec3};
    use rcad_gcs::rigid3d::{Constraint as Mate, ConstraintSolver, ConstraintSystem};

    let mut system = ConstraintSystem::with_components(2);
    system.add_constraint(Mate::fixed_with_orientation(
        0,
        0,
        DVec3::ZERO,
        DQuat::IDENTITY,
    ));
    // Two coincident pairs pin the body's position and all but one turn;
    // the parallel mate is the one that holds either way round.
    system.add_constraint(Mate::coincident(
        1,
        0,
        DVec3::new(1.0, 0.0, 0.0),
        1,
        DVec3::new(1.0, 0.0, 0.0),
    ));
    system.add_constraint(Mate::coincident(
        2,
        0,
        DVec3::new(0.0, 1.0, 0.0),
        1,
        DVec3::new(0.0, 1.0, 0.0),
    ));
    system.add_constraint(Mate::parallel(
        3,
        0,
        DVec3::new(1.0, 0.0, 0.0),
        1,
        DVec3::new(1.0, 0.0, 0.0),
    ));
    let solver = ConstraintSolver::default();

    // From the identity pose: already assembled, nothing turns round.
    let straight = solver.solve_from(
        &system,
        vec![DVec3::ZERO, DVec3::ZERO],
        vec![DQuat::IDENTITY, DQuat::IDENTITY],
    );
    assert!(straight.success, "the fixture: {:?}", straight.status);
    assert_eq!(straight.status, SolveStatus::Converged);

    // From five sixths of a turn about Z, where the body own x axis
    // points the other way: the coincidences turn it back, and the
    // parallel mate side turns over with it. (Exactly half a turn is a
    // singular start for the coincidences and the solve stalls there,
    // which is a different story -- the drag guard is what keeps a
    // configuration off such a place.)
    let turned = solver.solve_from(
        &system,
        vec![DVec3::ZERO, DVec3::ZERO],
        vec![
            DQuat::IDENTITY,
            DQuat::from_axis_angle(DVec3::Z, 5.0 * std::f64::consts::PI / 6.0),
        ],
    );
    assert!(
        turned.success,
        "it is still a solution: {:?} {:.3e}",
        turned.status, turned.error
    );
    assert_eq!(
        turned.status,
        SolveStatus::BranchFlip,
        "the body was turned round without a word"
    );
}
