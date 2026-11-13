//! Editing a 3D sketch in place: removing and suppressing constraints,
//! removing entities, changing dimension values, measuring, checking.
//!
//! The 3D counterpart of `tests/editing.rs`. Indices are stable
//! throughout: nothing an edit does shifts the index of an entity or
//! constraint that is still there.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

use glam::DVec3;
use rcad_gcs::space3d::{
    SpaceConstraint, SpaceEntityKind, SpaceFrame3, SpacePointRef, SpaceSketch, SpaceSketchError,
};
use rcad_gcs::{ConstraintState, SketchError};

/// A fully constrained bracket: a base line from the origin along X and an
/// upright from the end of the base along Z. Returns the sketch, the two
/// lines and the indices of the width (base length) and height (upright
/// length) dimensions.
fn bracket(width: f64, height: f64) -> (SpaceSketch, [usize; 2], usize, usize) {
    let mut sk = SpaceSketch::new();
    let base = sk.add_line(0.0, 0.0, 0.0, width, 0.0, 0.0);
    let upright = sk.add_line(width, 0.0, 0.0, width, 0.0, height);

    sk.add_constraint(SpaceConstraint::fix_point(
        SpacePointRef::LineStart(base),
        0.0,
        0.0,
        0.0,
    ));
    sk.add_constraint(SpaceConstraint::LineParallelToDir {
        line: base,
        dx: 1.0,
        dy: 0.0,
        dz: 0.0,
    });
    sk.add_constraint(SpaceConstraint::coincident(
        SpacePointRef::LineEnd(base),
        SpacePointRef::LineStart(upright),
    ));
    sk.add_constraint(SpaceConstraint::LineParallelToDir {
        line: upright,
        dx: 0.0,
        dy: 0.0,
        dz: 1.0,
    });
    let width_dim = sk.add_constraint(SpaceConstraint::LineLength {
        line: base,
        length: width,
    });
    let height_dim = sk.add_constraint(SpaceConstraint::LineLength {
        line: upright,
        length: height,
    });
    (sk, [base, upright], width_dim, height_dim)
}

/// Move the end of the base, and the upright with it, so that the bracket
/// is `width` wide, keeping the two lines joined.
fn widen_by_hand(sk: &mut SpaceSketch, lines: &[usize; 2], width: f64) {
    let base = sk.entities[lines[0]].param_start;
    let upright = sk.entities[lines[1]].param_start;
    let mut params = sk.params_snapshot();
    params[base + 3] = width; // base: end x
    params[upright] = width; // upright: start x
    params[upright + 3] = width; // upright: end x
    sk.set_params(params);
}

fn line_length(sk: &SpaceSketch, line: usize) -> f64 {
    let (x1, y1, z1, x2, y2, z2) = sk.line_endpoints(line).expect("a line");
    ((x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2)).sqrt()
}

#[test]
fn add_constraint_returns_stable_indices() {
    let mut sk = SpaceSketch::new();
    let a = sk.add_point(0.0, 0.0, 0.0);
    let b = sk.add_point(3.0, 4.0, 0.0);
    let first = sk.add_constraint(SpaceConstraint::fix_point(a, 0.0, 0.0, 0.0));
    let second = sk.add_constraint(SpaceConstraint::point_distance(a, b, 5.0));
    let third = sk.add_constraint(SpaceConstraint::PointDistanceAlongDir {
        p1: b.into(),
        p2: a.into(),
        dx: 1.0,
        dy: 0.0,
        dz: 0.0,
        distance: 3.0,
    });
    assert_eq!((first, second, third), (0, 1, 2));

    assert!(sk.remove_constraint(second));
    // The slots after the removed one keep their indices.
    assert_eq!(sk.constraints.len(), 3);
    assert_eq!(sk.constraint_state(first), ConstraintState::Active);
    assert_eq!(sk.constraint_state(second), ConstraintState::Removed);
    assert_eq!(sk.constraint_state(third), ConstraintState::Active);
    assert!(
        matches!(
            sk.constraints[third],
            SpaceConstraint::PointDistanceAlongDir { .. }
        ),
        "index 2 still names the same constraint"
    );
    // A new constraint gets a new index, not the freed one.
    let fourth = sk.add_constraint(SpaceConstraint::point_distance(a, b, 5.0));
    assert_eq!(fourth, 3);
    assert!(!sk.remove_constraint(second), "already removed");
}

#[test]
fn a_removed_constraint_is_no_longer_enforced_or_diagnosed() {
    let (mut sk, lines, width_dim, _) = bracket(4.0, 3.0);
    assert_eq!(sk.dof(), 0);
    assert_eq!(sk.structural_dof_estimate(), 0);
    assert!(sk.solve().converged);

    assert!(sk.remove_constraint(width_dim));
    assert_eq!(sk.dof(), 1, "the width is free again");
    assert_eq!(sk.structural_dof_estimate(), 1);
    assert_eq!(sk.dof_analysis().equations, 11);
    assert!(
        sk.components().iter().flatten().all(|&c| c != width_dim),
        "a removed constraint belongs to no component"
    );

    // Widen the bracket: nothing pulls it back to 4.
    widen_by_hand(&mut sk, &lines, 6.5);
    let result = sk.solve();
    assert!(result.converged);
    assert_eq!(result.iterations, 0, "the wider bracket already holds");
    assert!((line_length(&sk, lines[0]) - 6.5).abs() < 1e-9);
    assert!(sk.diagnose().redundant.is_empty());
}

#[test]
fn a_suppressed_dimension_is_a_reference_dimension() {
    let (mut sk, lines, width_dim, _) = bracket(4.0, 3.0);
    assert!(sk.solve().converged);

    // Stop driving the width and widen the bracket by hand.
    assert!(sk.set_constraint_enabled(width_dim, false));
    assert_eq!(sk.constraint_state(width_dim), ConstraintState::Suppressed);
    widen_by_hand(&mut sk, &lines, 7.0);
    assert!(sk.solve().converged);

    // The dimension still says 4, but it measures 7 and holds nothing.
    assert_eq!(sk.constraint_value(width_dim), Some(4.0));
    let measured = sk.measure(width_dim).expect("a length can be measured");
    assert!((measured - 7.0).abs() < 1e-9, "measured {measured}");

    // Driving it again pulls the width back.
    assert!(sk.set_constraint_enabled(width_dim, true));
    assert!(sk.solve().converged);
    assert!((line_length(&sk, lines[0]) - 4.0).abs() < 1e-9);

    // A removed constraint cannot be re-enabled.
    sk.remove_constraint(width_dim);
    assert!(!sk.set_constraint_enabled(width_dim, true));
    assert_eq!(sk.constraint_state(width_dim), ConstraintState::Removed);
}

#[test]
fn changing_a_dimension_re_solves_from_where_the_sketch_stands() {
    let (mut sk, lines, width_dim, height_dim) = bracket(4.0, 3.0);
    assert!(sk.solve().converged);
    let original = sk.params_snapshot();
    let revision = sk.topology_revision();

    // Drive the width up and the height down, then both back: a round trip
    // through three warm-started solves, with no sketch rebuilt in between.
    assert!(sk.set_constraint_value(width_dim, 9.0));
    assert!(sk.set_constraint_value(height_dim, 1.5));
    assert!(sk.solve().converged);
    assert!((line_length(&sk, lines[0]) - 9.0).abs() < 1e-9);
    assert!((line_length(&sk, lines[1]) - 1.5).abs() < 1e-9);

    assert!(sk.set_constraint_value(width_dim, 6.0));
    assert!(sk.solve().converged);

    assert!(sk.set_constraint_value(width_dim, 4.0));
    assert!(sk.set_constraint_value(height_dim, 3.0));
    assert!(sk.solve().converged);

    let worst = original
        .iter()
        .zip(sk.params_snapshot())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(worst <= 1e-9, "round-trip drift {worst:.3e}");
    assert_eq!(
        sk.topology_revision(),
        revision,
        "changing values is not a structural change"
    );

    // A relation has no value to set.
    assert!(!sk.set_constraint_value(2, 1.0), "Coincident has no value");
}

#[test]
fn removing_an_entity_takes_its_constraints_with_it() {
    let (mut sk, lines, width_dim, height_dim) = bracket(4.0, 3.0);
    let free_before = sk.dof_analysis().free_params;

    // Remove the upright: the coincidence that joins it to the base, its
    // direction and its length go too; the base stays as it was.
    let gone = sk.remove_entity(lines[1]);
    assert_eq!(gone, vec![2, 3, height_dim]);
    assert!(sk.is_entity_removed(lines[1]));
    for &index in &gone {
        assert_eq!(sk.constraint_state(index), ConstraintState::Removed);
    }
    assert_eq!(sk.constraint_state(width_dim), ConstraintState::Active);

    // Its six parameters are no longer freedoms.
    assert_eq!(sk.dof_analysis().free_params, free_before - 6);
    assert_eq!(sk.dof(), 0);
    assert!(sk.validate().is_ok());
    assert!(sk.solve().converged);
    assert!(sk.remove_entity(lines[1]).is_empty(), "already removed");

    // Entities added afterwards get new ids and work as usual.
    let post = sk.add_line(4.0, 0.0, 0.0, 4.0, 0.2, 2.0);
    assert_eq!(post, 2);
    sk.add_constraint(SpaceConstraint::coincident(
        SpacePointRef::LineEnd(lines[0]),
        SpacePointRef::LineStart(post),
    ));
    sk.add_constraint(SpaceConstraint::LineParallelToDir {
        line: post,
        dx: 0.0,
        dy: 0.0,
        dz: 1.0,
    });
    sk.add_constraint(SpaceConstraint::LineLength {
        line: post,
        length: 2.0,
    });
    assert!(sk.solve().converged);
    assert_eq!(sk.dof(), 0, "the bracket is whole again");
}

/// A removed plane or cylinder takes its representation gauge with it: the
/// freedoms that remain are exactly those of what is left.
#[test]
fn a_removed_plane_or_cylinder_counts_no_freedom() {
    let mut sk = SpaceSketch::new();
    sk.add_point(1.0, 2.0, 3.0);
    let plane = sk.add_plane(0.0, 0.0, 1.0, 2.0);
    let cylinder = sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5);
    assert_eq!(sk.dof(), 3 + 3 + 5, "point, plane, cylinder");

    sk.remove_entity(plane);
    assert_eq!(sk.dof(), 3 + 5);
    sk.remove_entity(cylinder);
    assert_eq!(sk.dof(), 3);
    assert!(sk.solve().converged);
}

#[test]
fn measure_reads_each_kind_of_dimension() {
    let mut sk = SpaceSketch::new();
    let a = sk.add_point(1.0, 1.0, 1.0);
    let b = sk.add_point(4.0, 5.0, 1.0);
    let base = sk.add_line(0.0, 0.0, 0.0, 2.0, 0.0, 0.0);
    let slanted = sk.add_line(0.0, 0.0, 0.0, 0.0, 3.0, 3.0);
    let diagonal = sk.add_line(0.0, 0.0, 0.0, 1.0, 1.0, 0.0);
    // z = 2, given with a normal that is not unit length.
    let level = sk.add_plane(0.0, 0.0, 2.0, 4.0);
    let wall = sk.add_plane(1.0, 0.0, 0.0, 0.0);
    let sphere = sk.add_sphere(0.0, 0.0, 0.0, 2.5);
    let cylinder = sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.75);
    let circle = sk.add_circle(0.0, 0.0, 0.0, 1.25, SpaceFrame3::XY);

    // The values asked for are deliberately not the ones the geometry has.
    let wrong = 99.0;
    let distance = sk.add_constraint(SpaceConstraint::point_distance(a, b, wrong));
    let along = |p1: usize, p2: usize| SpaceConstraint::PointDistanceAlongDir {
        p1: p1.into(),
        p2: p2.into(),
        dx: 0.0,
        dy: 2.0,
        dz: 0.0,
        distance: wrong,
    };
    let up = sk.add_constraint(along(b, a));
    let down = sk.add_constraint(along(a, b));
    let length = sk.add_constraint(SpaceConstraint::LineLength {
        line: slanted,
        length: wrong,
    });
    let right_angle = sk.add_constraint(SpaceConstraint::AngleLineLine {
        l1: base,
        l2: slanted,
        angle_rad: 0.1,
    });
    let half_right = sk.add_constraint(SpaceConstraint::AngleLineLine {
        l1: base,
        l2: diagonal,
        angle_rad: 0.1,
    });
    let planes = sk.add_constraint(SpaceConstraint::AnglePlanePlane {
        p1: level,
        p2: wall,
        angle_rad: 0.1,
    });
    let incline = sk.add_constraint(SpaceConstraint::AngleLinePlane {
        line: slanted,
        plane: level,
        angle_rad: 0.1,
    });
    let height = sk.add_constraint(SpaceConstraint::PointPlaneDistance {
        point: a.into(),
        plane: level,
        distance: wrong,
    });
    let sphere_radius = sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere,
        radius: wrong,
    });
    let cylinder_radius = sk.add_constraint(SpaceConstraint::CylinderRadius {
        cylinder,
        radius: wrong,
    });
    let circle_radius = sk.add_constraint(SpaceConstraint::CircleRadius {
        circle,
        radius: wrong,
    });
    let relation = sk.add_constraint(SpaceConstraint::LineParallelLine(base, diagonal));
    let direction = sk.add_constraint(SpaceConstraint::LineParallelToDir {
        line: base,
        dx: 1.0,
        dy: 0.0,
        dz: 0.0,
    });

    let close = |index: usize, expected: f64| {
        let got = sk.measure(index).expect("measurable");
        assert!(
            (got - expected).abs() < 1e-12,
            "{index}: {got} vs {expected}"
        );
    };
    close(distance, 5.0);
    close(up, 4.0);
    close(down, -4.0);
    close(length, 18.0_f64.sqrt());
    close(right_angle, FRAC_PI_2);
    close(half_right, FRAC_PI_4);
    close(planes, FRAC_PI_2);
    close(incline, FRAC_PI_4);
    close(height, 1.0);
    close(sphere_radius, 2.5);
    close(cylinder_radius, 0.75);
    close(circle_radius, 1.25);
    for index in [relation, direction] {
        assert_eq!(sk.measure(index), None);
        assert_eq!(sk.constraint_value(index), None);
    }
}

#[test]
fn structural_edits_bump_the_topology_revision() {
    let mut sk = SpaceSketch::new();
    let mut last = sk.topology_revision();
    let mut bumped = |sk: &SpaceSketch, what: &str| {
        let now = sk.topology_revision();
        assert!(now > last, "{what} must change the revision");
        last = now;
    };

    let a = sk.add_point(0.0, 0.0, 0.0);
    bumped(&sk, "adding a point");
    let b = sk.add_point(1.0, 0.0, 0.0);
    sk.add_line(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    bumped(&sk, "adding a line");
    sk.add_plane(0.0, 0.0, 1.0, 0.0);
    bumped(&sk, "adding a plane");
    sk.add_sphere(0.0, 0.0, 0.0, 1.0);
    bumped(&sk, "adding a sphere");
    sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    bumped(&sk, "adding a cylinder");
    sk.add_circle(0.0, 0.0, 0.0, 1.0, SpaceFrame3::XY);
    bumped(&sk, "adding a circle");
    sk.add_spline(&[
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (2.0, 1.0, 0.0),
        (3.0, 0.0, 0.0),
    ]);
    bumped(&sk, "adding a spline");

    let c = sk.add_constraint(SpaceConstraint::point_distance(a, b, 1.0));
    bumped(&sk, "adding a constraint");
    sk.set_constraint_enabled(c, false);
    bumped(&sk, "suppressing");
    sk.set_constraint_enabled(c, true);
    bumped(&sk, "enabling");
    sk.replace_constraint(c, SpaceConstraint::coincident(a, b));
    bumped(&sk, "replacing");
    sk.fix_entity(a);
    bumped(&sk, "fixing an entity");
    sk.fix_param(3);
    bumped(&sk, "fixing a parameter");
    sk.remove_constraint(c);
    bumped(&sk, "removing a constraint");
    sk.remove_entity(b);
    bumped(&sk, "removing an entity");

    // Enabling what is already enabled is not a change, and neither is
    // moving geometry.
    let c2 = sk.add_constraint(SpaceConstraint::fix_point(a, 0.0, 0.0, 0.0));
    let before = sk.topology_revision();
    sk.set_constraint_enabled(c2, true);
    sk.drift_free_params(0.1, 7);
    sk.solve();
    assert_eq!(sk.topology_revision(), before);
}

#[test]
fn validate_reports_what_is_wrong_with_the_model() {
    let mut sk = SpaceSketch::new();
    let point = sk.add_point(0.0, 0.0, 0.0);
    let line = sk.add_line(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    let sphere = sk.add_sphere(0.0, 0.0, 0.0, 1.0);
    let circle = sk.add_circle(0.0, 0.0, 0.0, 1.0, SpaceFrame3::XY);
    let arc = sk.add_arc(
        0.0,
        0.0,
        0.0,
        1.0,
        SpaceFrame3::sanitized(DVec3::Z, DVec3::X, 0.0, FRAC_PI_2),
    );
    let spline = sk.add_spline(&[
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (2.0, 1.0, 0.0),
        (3.0, 0.0, 0.0),
    ]);

    // Everything used the way it is meant to be is accepted, including a
    // circle centre read off an arc and a concentric pair of curves.
    sk.add_constraint(SpaceConstraint::fix_point(
        SpacePointRef::CircleCenter(arc),
        0.0,
        0.0,
        0.0,
    ));
    sk.add_constraint(SpaceConstraint::Concentric3D(circle, arc));
    sk.add_constraint(SpaceConstraint::PointAtArcEndpoint {
        point: point.into(),
        arc,
        at_start: true,
    });
    sk.add_constraint(SpaceConstraint::SplineControlPointMidpoint {
        spline,
        midpoint_index: 1,
        first_index: 0,
        second_index: 3,
    });
    assert_eq!(sk.validate(), Ok(()));

    // A sphere radius on a line, and the start point of a sphere.
    let on_line = sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere: line,
        radius: 1.0,
    });
    let bad_ref = sk.add_constraint(SpaceConstraint::fix_point(
        SpacePointRef::LineStart(sphere),
        0.0,
        0.0,
        0.0,
    ));
    // A full circle has no endpoint.
    let no_endpoint = sk.add_constraint(SpaceConstraint::PointAtArcEndpoint {
        point: point.into(),
        arc: circle,
        at_start: false,
    });
    // An entity that was never added.
    let missing = sk.add_constraint(SpaceConstraint::LineLength {
        line: 42,
        length: 1.0,
    });
    // Numbers that are not numbers: a value, and a direction.
    let nan_value = sk.add_constraint(SpaceConstraint::LineLength {
        line,
        length: f64::NAN,
    });
    let nan_direction = sk.add_constraint(SpaceConstraint::LineParallelToDir {
        line,
        dx: 1.0,
        dy: f64::INFINITY,
        dz: 0.0,
    });
    // A control point the spline does not have.
    let no_such_control_point = sk.add_constraint(SpaceConstraint::SplineControlPointMidpoint {
        spline,
        midpoint_index: 1,
        first_index: 0,
        second_index: 7,
    });

    let errors = sk.validate().expect_err("seven defects");
    let expected: Vec<SpaceSketchError> = vec![
        SketchError::WrongEntityKind {
            constraint: on_line,
            entity: line,
            found: SpaceEntityKind::SpaceLine,
            expected: &[SpaceEntityKind::Sphere],
        },
        SketchError::WrongEntityKind {
            constraint: bad_ref,
            entity: sphere,
            found: SpaceEntityKind::Sphere,
            expected: &[SpaceEntityKind::SpaceLine],
        },
        SketchError::WrongEntityKind {
            constraint: no_endpoint,
            entity: circle,
            found: SpaceEntityKind::SpaceCircle,
            expected: &[SpaceEntityKind::SpaceArc],
        },
        SketchError::EntityOutOfRange {
            constraint: missing,
            entity: 42,
        },
        SketchError::NonFiniteValue {
            constraint: nan_value,
        },
        SketchError::NonFiniteValue {
            constraint: nan_direction,
        },
        SketchError::ControlPointOutOfRange {
            constraint: no_such_control_point,
            entity: spline,
            index: 7,
            count: 4,
        },
    ];
    assert_eq!(errors, expected);
    assert!(errors[0].to_string().contains("SpaceLine"));
    assert!(errors[6].to_string().contains("control point 7"));

    // Removing the offenders clears the report.
    for index in [
        on_line,
        bad_ref,
        no_endpoint,
        missing,
        nan_value,
        nan_direction,
        no_such_control_point,
    ] {
        sk.remove_constraint(index);
    }
    assert_eq!(sk.validate(), Ok(()));

    // A constraint added against an entity that is already gone, and a
    // parameter that is not a number.
    let gone = sk.remove_entity(point);
    assert_eq!(gone.len(), 1, "the arc endpoint constraint used the point");
    let stale = sk.add_constraint(SpaceConstraint::fix_point(point, 0.0, 0.0, 0.0));
    let mut params = sk.params_snapshot();
    let poisoned = sk.entities[line].param_start + 4;
    params[poisoned] = f64::NAN;
    // The removed point's own parameters may hold anything.
    params[sk.entities[point].param_start] = f64::NAN;
    sk.set_params(params);
    assert_eq!(
        sk.validate(),
        Err(vec![
            SketchError::NonFiniteParameter { index: poisoned },
            SketchError::EntityRemoved {
                constraint: stale,
                entity: point,
            },
        ])
    );
}

/// Suppressing a constraint takes it out of the decomposition, so an
/// incremental solve that names it has nothing to iterate, and what it
/// reports is the state of the constraints that are enforced.
#[test]
fn suppressed_constraints_stay_out_of_incremental_solves() {
    let mut sk = SpaceSketch::new();
    let a = sk.add_point(0.0, 0.0, 0.0);
    let b = sk.add_point(3.0, 4.0, 0.0);
    let sphere = sk.add_sphere(9.0, 9.0, 9.0, 2.0);
    let pin = sk.add_constraint(SpaceConstraint::fix_point(a, 0.0, 0.0, 0.0));
    let reach = sk.add_constraint(SpaceConstraint::point_distance(a, b, 5.0));
    let radius = sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere,
        radius: 6.0,
    });
    assert_eq!(sk.components(), vec![vec![pin, reach], vec![radius]]);

    assert!(sk.set_constraint_enabled(radius, false));
    assert_eq!(sk.components(), vec![vec![pin, reach]]);

    let before = sk.params_snapshot();
    let result = sk.solve_incremental(&[radius]);
    assert_eq!(result.iterations, 0);
    assert!(
        result.converged,
        "the unmet radius is not enforced: {result:?}"
    );
    assert_eq!(sk.params_snapshot(), before, "nothing moved");
    assert_eq!(sk.measure(radius), Some(2.0));

    // Enforced again, the same call drives it.
    assert!(sk.set_constraint_enabled(radius, true));
    assert!(sk.solve_incremental(&[radius]).converged);
    assert!((sk.measure(radius).unwrap() - 6.0).abs() < 1e-9);
}

#[cfg(feature = "serde")]
#[test]
fn sketches_saved_before_the_edit_state_existed_still_load() {
    // What 0.1.0 wrote: no `states`, `removed` or `revision`.
    let old = r#"{
        "params": [0.0, 0.0, 0.0, 3.0, 4.0, 0.0],
        "entities": [
            {"kind": "SpacePoint", "param_start": 0},
            {"kind": "SpacePoint", "param_start": 3}
        ],
        "constraints": [
            {"Fixed": {"point": {"Point": 0}, "x": 0.0, "y": 0.0, "z": 0.0}},
            {"PointDistance": {"p1": {"Point": 0}, "p2": {"Point": 1}, "distance": 10.0}}
        ],
        "fixed": [false, false, false, false, false, false]
    }"#;
    let mut sk: SpaceSketch = serde_json::from_str(old).expect("old sketches deserialize");
    assert_eq!(sk.constraint_state(0), ConstraintState::Active);
    assert_eq!(sk.constraint_state(1), ConstraintState::Active);
    assert!(!sk.is_entity_removed(1));
    assert!(sk.solve().converged);
    assert!((sk.measure(1).unwrap() - 10.0).abs() < 1e-9);

    // And the edit state round-trips.
    sk.set_constraint_enabled(1, false);
    let json = serde_json::to_string(&sk).expect("serialize");
    let back: SpaceSketch = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(back.constraint_state(1), ConstraintState::Suppressed);
    assert_eq!(back.topology_revision(), sk.topology_revision());
}
