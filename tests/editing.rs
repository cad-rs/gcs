//! Editing a sketch in place: removing and suppressing constraints,
//! removing entities, changing dimension values, measuring, checking.
//!
//! Indices are stable throughout: nothing an edit does shifts the index of
//! an entity or constraint that is still there.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::entity::{EntityKind, PointRef};
use rcad_gcs::{ConstraintState, Sketch, SketchError};

/// A fully constrained rectangle with its lower-left corner pinned at the
/// origin. Returns the sketch, its four lines (bottom, right, top, left) and
/// the indices of the width and height dimensions.
fn rectangle(width: f64, height: f64) -> (Sketch, [usize; 4], usize, usize) {
    let mut sk = Sketch::new();
    let bottom = sk.add_line(0.0, 0.0, width, 0.0);
    let right = sk.add_line(width, 0.0, width, height);
    let top = sk.add_line(width, height, 0.0, height);
    let left = sk.add_line(0.0, height, 0.0, 0.0);

    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(bottom), 0.0, 0.0));
    for (end, start) in [(bottom, right), (right, top), (top, left), (left, bottom)] {
        sk.add_constraint(Constraint::coincident(
            PointRef::LineEnd(end),
            PointRef::LineStart(start),
        ));
    }
    sk.add_constraint(Constraint::Horizontal(bottom));
    sk.add_constraint(Constraint::Vertical(right));
    sk.add_constraint(Constraint::Horizontal(top));
    sk.add_constraint(Constraint::Vertical(left));
    let width_dim = sk.add_constraint(Constraint::LineLength {
        line: bottom,
        length: width,
    });
    let height_dim = sk.add_constraint(Constraint::LineLength {
        line: right,
        length: height,
    });
    (sk, [bottom, right, top, left], width_dim, height_dim)
}

/// Move the right-hand side of the rectangle so that it is `width` wide,
/// keeping every corner joined: the bottom line's end, the right line and
/// the top line's start all move together.
fn widen_by_hand(sk: &mut Sketch, lines: &[usize; 4], width: f64) {
    let mut params = sk.params_snapshot();
    let start_of = |line: usize| sk.entity_param_range(line).start;
    params[start_of(lines[0]) + 2] = width; // bottom: end x
    params[start_of(lines[1])] = width; // right: start x
    params[start_of(lines[1]) + 2] = width; // right: end x
    params[start_of(lines[2])] = width; // top: start x
    sk.set_params(params);
}

fn line_length(sk: &Sketch, line: usize) -> f64 {
    (sk.point_coords(PointRef::LineEnd(line)) - sk.point_coords(PointRef::LineStart(line))).length()
}

#[test]
fn add_constraint_returns_stable_indices() {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(3.0, 4.0);
    let first = sk.add_constraint(Constraint::fix_point(a, 0.0, 0.0));
    let second = sk.add_constraint(Constraint::point_distance(a, b, 5.0));
    let third = sk.add_constraint(Constraint::distance_x(b, a, 3.0));
    assert_eq!((first, second, third), (0, 1, 2));

    assert!(sk.remove_constraint(second));
    // The slots after the removed one keep their indices.
    assert_eq!(sk.constraints_len(), 3);
    assert_eq!(sk.constraint_state(first), ConstraintState::Active);
    assert_eq!(sk.constraint_state(second), ConstraintState::Removed);
    assert_eq!(sk.constraint_state(third), ConstraintState::Active);
    assert!(
        matches!(sk.constraints()[third], Constraint::DistanceX { .. }),
        "index 2 still names the same constraint"
    );
    // A new constraint gets a new index, not the freed one.
    let fourth = sk.add_constraint(Constraint::distance_y(b, a, 4.0));
    assert_eq!(fourth, 3);
    assert!(!sk.remove_constraint(second), "already removed");
}

#[test]
fn a_removed_constraint_is_no_longer_enforced_or_diagnosed() {
    let (mut sk, lines, width_dim, _) = rectangle(4.0, 3.0);
    assert_eq!(sk.dof(), 0);
    assert!(sk.solve().converged);

    assert!(sk.remove_constraint(width_dim));
    assert_eq!(sk.dof(), 1, "the width is free again");
    assert!(
        sk.components().iter().flatten().all(|&c| c != width_dim),
        "a removed constraint belongs to no component"
    );

    // Widen the rectangle: nothing pulls it back to 4.
    widen_by_hand(&mut sk, &lines, 6.5);
    let result = sk.solve();
    assert!(result.converged);
    assert_eq!(result.iterations, 0, "the wider rectangle already holds");
    assert!((line_length(&sk, lines[0]) - 6.5).abs() < 1e-9);
    assert!(sk.diagnose().redundant.is_empty());
}

#[test]
fn a_suppressed_dimension_is_a_reference_dimension() {
    let (mut sk, lines, width_dim, _) = rectangle(4.0, 3.0);
    assert!(sk.solve().converged);

    // Stop driving the width and widen the rectangle by hand.
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
    let (mut sk, lines, width_dim, height_dim) = rectangle(4.0, 3.0);
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
    assert!(!sk.set_constraint_value(5, 1.0), "Horizontal has no value");
}

#[test]
fn removing_an_entity_takes_its_constraints_with_it() {
    let (mut sk, lines, width_dim, height_dim) = rectangle(4.0, 3.0);
    let free_before = sk.dof_analysis().free_params;

    // Remove the top line: the two corner coincidences that touch it and
    // its Horizontal go too; the rest stays.
    let gone = sk.remove_entity(lines[2]);
    assert_eq!(gone.len(), 3, "{gone:?}");
    assert!(sk.is_entity_removed(lines[2]));
    for &index in &gone {
        assert_eq!(sk.constraint_state(index), ConstraintState::Removed);
    }
    assert_eq!(sk.constraint_state(width_dim), ConstraintState::Active);
    assert_eq!(sk.constraint_state(height_dim), ConstraintState::Active);

    // Its four parameters are no longer freedoms.
    assert_eq!(sk.dof_analysis().free_params, free_before - 4);
    assert!(sk.validate().is_ok());
    assert!(sk.solve().converged);
    assert!(sk.remove_entity(lines[2]).is_empty(), "already removed");

    // Entities added afterwards get new ids and work as usual.
    let lid = sk.add_line(4.0, 3.0, 0.0, 3.0);
    assert_eq!(lid, 4);
    sk.add_constraint(Constraint::coincident(
        PointRef::LineStart(lid),
        PointRef::LineEnd(lines[1]),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(lid),
        PointRef::LineStart(lines[3]),
    ));
    sk.add_constraint(Constraint::Horizontal(lid));
    assert!(sk.solve().converged);
    assert_eq!(sk.dof(), 0, "the rectangle is closed again");
}

#[test]
fn measure_reads_each_kind_of_dimension() {
    let mut sk = Sketch::new();
    let a = sk.add_point(1.0, 1.0);
    let b = sk.add_point(4.0, 5.0);
    let base = sk.add_line(0.0, 0.0, 2.0, 0.0);
    let slanted = sk.add_line(0.0, 0.0, 0.0, 3.0);
    let circle = sk.add_circle(0.0, 0.0, 2.5);

    // The values asked for are deliberately not the ones the geometry has.
    let distance = sk.add_constraint(Constraint::point_distance(a, b, 99.0));
    let dx = sk.add_constraint(Constraint::distance_x(b, a, 99.0));
    let dy = sk.add_constraint(Constraint::distance_y(a, b, 99.0));
    let length = sk.add_constraint(Constraint::LineLength {
        line: slanted,
        length: 99.0,
    });
    let angle = sk.add_constraint(Constraint::Angle {
        l1: base,
        l2: slanted,
        angle_rad: 0.1,
    });
    let radius = sk.add_constraint(Constraint::Radius {
        circle,
        radius: 99.0,
    });
    let diameter = sk.add_constraint(Constraint::Diameter {
        circle,
        diameter: 99.0,
    });
    let relation = sk.add_constraint(Constraint::Horizontal(base));

    let close = |index: usize, expected: f64| {
        let got = sk.measure(index).expect("measurable");
        assert!(
            (got - expected).abs() < 1e-12,
            "{index}: {got} vs {expected}"
        );
    };
    close(distance, 5.0);
    close(dx, 3.0);
    close(dy, -4.0);
    close(length, 3.0);
    close(angle, std::f64::consts::FRAC_PI_2);
    close(radius, 2.5);
    close(diameter, 5.0);
    assert_eq!(sk.measure(relation), None);
    assert_eq!(sk.constraint_value(relation), None);
}

#[test]
fn structural_edits_bump_the_topology_revision() {
    let mut sk = Sketch::new();
    let mut last = sk.topology_revision();
    let mut bumped = |sk: &Sketch, what: &str| {
        let now = sk.topology_revision();
        assert!(now > last, "{what} must change the revision");
        last = now;
    };

    let a = sk.add_point(0.0, 0.0);
    bumped(&sk, "adding an entity");
    let b = sk.add_point(1.0, 0.0);
    let c = sk.add_constraint(Constraint::point_distance(a, b, 1.0));
    bumped(&sk, "adding a constraint");
    sk.set_constraint_enabled(c, false);
    bumped(&sk, "suppressing");
    sk.set_constraint_enabled(c, true);
    bumped(&sk, "enabling");
    sk.replace_constraint(c, Constraint::distance_x(b, a, 1.0));
    bumped(&sk, "replacing");
    sk.fix_entity(a);
    bumped(&sk, "fixing");
    sk.remove_constraint(c);
    bumped(&sk, "removing a constraint");
    sk.remove_entity(b);
    bumped(&sk, "removing an entity");

    // Enabling what is already enabled is not a change.
    let c2 = sk.add_constraint(Constraint::fix_point(a, 0.0, 0.0));
    let before = sk.topology_revision();
    sk.set_constraint_enabled(c2, true);
    assert_eq!(sk.topology_revision(), before);
}

#[test]
fn validate_reports_what_is_wrong_with_the_model() {
    let mut sk = Sketch::new();
    let point = sk.add_point(0.0, 0.0);
    let line = sk.add_line(0.0, 0.0, 1.0, 0.0);
    let circle = sk.add_circle(0.0, 0.0, 1.0);
    assert!(sk.validate().is_ok());

    // A radius on a line, and the start point of a circle.
    let on_line = sk.add_constraint(Constraint::Radius {
        circle: line,
        radius: 1.0,
    });
    let bad_ref = sk.add_constraint(Constraint::fix_point(PointRef::LineStart(circle), 0.0, 0.0));
    // An entity that was never added.
    let missing = sk.add_constraint(Constraint::Horizontal(42));
    // A value that is not a number.
    let nan = sk.add_constraint(Constraint::LineLength {
        line,
        length: f64::NAN,
    });

    let errors = sk.validate().expect_err("four defects");
    assert!(errors.contains(&SketchError::WrongEntityKind {
        constraint: on_line,
        entity: line,
        found: EntityKind::Line,
        expected: &[EntityKind::Circle, EntityKind::Arc],
    }));
    assert!(errors.contains(&SketchError::WrongEntityKind {
        constraint: bad_ref,
        entity: circle,
        found: EntityKind::Circle,
        expected: &[EntityKind::Line],
    }));
    assert!(errors.contains(&SketchError::EntityOutOfRange {
        constraint: missing,
        entity: 42,
    }));
    assert!(errors.contains(&SketchError::NonFiniteValue { constraint: nan }));
    assert_eq!(errors.len(), 4, "{errors:?}");
    assert!(errors[0].to_string().contains("constraint"));

    // Removing the offenders clears the report.
    for index in [on_line, bad_ref, missing, nan] {
        sk.remove_constraint(index);
    }
    assert!(sk.validate().is_ok());

    // A constraint added against an entity that is already gone.
    sk.remove_entity(point);
    let stale = sk.add_constraint(Constraint::fix_point(point, 0.0, 0.0));
    assert_eq!(
        sk.validate(),
        Err(vec![SketchError::EntityRemoved {
            constraint: stale,
            entity: point,
        }])
    );
}

/// A `Fixed` has no single value, but its target is checked all the same.
#[test]
fn validate_reports_a_fixed_target_that_is_not_finite() {
    let mut sk = Sketch::new();
    let point = sk.add_point(0.0, 0.0);
    let good = sk.add_constraint(Constraint::fix_point(point, 1.0, 2.0));
    let bad = sk.add_constraint(Constraint::fix_point(point, 1.0, f64::NAN));
    assert!(sk.constraints()[good].has_finite_values());
    assert_eq!(
        sk.validate(),
        Err(vec![SketchError::NonFiniteValue { constraint: bad }])
    );
}

#[cfg(feature = "serde")]
#[test]
fn sketches_saved_before_the_edit_state_existed_still_load() {
    // What 0.1.0 wrote: no `states`, `removed` or `revision`.
    let old = r#"{
        "params": [0.0, 0.0, 3.0, 4.0],
        "fixed": [false, false, false, false],
        "entities": [
            {"kind": "Point", "param_start": 0},
            {"kind": "Point", "param_start": 2}
        ],
        "constraints": [
            {"Fixed": {"point": {"Point": 0}, "x": 0.0, "y": 0.0}},
            {"PointDistance": {"p1": {"Point": 0}, "p2": {"Point": 1}, "distance": 10.0}}
        ]
    }"#;
    let mut sk: Sketch = serde_json::from_str(old).expect("old sketches deserialize");
    assert_eq!(sk.constraint_state(0), ConstraintState::Active);
    assert_eq!(sk.constraint_state(1), ConstraintState::Active);
    assert!(!sk.is_entity_removed(1));
    assert!(sk.solve().converged);
    assert!((sk.measure(1).unwrap() - 10.0).abs() < 1e-9);

    // And the edit state round-trips.
    sk.set_constraint_enabled(1, false);
    let json = serde_json::to_string(&sk).expect("serialize");
    let back: Sketch = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(back.constraint_state(1), ConstraintState::Suppressed);
    assert_eq!(back.topology_revision(), sk.topology_revision());
}
