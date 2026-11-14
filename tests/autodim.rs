//! Automatic dimensioning (P4.6). The corpus gate
//! (`auto_dimensions_define_every_case` in `tests/corpus.rs`) defines
//! every solvable case of the 2D corpus and checks the under-constrained
//! ones against what their construction lacks; these cover free geometry,
//! a given datum and the 3D sketch.

use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::space3d::{SpaceConstraint, SpacePointRef, SpaceSketch};
use rcad_gcs::{Constraint, ConstraintStatus, Sketch};

/// Adds `suggestions` to a copy of `sk`; the copy's dof and whether every
/// added constraint is met and independent.
fn define(sk: &Sketch, suggestions: &[Constraint]) -> (i64, bool) {
    let mut defined = sk.clone();
    let first = defined.constraints_len();
    for c in suggestions {
        defined.add_constraint(c.clone());
    }
    let status = defined.status();
    let clean =
        (first..defined.constraints_len()).all(|c| status.constraints[c] == ConstraintStatus::Met);
    (defined.diagnose().dof, clean)
}

#[test]
fn free_geometry_is_pinned_at_its_first_point_and_dimensioned_from_it() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 1.0);
    let l = sk.add_line(2.0, 0.0, 5.0, 4.0);
    let c = sk.add_circle(-1.0, 3.0, 1.5);
    let a = sk.add_arc(4.0, -2.0, 2.0, 0.3, 2.0);
    let suggestions = sk.auto_dimensions(None);
    // Fixed point, then length + 4 end coordinates − 1 (the length fixes
    // one), radius + centre, radius + centre + the two ends' angles.
    let equations: usize = suggestions.iter().map(|c| c.equation_count()).sum();
    assert_eq!(equations as i64, sk.diagnose().dof);
    assert!(
        matches!(suggestions[0], Fixed { point: PointRef::Point(q), .. } if q == p),
        "{suggestions:?}"
    );
    assert!(suggestions
        .iter()
        .any(|s| matches!(s, LineLength { line, .. } if *line == l)));
    assert!(suggestions
        .iter()
        .any(|s| matches!(s, Radius { circle, .. } if *circle == c)));
    assert!(suggestions
        .iter()
        .any(|s| matches!(s, Radius { circle, .. } if *circle == a)));
    assert_eq!(define(&sk, &suggestions), (0, true));
}

#[test]
fn a_given_datum_is_measured_from_and_pinned_if_it_can_move() {
    let mut sk = Sketch::new();
    let q = sk.add_point(0.0, 0.0);
    let p = sk.add_point(3.0, 2.0);
    sk.add_constraint(Constraint::point_distance(q, p, 13f64.sqrt()));
    assert!(sk.solve().converged);
    let suggestions = sk.auto_dimensions(Some(PointRef::Point(p)));
    assert!(
        matches!(suggestions[0], Fixed { point: PointRef::Point(x), .. } if x == p),
        "{suggestions:?}"
    );
    // Positions are measured from the datum.
    assert!(suggestions.iter().skip(1).all(|s| matches!(
        s,
        DistanceX { p2: PointRef::Point(x), .. } | DistanceY { p2: PointRef::Point(x), .. } if *x == p
    )));
    assert_eq!(define(&sk, &suggestions), (0, true));
}

#[test]
fn a_fully_defined_sketch_is_offered_nothing() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 2.0);
    sk.add_constraint(Constraint::fix_point(p, 1.0, 2.0));
    assert!(sk.auto_dimensions(None).is_empty());
}

#[test]
fn free_space_geometry_is_dimensioned_but_a_cylinder_axis_is_left() {
    let mut sk = SpaceSketch::new();
    let p = sk.add_point(1.0, 2.0, 3.0);
    sk.add_line(0.0, 0.0, 0.0, 1.0, 2.0, 2.0);
    sk.add_sphere(5.0, 5.0, 5.0, 1.0);
    sk.add_plane(0.3, -0.4, 0.5, 1.2);
    let suggestions = sk.auto_dimensions(None);
    assert!(
        matches!(suggestions[0], SpaceConstraint::Fixed { point: SpacePointRef::Point(x), .. } if x == p),
        "{suggestions:?}"
    );
    let equations: usize = suggestions.iter().map(|c| c.equation_count()).sum();
    assert_eq!(equations as i64, sk.diagnose().dof);
    let mut defined = sk.clone();
    for c in &suggestions {
        defined.add_constraint(c.clone());
    }
    assert_eq!(defined.diagnose().dof, 0);
    assert!(defined.diagnose().redundant.is_empty());

    // A cylinder: its radius is dimensioned, its axis (four freedoms) is
    // not — nothing in the vocabulary dimensions it.
    let mut sk = SpaceSketch::new();
    sk.add_point(0.0, 0.0, 0.0);
    sk.add_cylinder(1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.5);
    let suggestions = sk.auto_dimensions(None);
    let mut defined = sk.clone();
    for c in &suggestions {
        defined.add_constraint(c.clone());
    }
    assert!(suggestions
        .iter()
        .any(|s| matches!(s, SpaceConstraint::CylinderRadius { .. })));
    assert_eq!(defined.diagnose().dof, 4);
    assert!(defined.diagnose().redundant.is_empty());
}

#[test]
fn a_datum_held_in_one_coordinate_is_measured_from_a_held_point() {
    // The datum's x is held (3 from a fixed point), its y is not: a Fixed
    // on it would repeat the x, so it is not offered; the fixed point's y
    // distance from the datum measures the datum's y.
    let mut sk = Sketch::new();
    let held = sk.add_point(0.0, 0.0);
    let datum = sk.add_point(3.0, 1.0);
    let other = sk.add_point(4.0, 4.0);
    sk.add_constraint(Constraint::fix_point(held, 0.0, 0.0));
    sk.add_constraint(DistanceX {
        p1: PointRef::Point(datum),
        p2: PointRef::Point(held),
        distance: 3.0,
    });
    let suggestions = sk.auto_dimensions(Some(PointRef::Point(datum)));
    assert!(
        !suggestions.iter().any(|s| matches!(s, Fixed { .. })),
        "{suggestions:?}"
    );
    assert!(suggestions.iter().any(|s| matches!(
        s,
        DistanceY { p1: PointRef::Point(h), p2: PointRef::Point(d), .. } if *d == datum && *h == held
    )));
    assert!(suggestions.iter().any(|s| matches!(
        s,
        DistanceX { p1: PointRef::Point(o), .. } if *o == other
    )));
    assert_eq!(define(&sk, &suggestions), (0, true));
}

#[test]
fn the_default_datum_is_the_first_point_that_cannot_move() {
    let mut sk = Sketch::new();
    let free = sk.add_point(5.0, 5.0);
    let pinned = sk.add_point(1.0, 2.0);
    sk.add_constraint(Constraint::fix_point(pinned, 1.0, 2.0));
    let suggestions = sk.auto_dimensions(None);
    // Nothing to pin: the free point is measured from the pinned one.
    assert_eq!(suggestions.len(), 2, "{suggestions:?}");
    for s in &suggestions {
        assert!(
            matches!(s, DistanceX { p1: PointRef::Point(a), p2: PointRef::Point(b), .. }
                | DistanceY { p1: PointRef::Point(a), p2: PointRef::Point(b), .. }
                if *a == free && *b == pinned),
            "{suggestions:?}"
        );
    }
    assert_eq!(define(&sk, &suggestions), (0, true));
}
