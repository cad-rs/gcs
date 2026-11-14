//! How an entity can still move, in geometric terms (P4.2). The corpus
//! gates (`entity_freedoms_match_the_labels` in `tests/corpus.rs`,
//! `body_freedoms_match_the_labels` in `tests/assembly_corpus.rs`) check
//! every entity and body of those corpora against labels; these cover the
//! kinds and shapes of motion the corpora do not reach.

use glam::{DVec2, DVec3};
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::space3d::{SpaceConstraint, SpaceSketch};
use rcad_gcs::{EntityFreedom, EntityStatus, Measure, Motion2, Motion3, Sketch};

const TOL: f64 = 1e-6;

fn freedom(sk: &Sketch, e: usize) -> EntityFreedom<Motion2> {
    sk.freedoms()[e].clone().expect("not removed")
}

fn space_freedom(sk: &SpaceSketch, e: usize) -> EntityFreedom<Motion3> {
    sk.freedoms()[e].clone().expect("not removed")
}

/// A slide along `want` (either sign).
fn slides_along(m: &Motion2, want: DVec2) -> bool {
    matches!(m, Motion2::Translate { direction }
        if direction.dot(want.normalize()).abs() > 1.0 - TOL)
}

fn turns_about(m: &Motion2, want: DVec2) -> bool {
    matches!(m, Motion2::Rotate { center } if (*center - want).length() < TOL)
}

fn resizes(m: &Motion2, want: Measure) -> bool {
    matches!(m, Motion2::Resize { what } if *what == want)
}

#[test]
fn a_free_point_slides_along_both_axes() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 2.0);
    let f = freedom(&sk, p);
    assert_eq!(f.motions.len(), 2, "{f:?}");
    assert!(slides_along(&f.motions[0], DVec2::X), "{f:?}");
    assert!(slides_along(&f.motions[1], DVec2::Y), "{f:?}");
    assert_eq!(f.motions, f.own);
}

#[test]
fn free_lines_and_arcs_slide_turn_and_change_shape() {
    let mut sk = Sketch::new();
    let line = sk.add_line(0.0, 0.0, 3.0, 4.0);
    let arc = sk.add_arc(1.0, 2.0, 3.0, 0.2, 1.4);
    let f = freedom(&sk, line);
    assert_eq!(f.motions.len(), 4, "{f:?}");
    assert!(slides_along(&f.motions[0], DVec2::X), "{f:?}");
    assert!(slides_along(&f.motions[1], DVec2::Y), "{f:?}");
    assert!(turns_about(&f.motions[2], DVec2::new(1.5, 2.0)), "{f:?}");
    assert!(resizes(&f.motions[3], Measure::Length), "{f:?}");
    let f = freedom(&sk, arc);
    assert_eq!(f.motions.len(), 5, "{f:?}");
    assert!(slides_along(&f.motions[0], DVec2::X), "{f:?}");
    assert!(slides_along(&f.motions[1], DVec2::Y), "{f:?}");
    assert!(turns_about(&f.motions[2], DVec2::new(1.0, 2.0)), "{f:?}");
    assert!(resizes(&f.motions[3], Measure::Radius), "{f:?}");
    assert!(resizes(&f.motions[4], Measure::Sweep), "{f:?}");
}

#[test]
fn a_line_held_at_one_end_turns_about_it_and_stretches() {
    let mut sk = Sketch::new();
    let line = sk.add_line(1.0, 1.0, 4.0, 5.0);
    sk.add_constraint(Fixed {
        point: PointRef::LineStart(line),
        x: 1.0,
        y: 1.0,
    });
    let f = freedom(&sk, line);
    assert_eq!(f.motions.len(), 2, "{f:?}");
    assert!(turns_about(&f.motions[0], DVec2::new(1.0, 1.0)), "{f:?}");
    assert!(resizes(&f.motions[1], Measure::Length), "{f:?}");
    assert_eq!(f.motions, f.own);

    // Kept horizontal, it can only stretch.
    sk.add_constraint(Horizontal(line));
    assert!(sk.solve().converged);
    let f = freedom(&sk, line);
    assert_eq!(f.motions.len(), 1, "{f:?}");
    assert!(resizes(&f.motions[0], Measure::Length), "{f:?}");
}

#[test]
fn a_circle_whose_centre_rides_a_line_slides_along_it_and_resizes() {
    let mut sk = Sketch::new();
    let line = sk.add_line(0.0, 0.0, 2.0, 1.0);
    sk.fix_entity(line);
    let c = sk.add_circle(4.0, 2.0, 1.0);
    sk.add_constraint(PointOnLine {
        point: PointRef::Center(c),
        line,
    });
    let f = freedom(&sk, c);
    assert_eq!(f.motions.len(), 2, "{f:?}");
    assert!(slides_along(&f.motions[0], DVec2::new(2.0, 1.0)), "{f:?}");
    assert!(resizes(&f.motions[1], Measure::Radius), "{f:?}");
    assert_eq!(freedom(&sk, line), EntityFreedom::default());
}

#[test]
fn an_arc_through_two_fixed_ends_changes_radius_and_sweep_together() {
    // Both ends held: the centre moves along the bisector while radius and
    // sweep change together — one motion that changes no single measure.
    let mut sk = Sketch::new();
    let arc = sk.add_arc(0.0, 0.0, 2.0, 0.3, 1.9);
    let start = sk.point_coords(PointRef::ArcStart(arc));
    let end = sk.point_coords(PointRef::ArcEnd(arc));
    sk.add_constraint(Fixed {
        point: PointRef::ArcStart(arc),
        x: start.x,
        y: start.y,
    });
    sk.add_constraint(Fixed {
        point: PointRef::ArcEnd(arc),
        x: end.x,
        y: end.y,
    });
    let f = freedom(&sk, arc);
    assert_eq!(f.motions.len(), 1, "{f:?}");
    let Motion2::Other { velocity } = &f.motions[0] else {
        panic!("{f:?}");
    };
    assert_eq!(velocity.len(), 5);
    let norm = velocity.iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!((norm - 1.0).abs() < 1e-9, "{velocity:?}");
    // Radius and sweep both change, in opposite senses.
    let sweep = velocity[4] - velocity[3];
    assert!(velocity[2].abs() > TOL && sweep.abs() > TOL, "{velocity:?}");
    assert!(velocity[2] * sweep < 0.0, "{velocity:?}");
}

#[test]
fn a_bent_chain_moves_through_its_middle_point_but_not_that_point_alone() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    sk.add_constraint(rcad_gcs::Constraint::fix_point(p0, 0.0, 0.0));
    sk.add_constraint(rcad_gcs::Constraint::point_distance(p0, p1, 1.0));
    sk.add_constraint(rcad_gcs::Constraint::point_distance(p1, p2, 1.0));
    let status = sk.status().entities;
    // The middle point swings about the anchor, carrying the end with it;
    // with the anchor and the end held it cannot move.
    let f = freedom(&sk, p1);
    assert_eq!(status[p1], Some(EntityStatus::UnderDefined));
    assert_eq!(f.motions.len(), 1, "{f:?}");
    assert!(slides_along(&f.motions[0], DVec2::Y), "{f:?}");
    assert!(f.own.is_empty(), "{f:?}");
    // The end goes anywhere nearby, or swings about the middle point.
    let f = freedom(&sk, p2);
    assert_eq!(f.motions.len(), 2, "{f:?}");
    assert_eq!(f.own.len(), 1, "{f:?}");
    assert!(slides_along(&f.own[0], DVec2::X), "{f:?}");
    assert_eq!(freedom(&sk, p0), EntityFreedom::default());
}

#[test]
fn equal_radii_move_together_and_neither_alone() {
    let mut sk = Sketch::new();
    let a = sk.add_circle(0.0, 0.0, 1.0);
    let b = sk.add_circle(5.0, 0.0, 1.0);
    for c in [a, b] {
        sk.add_constraint(Fixed {
            point: PointRef::Center(c),
            x: if c == a { 0.0 } else { 5.0 },
            y: 0.0,
        });
    }
    sk.add_constraint(EqualRadius(a, b));
    for c in [a, b] {
        let f = freedom(&sk, c);
        assert_eq!(f.motions.len(), 1, "{f:?}");
        assert!(resizes(&f.motions[0], Measure::Radius), "{f:?}");
        assert!(f.own.is_empty(), "{f:?}");
    }
}

#[test]
fn a_removed_entity_has_no_freedom_and_a_held_one_an_empty_one() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 1.0);
    let q = sk.add_point(2.0, 2.0);
    sk.add_constraint(rcad_gcs::Constraint::fix_point(q, 2.0, 2.0));
    sk.remove_entity(p);
    let freedoms = sk.freedoms();
    assert_eq!(freedoms[p], None);
    assert_eq!(freedoms[q], Some(EntityFreedom::default()));
}

// ── 3D ──────────────────────────────────────────────────────────────────────

fn slides3(m: &Motion3, want: DVec3) -> bool {
    matches!(m, Motion3::Translate { direction }
        if direction.dot(want.normalize()).abs() > 1.0 - TOL)
}

/// A pure turn whose axis is perpendicular to `normal_to`.
fn turns_across(m: &Motion3, normal_to: DVec3) -> bool {
    matches!(m, Motion3::Rotate { axis, .. } if axis.dot(normal_to.normalize()).abs() < TOL)
}

fn resizes3(m: &Motion3, want: Measure) -> bool {
    matches!(m, Motion3::Resize { what } if *what == want)
}

#[test]
fn free_space_entities_move_as_their_kind_allows() {
    let mut sk = SpaceSketch::new();
    let point = sk.add_point(1.0, 2.0, 3.0);
    let line = sk.add_line(0.0, 0.0, 0.0, 1.0, 2.0, 2.0);
    let sphere = sk.add_sphere(1.0, 1.0, 1.0, 2.0);
    let plane = sk.add_plane(0.3, -0.4, 0.5, 1.2);
    let cylinder = sk.add_cylinder(1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.5);
    let spline = sk.add_spline(&[(0.0, 0.0, 0.0), (1.0, 1.0, 0.0), (2.0, 0.0, 1.0)]);

    let f = space_freedom(&sk, point);
    assert_eq!(f.motions.len(), 3, "{f:?}");
    for (m, axis) in f.motions.iter().zip([DVec3::X, DVec3::Y, DVec3::Z]) {
        assert!(slides3(m, axis), "{f:?}");
    }

    // Three slides, two turns across the line about its midpoint (the spin
    // about itself is no motion), and its length.
    let f = space_freedom(&sk, line);
    let along = DVec3::new(1.0, 2.0, 2.0);
    assert_eq!(f.motions.len(), 6, "{f:?}");
    for m in &f.motions[3..5] {
        assert!(turns_across(m, along), "{f:?}");
        let Motion3::Rotate { point, .. } = m else {
            unreachable!()
        };
        assert!((*point - along * 0.5).length() < TOL, "{f:?}");
    }
    assert!(resizes3(&f.motions[5], Measure::Length), "{f:?}");

    let f = space_freedom(&sk, sphere);
    assert_eq!(f.motions.len(), 4, "{f:?}");
    assert!(resizes3(&f.motions[3], Measure::Radius), "{f:?}");

    // A plane slides along its normal and tilts two ways; its normal
    // length is how it is stored.
    let normal = DVec3::new(0.3, -0.4, 0.5);
    let f = space_freedom(&sk, plane);
    assert_eq!(f.motions.len(), 3, "{f:?}");
    assert!(slides3(&f.motions[0], normal), "{f:?}");
    assert!(turns_across(&f.motions[1], normal), "{f:?}");
    assert!(turns_across(&f.motions[2], normal), "{f:?}");

    // A cylinder slides and tilts across its axis and resizes; sliding
    // along it and spinning about it are no motions.
    let axis = DVec3::new(0.0, 1.0, 1.0);
    let f = space_freedom(&sk, cylinder);
    assert_eq!(f.motions.len(), 5, "{f:?}");
    for m in &f.motions[..2] {
        assert!(matches!(m, Motion3::Translate { direction } if direction.dot(axis).abs() < TOL));
    }
    for m in &f.motions[2..4] {
        assert!(turns_across(m, axis), "{f:?}");
    }
    assert!(resizes3(&f.motions[4], Measure::Radius), "{f:?}");

    // A spline's control points: six rigid motions and three changes of
    // shape that name no single measure.
    let f = space_freedom(&sk, spline);
    assert_eq!(f.motions.len(), 9, "{f:?}");
    assert_eq!(
        f.motions
            .iter()
            .filter(|m| matches!(m, Motion3::Other { .. }))
            .count(),
        3,
        "{f:?}"
    );
    for e in [point, line, sphere, plane, cylinder, spline] {
        let f = space_freedom(&sk, e);
        assert_eq!(f.motions, f.own);
    }
}

#[test]
fn a_plane_with_its_offset_held_still_slides_by_its_normal_length() {
    // With d held, scaling the normal moves the plane: a slide along the
    // normal, not a gauge freedom.
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(0.0, 0.0, 2.0, 1.0);
    let start = sk.entities[plane].param_start;
    sk.fix_param(start + 3);
    let f = space_freedom(&sk, plane);
    assert_eq!(f.motions.len(), 3, "{f:?}");
    assert!(slides3(&f.motions[0], DVec3::Z), "{f:?}");
}

#[test]
fn a_point_on_a_fixed_plane_slides_in_it() {
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(0.0, 0.0, 1.0, 2.0);
    sk.fix_entity(plane);
    let p = sk.add_point(1.0, 1.0, 2.0);
    sk.add_constraint(SpaceConstraint::PointOnPlane {
        point: p.into(),
        plane,
    });
    let f = space_freedom(&sk, p);
    assert_eq!(f.motions.len(), 2, "{f:?}");
    assert!(slides3(&f.motions[0], DVec3::X), "{f:?}");
    assert!(slides3(&f.motions[1], DVec3::Y), "{f:?}");
}
