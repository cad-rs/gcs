//! What a 2D solve keeps of the structure between solves (P3.4).
//!
//! A sketch keeps its components and, per component, what its solve needs
//! that depends on the structure alone, while the topology revision and the
//! entity list stay the same. Every test here compares a sketch that was
//! solved before the edit (and so may reuse what it kept) with one built
//! directly in the edited state and never solved: the two must agree bit
//! for bit.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::Sketch;

fn bits(sk: &Sketch) -> Vec<u64> {
    sk.params_snapshot().iter().map(|v| v.to_bits()).collect()
}

/// Two rectangles, one fixed at the origin and one hanging on its corner;
/// `width` is the bottom length of the first.
fn rectangles(width: f64) -> Sketch {
    let mut sk = Sketch::new();
    let mut anchor: Option<PointRef> = None;
    for i in 0..2 {
        let x = i as f64 * 7.0;
        let origin = sk.add_point(x, 0.0);
        let b = sk.add_line(x, 0.0, x + 6.0, 0.1);
        let r = sk.add_line(x + 6.0, 0.0, x + 6.2, 4.0);
        let t = sk.add_line(x + 6.0, 4.0, x, 4.1);
        let l = sk.add_line(x, 4.0, x, 0.0);
        sk.add_constraint(match anchor {
            None => Constraint::fix_point(origin, 0.0, 0.0),
            Some(at) => Coincident(PointRef::Point(origin), at),
        });
        sk.add_constraint(Coincident(PointRef::LineStart(b), PointRef::Point(origin)));
        for (a, z) in [(b, r), (r, t), (t, l), (l, b)] {
            sk.add_constraint(Coincident(PointRef::LineEnd(a), PointRef::LineStart(z)));
        }
        sk.add_constraint(Horizontal(b));
        sk.add_constraint(Horizontal(t));
        sk.add_constraint(Vertical(r));
        sk.add_constraint(Vertical(l));
        sk.add_constraint(LineLength {
            line: b,
            length: if i == 0 { width } else { 5.0 },
        });
        sk.add_constraint(LineLength {
            line: r,
            length: 4.0,
        });
        anchor = Some(PointRef::LineEnd(b));
    }
    sk
}

#[test]
fn a_dimension_edit_resolves_as_a_fresh_sketch_would() {
    let mut kept = rectangles(6.0);
    assert!(kept.solve().converged);
    let start = kept.params_snapshot();
    // The first rectangle's width is constraint 10.
    assert!(kept.set_constraint_value(10, 8.5));
    assert!(kept.solve().converged);

    let mut fresh = rectangles(8.5);
    fresh.set_params(start);
    assert!(fresh.solve().converged);
    assert_eq!(bits(&kept), bits(&fresh));
}

#[test]
fn structural_edits_after_a_solve_are_seen() {
    // Each edit, made after a solve, against the same edit made before any.
    type Edit = fn(&mut Sketch);
    let edits: [(&str, Edit); 4] = [
        ("suppress a length", |sk| {
            sk.set_constraint_enabled(11, false);
        }),
        ("add a distance", |sk| {
            sk.add_constraint(Constraint::point_distance(
                PointRef::Point(0),
                PointRef::LineEnd(6),
                9.0,
            ));
        }),
        ("fix a coordinate", |sk| {
            let x = sk.entity_param_range(7).start;
            sk.fix_param(x);
        }),
        ("remove a constraint", |sk| {
            sk.remove_constraint(9);
        }),
    ];
    for (name, edit) in edits {
        let mut kept = rectangles(6.0);
        assert!(kept.solve().converged, "{name}");
        let start = kept.params_snapshot();
        edit(&mut kept);
        let a = kept.solve();

        let mut fresh = rectangles(6.0);
        fresh.set_params(start);
        edit(&mut fresh);
        let b = fresh.solve();
        assert_eq!(a.status, b.status, "{name}");
        assert_eq!(a.iterations, b.iterations, "{name}");
        assert_eq!(bits(&kept), bits(&fresh), "{name}");
    }
}

#[test]
fn entities_edited_directly_are_seen() {
    // The entity list is public: swapping two points re-points the
    // constraints at the other's parameters without any revision change.
    let build = || {
        let mut sk = Sketch::new();
        let p = sk.add_point(0.3, 0.2);
        let q = sk.add_point(4.0, 5.5);
        let r = sk.add_point(1.0, 1.0);
        sk.add_constraint(Constraint::fix_point(p, 0.0, 0.0));
        sk.add_constraint(Constraint::point_distance(p, r, 2.0));
        sk.add_constraint(Constraint::fix_point(q, 5.0, 5.0));
        sk
    };
    let mut kept = build();
    assert!(kept.solve().converged);
    let start = kept.params_snapshot();
    kept.entities.swap(0, 1);
    let a = kept.solve();

    let mut fresh = build();
    fresh.set_params(start);
    fresh.entities.swap(0, 1);
    let b = fresh.solve();
    assert!(b.converged, "{b:?}");
    assert_eq!(a.status, b.status);
    assert_eq!(bits(&kept), bits(&fresh));
}

#[test]
fn a_clone_edited_on_its_own_does_not_disturb_the_original() {
    let mut original = rectangles(6.0);
    assert!(original.solve().converged);
    let mut copy = original.clone();
    copy.replace_constraint(
        10,
        LineLength {
            line: 1,
            length: 3.0,
        },
    );
    assert!(copy.solve().converged);

    let start = original.params_snapshot();
    assert!(original.set_constraint_value(11, 4.5));
    assert!(original.solve().converged);
    let mut fresh = rectangles(6.0);
    fresh.set_params(start);
    assert!(fresh.set_constraint_value(11, 4.5));
    assert!(fresh.solve().converged);
    assert_eq!(bits(&original), bits(&fresh));
}
