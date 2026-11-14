//! Minimal conflicting sets on seeded faults (P4.4). The corpus gates
//! (`conflicts_match_the_labels` in `tests/corpus.rs` and
//! `tests/assembly_corpus.rs`, the conflict tests of
//! `tests/space3d_corpus.rs`) check every conflicting case of those corpora;
//! these cover what they do not: a conflict the linear diagnosis cannot
//! see, sets that are not minimal to start with, several conflicts at
//! once, suppressed constraints.

use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::{Conflict, Constraint, Sketch};

fn fixed(p: usize, x: f64, y: f64) -> Constraint {
    Fixed {
        point: PointRef::Point(p),
        x,
        y,
    }
}

fn conflict(constraints: &[usize]) -> Conflict {
    Conflict {
        constraints: constraints.to_vec(),
    }
}

/// Taking out the suggestion of every conflict lets the sketch solve.
fn suggestions_resolve(sk: &Sketch, found: &[Conflict]) {
    let mut fixed = sk.clone();
    for c in found {
        fixed.set_constraint_enabled(c.suggested(), false);
    }
    assert!(fixed.solve().converged, "{found:?}");
}

/// Two fixed points 5 apart and a point 1 from each: the distances cannot
/// bridge the gap.
fn short_bridge() -> Sketch {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(5.0, 0.0);
    let c = sk.add_point(2.0, 1.0);
    sk.add_constraint(fixed(a, 0.0, 0.0));
    sk.add_constraint(fixed(b, 5.0, 0.0));
    sk.add_constraint(Constraint::point_distance(a, c, 1.0));
    sk.add_constraint(Constraint::point_distance(b, c, 1.0));
    sk
}

#[test]
fn a_conflict_the_linear_diagnosis_cannot_see_is_found() {
    let sk = short_bridge();
    // Where the geometry stands the two distance rows are independent: no
    // dependency, so no conflict by linear algebra.
    let d = sk.diagnose();
    assert!(
        d.conflicting.is_empty() && d.dependency_groups.is_empty(),
        "{d:?}"
    );
    let before = sk.params_snapshot();
    let found = sk.conflicts();
    assert_eq!(found, vec![conflict(&[3, 2, 1, 0])]);
    assert_eq!(found[0].suggested(), 3);
    assert_eq!(sk.params_snapshot(), before, "the sketch is not touched");
    suggestions_resolve(&sk, &found);
}

#[test]
fn a_constraint_outside_the_conflict_is_filtered_out() {
    // The same bridge with one more point hanging off the first: the same
    // component, no part of the conflict.
    let mut sk = short_bridge();
    let d = sk.add_point(0.0, 2.0);
    sk.add_constraint(Constraint::point_distance(0, d, 2.0));
    let found = sk.conflicts();
    assert_eq!(found, vec![conflict(&[3, 2, 1, 0])]);
    suggestions_resolve(&sk, &found);
}

#[test]
fn of_a_consistent_pair_and_a_third_the_oldest_one_needed_is_kept() {
    // Pinned twice at the origin, once at (1, 0): {0, 2} and {1, 2} both
    // conflict; the newest-first filter keeps the oldest one it needs.
    let mut sk = Sketch::new();
    let p = sk.add_point(0.0, 0.0);
    sk.add_constraint(fixed(p, 0.0, 0.0));
    sk.add_constraint(fixed(p, 0.0, 0.0));
    sk.add_constraint(fixed(p, 1.0, 0.0));
    let found = sk.conflicts();
    assert_eq!(found, vec![conflict(&[2, 0])]);
    suggestions_resolve(&sk, &found);
}

#[test]
fn a_conflict_next_to_a_consistent_duplicate_leaves_the_duplicate_out() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 4.0, 0.0);
    let m = sk.add_line(4.0, 0.0, 4.0, 2.0);
    sk.add_constraint(Fixed {
        point: PointRef::LineStart(l),
        x: 0.0,
        y: 0.0,
    });
    sk.add_constraint(Horizontal(l));
    sk.add_constraint(LineLength {
        line: l,
        length: 4.0,
    });
    sk.add_constraint(Coincident(PointRef::LineEnd(l), PointRef::LineStart(m)));
    sk.add_constraint(Vertical(m));
    sk.add_constraint(Vertical(m));
    sk.add_constraint(LineLength {
        line: m,
        length: 2.0,
    });
    // Three along x where the line is 4 long and horizontal.
    sk.add_constraint(DistanceX {
        p1: PointRef::LineEnd(l),
        p2: PointRef::LineStart(l),
        distance: 3.0,
    });
    let found = sk.conflicts();
    assert_eq!(found, vec![conflict(&[7, 2])]);
    suggestions_resolve(&sk, &found);
}

#[test]
fn two_conflicts_in_one_component_are_both_found() {
    // A point pinned at two places, and a line from it with two lengths.
    let mut sk = Sketch::new();
    let p = sk.add_point(0.0, 0.0);
    let l = sk.add_line(0.0, 0.0, 2.0, 0.0);
    sk.add_constraint(fixed(p, 0.0, 0.0));
    sk.add_constraint(fixed(p, 1.0, 0.0));
    sk.add_constraint(Coincident(PointRef::LineStart(l), PointRef::Point(p)));
    sk.add_constraint(LineLength {
        line: l,
        length: 2.0,
    });
    sk.add_constraint(LineLength {
        line: l,
        length: 3.0,
    });
    let mut found = sk.conflicts();
    found.sort_by_key(|c| c.constraints.clone());
    assert_eq!(found, vec![conflict(&[1, 0]), conflict(&[4, 3])]);
    suggestions_resolve(&sk, &found);
}

#[test]
fn conflicts_in_separate_components_are_found_apart() {
    let mut sk = Sketch::new();
    let p = sk.add_point(0.0, 0.0);
    let q = sk.add_point(5.0, 5.0);
    sk.add_constraint(fixed(p, 0.0, 0.0));
    sk.add_constraint(fixed(q, 5.0, 5.0));
    sk.add_constraint(fixed(p, 1.0, 0.0));
    sk.add_constraint(fixed(q, 5.0, 6.0));
    let mut found = sk.conflicts();
    found.sort_by_key(|c| c.constraints.clone());
    assert_eq!(found, vec![conflict(&[2, 0]), conflict(&[3, 1])]);
    suggestions_resolve(&sk, &found);
}

#[test]
fn a_suppressed_constraint_takes_no_part() {
    let mut sk = Sketch::new();
    let p = sk.add_point(0.0, 0.0);
    sk.add_constraint(fixed(p, 0.0, 0.0));
    let other = sk.add_constraint(fixed(p, 1.0, 0.0));
    assert_eq!(sk.conflicts(), vec![conflict(&[1, 0])]);
    sk.set_constraint_enabled(other, false);
    assert_eq!(sk.conflicts(), Vec::new());
    // Still suppressed afterwards.
    assert!(sk.solve().converged);
    assert_eq!(sk.point_coords(PointRef::Point(p)).x, 0.0);
}
