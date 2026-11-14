//! The elimination of trivial equations before a 2D solve (P3.1).
//!
//! Coincidences, fixed points, axis distances, horizontal and vertical
//! lines, equal radii, concentric circles, radii and diameters each say that
//! one parameter is another plus a constant, or a value. The solve puts the
//! parameters where those equations say and iterates only on what is left;
//! an equation that closes a cycle of them stays an ordinary equation, so a
//! contradiction among them is still met (or not) and still diagnosed.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::{Sketch, SolveStatus};

/// A rectangle with a fixed corner, chained coincidences, a horizontal and
/// a vertical pair and two lengths: fully constrained.
fn rectangle(w: f64, h: f64) -> Sketch {
    let mut sk = Sketch::new();
    let origin = sk.add_point(0.0, 0.0);
    let b = sk.add_line(0.0, 0.0, w, 0.0);
    let r = sk.add_line(w, 0.0, w, h);
    let t = sk.add_line(w, h, 0.0, h);
    let l = sk.add_line(0.0, h, 0.0, 0.0);
    sk.add_constraint(Constraint::fix_point(origin, 0.0, 0.0));
    sk.add_constraint(Coincident(PointRef::LineStart(b), PointRef::Point(origin)));
    for (a, z) in [(b, r), (r, t), (t, l), (l, b)] {
        sk.add_constraint(Coincident(PointRef::LineEnd(a), PointRef::LineStart(z)));
    }
    sk.add_constraint(Horizontal(b));
    sk.add_constraint(Horizontal(t));
    sk.add_constraint(Vertical(r));
    sk.add_constraint(Vertical(l));
    sk.add_constraint(LineLength { line: b, length: w });
    sk.add_constraint(LineLength { line: r, length: h });
    sk
}

/// Points on the x axis tied by axis distances, the first fixed: every
/// equation is an alias, so nothing is left to iterate on.
fn axis_chain(offsets: &[f64]) -> (Sketch, Vec<usize>) {
    let mut sk = Sketch::new();
    let ids: Vec<usize> = (0..=offsets.len())
        .map(|i| sk.add_point(0.37 * i as f64, 0.1 * i as f64))
        .collect();
    sk.add_constraint(Constraint::fix_point(ids[0], 1.0, 2.0));
    // Tie the points in an order that joins classes from both ends, so the
    // offsets have to be carried through merged classes.
    let mut order: Vec<usize> = (0..offsets.len()).collect();
    order.sort_by_key(|&i| (i % 2, std::cmp::Reverse(i)));
    for i in order {
        sk.add_constraint(Constraint::distance_x(ids[i + 1], ids[i], offsets[i]));
        sk.add_constraint(Constraint::distance_y(ids[i + 1], ids[i], -offsets[i]));
    }
    (sk, ids)
}

#[test]
fn a_system_of_aliases_alone_is_placed_without_iterating() {
    let offsets = [0.5, -1.25, 3.0, 0.125, -2.0, 1.5, 0.75];
    let (mut sk, ids) = axis_chain(&offsets);
    let result = sk.solve();
    assert_eq!(result.status, SolveStatus::Converged);
    assert_eq!(result.iterations, 0, "nothing should be left to iterate on");
    let (mut x, mut y) = (1.0, 2.0);
    for (i, &id) in ids.iter().enumerate() {
        if i > 0 {
            x += offsets[i - 1];
            y -= offsets[i - 1];
        }
        let p = sk.point_coords(PointRef::Point(id));
        assert!(
            (p.x - x).abs() < 1e-12 && (p.y - y).abs() < 1e-12,
            "point {i} at {p:?}, expected ({x}, {y})"
        );
    }
}

#[test]
fn a_drifted_rectangle_comes_back_to_the_branch_it_was_drawn_on() {
    // Drifted by up to 3 on every free coordinate (half the width), the
    // solve used to come back mirrored in most of these starts. With the
    // coincidences, the pairs and the fixed corner eliminated, what is left
    // is the two lengths, and the side they act on is placed where the
    // drifted corners say on average.
    for seed in 0..16u64 {
        let mut sk = rectangle(6.0, 4.0);
        let authored = sk.params_snapshot();
        sk.drift_free_params(3.0, 0x5EED + seed);
        let result = sk.solve();
        assert!(result.converged, "seed {seed}: {result:?}");
        let worst = authored
            .iter()
            .zip(sk.params_snapshot())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(worst < 1e-9, "seed {seed}: came back {worst} away");
    }
}

#[test]
fn a_contradiction_among_aliases_is_not_solved_away_and_is_diagnosed() {
    // x1 − x0 = 1, x2 − x1 = 1, x2 − x0 = 3: the third closes a cycle that
    // the first two say should be 2.
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    // Started away from where the aliases put them, so that the elimination
    // moves them before it finds it cannot solve.
    let p1 = sk.add_point(1.3, 0.2);
    let p2 = sk.add_point(2.4, -0.1);
    sk.add_constraint(Constraint::fix_point(p0, 0.0, 0.0));
    sk.add_constraint(Constraint::distance_x(p1, p0, 1.0));
    sk.add_constraint(Constraint::distance_x(p2, p1, 1.0));
    let closing = sk.add_constraint(Constraint::distance_x(p2, p0, 3.0));
    let before = sk.params_snapshot();

    let result = sk.solve();
    assert!(!result.converged, "{result:?}");
    assert_eq!(
        sk.params_snapshot(),
        before,
        "a failed solve puts the start back"
    );
    let diagnosis = sk.diagnose();
    assert!(
        diagnosis.conflicting.contains(&closing),
        "the closing constraint is named: {diagnosis:?}"
    );
}

#[test]
fn a_consistent_cycle_of_aliases_solves_and_is_diagnosed_redundant() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.2, 0.0);
    let p2 = sk.add_point(1.7, 0.0);
    sk.add_constraint(Constraint::fix_point(p0, 0.0, 0.0));
    sk.add_constraint(Constraint::distance_x(p1, p0, 0.1));
    sk.add_constraint(Constraint::distance_x(p2, p1, 0.2));
    // 0.1 + 0.2 is not 0.3 in floating point: the cycle holds to rounding.
    let closing = sk.add_constraint(Constraint::distance_x(p2, p0, 0.3));
    let cycle = [closing - 2, closing - 1, closing];

    let result = sk.solve();
    assert!(result.converged, "{result:?}");
    let diagnosis = sk.diagnose();
    assert!(diagnosis.conflicting.is_empty(), "{diagnosis:?}");
    assert!(
        !diagnosis.redundant.is_empty() && diagnosis.redundant.iter().all(|c| cycle.contains(c)),
        "the cycle is reported redundant: {diagnosis:?}"
    );
}

#[test]
fn eliminated_coordinates_are_equal_exactly() {
    // A free point coincident with a line end and the circle centre, the
    // line horizontal; only a distance from a fixed point is nonlinear.
    let mut sk = Sketch::new();
    let anchor = sk.add_point(0.0, 0.0);
    let p = sk.add_point(2.9, 0.4);
    let line = sk.add_line(-1.0, 0.3, 3.2, -0.2);
    let circle = sk.add_circle(3.1, 0.1, 0.5);
    sk.add_constraint(Constraint::fix_point(anchor, 0.0, 0.0));
    sk.add_constraint(Coincident(PointRef::Point(p), PointRef::LineEnd(line)));
    sk.add_constraint(Coincident(PointRef::Center(circle), PointRef::Point(p)));
    sk.add_constraint(Horizontal(line));
    sk.add_constraint(Constraint::point_distance(anchor, p, 3.0));
    let result = sk.solve();
    assert!(result.converged, "{result:?}");
    let at = sk.point_coords(PointRef::Point(p));
    assert_eq!(at, sk.point_coords(PointRef::LineEnd(line)));
    assert_eq!(at, sk.point_coords(PointRef::Center(circle)));
    assert_eq!(
        sk.point_coords(PointRef::LineStart(line)).y,
        sk.point_coords(PointRef::LineEnd(line)).y
    );
    assert!((at.length() - 3.0).abs() < 1e-10);
}

#[test]
fn a_coordinate_that_is_not_a_number_still_diverges() {
    // A fixed point would put the NaN coordinate in place; the solve still
    // reports the input as it is rather than mend it.
    let mut sk = Sketch::new();
    let p = sk.add_point(f64::NAN, 1.0);
    let q = sk.add_point(2.0, 1.0);
    sk.add_constraint(Constraint::fix_point(p, 0.0, 0.0));
    sk.add_constraint(Constraint::point_distance(p, q, 1.0));
    let result = sk.solve();
    assert_eq!(result.status, SolveStatus::Diverged, "{result:?}");
    assert!(sk.point_coords(PointRef::Point(p)).x.is_nan());
}
