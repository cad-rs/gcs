//! Solving a component block by block (P3.2).
//!
//! A component whose equations fall into several blocks — sets of
//! equations that can be solved for as many unknowns, each using only the
//! unknowns of the blocks before it — is solved one block after the other.
//! A block that does not converge sends the whole component back to where
//! it started and to the solve of the whole system, so a failure ends as
//! it would without the blocks.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::entity::PointRef;
use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::SpaceSketch;
use rcad_gcs::Sketch;

#[test]
fn a_block_that_cannot_be_met_puts_the_3d_component_back() {
    // Three held points; d hangs on all three (a block of three unknowns),
    // e on a, b and d, 50 away from a point 1 from a: no solution.
    let mut sk = SpaceSketch::new();
    let a = sk.add_point(0.0, 0.0, 0.0);
    let b = sk.add_point(1.0, 0.0, 0.0);
    let c = sk.add_point(0.0, 1.0, 0.0);
    for id in [a, b, c] {
        sk.fix_entity(id);
    }
    let d = sk.add_point(0.4, 0.3, 0.9);
    let e = sk.add_point(1.2, 1.1, 1.3);
    for (from, distance) in [(a, 1.0), (b, 1.2), (c, 1.1)] {
        sk.add_constraint(SpaceConstraint::point_distance(d, from, distance));
    }
    for (from, distance) in [(a, 1.0), (b, 1.0), (d, 50.0)] {
        sk.add_constraint(SpaceConstraint::point_distance(e, from, distance));
    }
    let start = sk.params_snapshot();
    let result = sk.solve();
    assert!(!result.converged, "{result:?}");
    assert_eq!(sk.params_snapshot(), start, "the start is put back");
}

#[test]
fn a_block_that_cannot_be_met_puts_the_2d_component_back() {
    // Fixed points (eliminated as aliases) a and b; c 1 from each (a block
    // of two unknowns); d 1 from a and 30 from c: no solution.
    let mut sk = Sketch::new();
    let a = sk.add_point(0.1, -0.1);
    let b = sk.add_point(1.1, 0.2);
    let c = sk.add_point(0.4, 0.9);
    let d = sk.add_point(-0.6, 0.7);
    sk.add_constraint(Constraint::fix_point(a, 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(b, 1.0, 0.0));
    sk.add_constraint(Constraint::point_distance(c, a, 1.0));
    sk.add_constraint(Constraint::point_distance(c, b, 1.0));
    sk.add_constraint(Constraint::point_distance(d, a, 1.0));
    sk.add_constraint(Constraint::point_distance(d, c, 30.0));
    let start = sk.params_snapshot();
    let result = sk.solve();
    assert!(!result.converged, "{result:?}");
    assert_eq!(sk.params_snapshot(), start, "the start is put back");
}

#[test]
fn each_point_of_a_chain_settles_next_to_where_it_stands() {
    // A chain of points, each 1 from the one before, along a slight curve;
    // nudged, every point comes back to its own side of the chain.
    let mut sk = Sketch::new();
    let ids: Vec<usize> = (0..12)
        .map(|i| {
            let t = i as f64;
            sk.add_point(t * 0.98, 0.2 * (t * 0.5).sin())
        })
        .collect();
    sk.add_constraint(Constraint::fix_point(ids[0], 0.0, 0.0));
    for w in ids.windows(2) {
        sk.add_constraint(Constraint::point_distance(w[0], w[1], 1.0));
        sk.add_constraint(Constraint::distance_y(
            PointRef::Point(w[1]),
            PointRef::Point(w[0]),
            0.0,
        ));
    }
    sk.drift_free_params(0.3, 11);
    let result = sk.solve();
    assert!(result.converged, "{result:?}");
    for (i, &id) in ids.iter().enumerate() {
        let p = sk.point_coords(PointRef::Point(id));
        assert!(
            (p.x - i as f64).abs() < 1e-9 && p.y.abs() < 1e-9,
            "point {i} at {p:?}"
        );
    }
}
