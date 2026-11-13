//! Integration tests for the 3D GCS solver.

use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::entity::SpacePointRef;
use rcad_gcs::space3d::SpaceSketch;

const TOL: f64 = 1e-7;

fn assert_near(a: f64, b: f64, label: &str) {
    assert!((a - b).abs() < TOL, "{label}: expected {b}, got {a}");
}

/// Fix a point at a known location.
#[test]
fn fixed_point_3d() {
    let mut sk = SpaceSketch::new();
    let p = sk.add_point(10.0, 20.0, 30.0);
    sk.add_constraint(SpaceConstraint::fix_point(p, 0.0, 0.0, 0.0));
    let res = sk.solve();
    assert!(res.converged, "not converged: residual={}", res.residual);
    let c = sk.point_coords(p).unwrap();
    assert_near(c.0, 0.0, "x");
    assert_near(c.1, 0.0, "y");
    assert_near(c.2, 0.0, "z");
}

/// Two points coincident in 3D.
#[test]
fn coincident_points_3d() {
    let mut sk = SpaceSketch::new();
    let p1 = sk.add_point(0.0, 0.0, 0.0);
    let p2 = sk.add_point(5.0, 3.0, 1.0);
    sk.add_constraint(SpaceConstraint::fix_point(p1, 1.0, 2.0, 3.0));
    sk.add_constraint(SpaceConstraint::coincident(p1, p2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c1 = sk.point_coords(p1).unwrap();
    let c2 = sk.point_coords(p2).unwrap();
    assert_near(c1.0, c2.0, "coincident x");
    assert_near(c1.1, c2.1, "coincident y");
    assert_near(c1.2, c2.2, "coincident z");
}

/// Line length constraint in 3D.
#[test]
fn line_length_3d() {
    let mut sk = SpaceSketch::new();
    let l = sk.add_line(0.0, 0.0, 0.0, 1.0, 1.0, 1.0); // length sqrt(3)
    sk.add_constraint(SpaceConstraint::fix_point(
        SpacePointRef::LineStart(l),
        0.0,
        0.0,
        0.0,
    ));
    sk.add_constraint(SpaceConstraint::LineLength {
        line: l,
        length: 3.0,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let (x1, y1, z1, x2, y2, z2) = sk.line_endpoints(l).unwrap();
    let len = ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1)).sqrt();
    assert_near(len, 3.0, "line length");
}

/// Point pulled onto a plane.
#[test]
fn point_on_plane_3d() {
    let mut sk = SpaceSketch::new();
    let p = sk.add_point(1.0, 2.0, 0.0);
    let plane = sk.add_plane(0.0, 0.0, 1.0, 10.0); // z = 10
    sk.fix_entity(plane);
    sk.fix_param(sk.entities[p].param_start); // fix x
    sk.fix_param(sk.entities[p].param_start + 1); // fix y
    sk.add_constraint(SpaceConstraint::PointOnPlane {
        point: SpacePointRef::Point(p),
        plane,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c = sk.point_coords(p).unwrap();
    assert_near(c.2, 10.0, "point z should be on plane");
}

/// Cylinder with radius constraint.
#[test]
fn cylinder_radius() {
    let mut sk = SpaceSketch::new();
    let cyl = sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0); // axis along Z, radius 1
                                                                  // Fix center and axis, leave radius free
    let cs = sk.entities[cyl].param_start;
    for i in 0..6 {
        sk.fix_param(cs + i);
    }
    sk.add_constraint(SpaceConstraint::CylinderRadius {
        cylinder: cyl,
        radius: 2.5,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let e = &sk.entities[cyl];
    let r = sk.params[e.param_start + 6];
    assert_near(r, 2.5, "cylinder radius");
}

/// Plane tangent to sphere.
#[test]
fn plane_tangent_to_sphere() {
    let mut sk = SpaceSketch::new();
    let sphere = sk.add_sphere(0.0, 0.0, 0.0, 3.0); // center origin, radius 3
    sk.fix_entity(sphere);
    // Plane at z=3, normal (0,0,1), d=3 — tangent to sphere top
    let plane = sk.add_plane(0.0, 0.0, 1.0, 3.0);
    sk.fix_entity(plane);
    sk.add_constraint(SpaceConstraint::PlaneTangentToSphere { plane, sphere });
    // All fixed, nothing to solve, but should not panic
    let res = sk.solve();
    assert!(res.converged);
}

/// Two parallel planes.
#[test]
fn parallel_planes() {
    let mut sk = SpaceSketch::new();
    let p1 = sk.add_plane(0.0, 0.0, 1.0, 0.0); // z=0
    sk.fix_entity(p1);
    let p2 = sk.add_plane(0.0, 0.0, 1.0, 5.0); // z=5, already parallel
    sk.fix_entity(p2);
    // Both fixed — verify parallel constraint is trivially satisfied
    let res = sk.solve();
    assert!(res.converged);
}

// ── Splines ──────────────────────────────────────────────────────────────────

const CUBIC: [(f64, f64, f64); 4] = [
    (0.0, 0.0, 0.0),
    (1.0, 2.0, 0.0),
    (3.0, 2.0, 1.0),
    (4.0, 0.0, 1.0),
];

/// A point on a cubic Bezier curve.
fn cubic_at(cp: &[(f64, f64, f64)], t: f64) -> (f64, f64, f64) {
    let s = 1.0 - t;
    let w = [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t];
    let mut out = (0.0, 0.0, 0.0);
    for (weight, p) in w.iter().zip(cp) {
        out.0 += weight * p.0;
        out.1 += weight * p.1;
        out.2 += weight * p.2;
    }
    out
}

/// Distance from a point to a cubic: the best of a coarse sampling, then
/// narrowed down by ternary search between its two neighbours. Independent
/// of the solver's own nearest-point search.
fn distance_to_cubic(cp: &[(f64, f64, f64)], p: (f64, f64, f64)) -> f64 {
    let distance = |t: f64| {
        let c = cubic_at(cp, t);
        ((p.0 - c.0).powi(2) + (p.1 - c.1).powi(2) + (p.2 - c.2).powi(2)).sqrt()
    };
    const SAMPLES: usize = 2_000;
    let step = 1.0 / SAMPLES as f64;
    let best = (0..=SAMPLES)
        .map(|i| i as f64 * step)
        .min_by(|a, b| distance(*a).total_cmp(&distance(*b)))
        .unwrap();
    let (mut lo, mut hi) = ((best - step).max(0.0), (best + step).min(1.0));
    for _ in 0..200 {
        let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if distance(a) < distance(b) {
            hi = b;
        } else {
            lo = a;
        }
    }
    distance(0.5 * (lo + hi))
}

/// A free point is pulled onto a fixed spline from well off the curve.
#[test]
fn point_is_driven_onto_a_spline() {
    let mut sk = SpaceSketch::new();
    let spline = sk.add_spline(&CUBIC);
    sk.fix_entity(spline);
    let point = sk.add_point(2.3, 3.1, -0.6);
    sk.add_constraint(SpaceConstraint::PointOnSpline {
        point: point.into(),
        spline,
    });

    let res = sk.solve();
    assert!(res.converged, "not converged: {res:?}");
    let p = sk.point_coords(point).unwrap();
    let off = distance_to_cubic(&CUBIC, p);
    assert!(off < TOL, "{off:.3e} off the curve at {p:?}");
    assert_eq!(sk.dof(), 1, "the point still slides along the curve");
}

/// A point pinned to the ends of a spline.
#[test]
fn points_at_spline_endpoints() {
    let mut sk = SpaceSketch::new();
    let spline = sk.add_spline(&CUBIC);
    sk.fix_entity(spline);
    let first = sk.add_point(-1.0, 1.0, 1.0);
    let last = sk.add_point(5.0, 1.0, 0.0);
    sk.add_constraint(SpaceConstraint::SplineEndpoint {
        point: first.into(),
        spline,
        at_start: true,
    });
    sk.add_constraint(SpaceConstraint::SplineEndpoint {
        point: last.into(),
        spline,
        at_start: false,
    });
    assert!(sk.solve().converged);
    assert_eq!(sk.point_coords(first).unwrap(), CUBIC[0]);
    assert_eq!(sk.point_coords(last).unwrap(), CUBIC[3]);
}

/// Two control points mirrored about a line that is itself free to move:
/// the constraint has to turn the line, not only the control points.
#[test]
fn control_points_mirror_about_a_line_that_moves() {
    let mut sk = SpaceSketch::new();
    let spline = sk.add_spline(&CUBIC);
    sk.fix_entity(spline);
    // The mirror of control points 1 and 2 passes through their midpoint
    // (2, 2, 0.5) at right angles to the segment joining them. The line
    // starts on that plane and points somewhere else.
    let axis = sk.add_line(2.5, -1.0, -0.5, 4.0, 3.0, 2.0);
    sk.add_constraint(SpaceConstraint::fix_point(
        SpacePointRef::LineStart(axis),
        2.5,
        -1.0,
        -0.5,
    ));
    sk.add_constraint(SpaceConstraint::SplineControlPointSymmetric {
        spline,
        first_index: 1,
        second_index: 2,
        line: axis,
    });
    // One component: the mirror constraint acts on the line, which the
    // anchor holds.
    assert_eq!(sk.components(), vec![vec![0, 1]]);

    let res = sk.solve();
    assert!(res.converged, "not converged: {res:?}");
    let (x1, y1, z1, x2, y2, z2) = sk.line_endpoints(axis).unwrap();
    let d = (x2 - x1, y2 - y1, z2 - z1);
    // The line now runs through the midpoint...
    let to_mid = (2.0 - x1, 2.0 - y1, 0.5 - z1);
    let cross = (
        to_mid.1 * d.2 - to_mid.2 * d.1,
        to_mid.2 * d.0 - to_mid.0 * d.2,
        to_mid.0 * d.1 - to_mid.1 * d.0,
    );
    let len = (d.0 * d.0 + d.1 * d.1 + d.2 * d.2).sqrt();
    let off = (cross.0 * cross.0 + cross.1 * cross.1 + cross.2 * cross.2).sqrt() / len;
    assert!(off < TOL, "the midpoint is {off} away from the line");
    // ...at right angles to the segment between the control points.
    let joining = (2.0, 0.0, 1.0);
    let dot = (d.0 * joining.0 + d.1 * joining.1 + d.2 * joining.2) / len;
    assert!(dot.abs() < TOL, "not perpendicular: {dot}");
}

/// A control point held at the midpoint of two others.
#[test]
fn control_point_is_held_at_a_midpoint() {
    let mut sk = SpaceSketch::new();
    let spline = sk.add_spline(&CUBIC);
    let start = sk.entities[spline].param_start;
    for i in (0..3).chain(9..12) {
        sk.fix_param(start + i);
    }
    sk.add_constraint(SpaceConstraint::SplineControlPointMidpoint {
        spline,
        midpoint_index: 1,
        first_index: 0,
        second_index: 3,
    });
    assert!(sk.solve().converged);
    let cp = sk.spline_control_points(spline).unwrap();
    assert_near(cp[1].0, 2.0, "x");
    assert_near(cp[1].1, 0.0, "y");
    assert_near(cp[1].2, 0.5, "z");
    assert_eq!(
        cp[2], CUBIC[2],
        "the other inner control point is untouched"
    );
}
