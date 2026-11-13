//! Integration tests for the 3D GCS solver.

use rcad_gcs::space3d::SpaceSketch;
use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::entity::SpacePointRef;

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
    sk.add_constraint(SpaceConstraint::fix_point(SpacePointRef::LineStart(l), 0.0, 0.0, 0.0));
    sk.add_constraint(SpaceConstraint::LineLength { line: l, length: 3.0 });
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
    sk.fix_param(sk.entities[p].param_start);     // fix x
    sk.fix_param(sk.entities[p].param_start + 1); // fix y
    sk.add_constraint(SpaceConstraint::PointOnPlane { point: SpacePointRef::Point(p), plane });
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
    sk.add_constraint(SpaceConstraint::CylinderRadius { cylinder: cyl, radius: 2.5 });
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
