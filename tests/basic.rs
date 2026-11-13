//! Integration tests for the GCS solver.

use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::{Constraint, Sketch};

const TOL: f64 = 1e-7;

// ─── helpers ──────────────────────────────────────────────────────────────────

fn assert_near(a: f64, b: f64, label: &str) {
    assert!((a - b).abs() < TOL, "{label}: expected {b}, got {a}");
}

// ─── point constraints ────────────────────────────────────────────────────────

/// Fix a point at the origin.
#[test]
fn fixed_point() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 7.0);
    sk.add_constraint(Constraint::fix_point(p, 0.0, 0.0));
    let res = sk.solve();
    assert!(res.converged, "not converged: residual={}", res.residual);
    let c = sk.point_coords(PointRef::Point(p));
    assert_near(c.x, 0.0, "x");
    assert_near(c.y, 0.0, "y");
}

/// Two points must be coincident.
#[test]
fn coincident_points() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(5.0, 3.0);
    sk.add_constraint(Constraint::fix_point(p1, 1.0, 2.0));
    sk.add_constraint(Constraint::coincident(p1, p2));
    let res = sk.solve();
    assert!(res.converged);
    let c1 = sk.point_coords(PointRef::Point(p1));
    let c2 = sk.point_coords(PointRef::Point(p2));
    assert_near(c1.x, c2.x, "coincident x");
    assert_near(c1.y, c2.y, "coincident y");
}

/// Distance constraint: p1 fixed at origin, p2 at distance 5.
#[test]
fn point_distance() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(3.0, 4.0); // initial distance 5 — already satisfied
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
    // fix y of p2 so the system is fully constrained
    let p2_y_param = sk.entities[p2].param_start + 1;
    sk.fix_param(p2_y_param);
    sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c2 = sk.point_coords(PointRef::Point(p2));
    let dist = (c2.x * c2.x + c2.y * c2.y).sqrt();
    assert_near(dist, 5.0, "distance");
}

// ─── line constraints ─────────────────────────────────────────────────────────

/// Horizontal constraint: line becomes horizontal.
#[test]
fn horizontal_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 2.0, 1.5); // not horizontal
    sk.add_constraint(Horizontal(l));
    // Fix start point to remove translation DOF
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(l);
    assert_near(p[1], p[3], "y1 == y2 (horizontal)");
}

/// Vertical constraint.
#[test]
fn vertical_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 1.5, 2.0);
    sk.add_constraint(Vertical(l));
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(l);
    assert_near(p[0], p[2], "x1 == x2 (vertical)");
}

/// Line length constraint.
#[test]
fn line_length() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 3.0, 0.0); // length 3, want 5
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    // fix direction: keep y2=0 (horizontal)
    sk.add_constraint(Horizontal(l));
    sk.add_constraint(LineLength {
        line: l,
        length: 5.0,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(l);
    let len = ((p[2] - p[0]).powi(2) + (p[3] - p[1]).powi(2)).sqrt();
    assert_near(len, 5.0, "length");
}

/// Parallel constraint: two lines become parallel.
#[test]
fn parallel_lines() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 2.0, 0.0); // horizontal
    let l2 = sk.add_line(0.0, 1.0, 2.0, 1.5); // slightly tilted
    sk.fix_entity(l1); // fix l1 as reference
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 0.0, 1.0));
    sk.add_constraint(Parallel(l1, l2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(l1);
    let p2 = sk.entity_params(l2);
    let dy1 = p1[3] - p1[1];
    let dx1 = p1[2] - p1[0];
    let dy2 = p2[3] - p2[1];
    let dx2 = p2[2] - p2[0];
    // cross product ≈ 0
    assert_near(dx1 * dy2 - dy1 * dx2, 0.0, "parallel cross product");
}

/// Perpendicular constraint.
#[test]
fn perpendicular_lines() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 2.0, 0.0); // horizontal
    let l2 = sk.add_line(1.0, 0.0, 1.5, 1.5); // almost vertical
    sk.fix_entity(l1);
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 1.0, 0.0));
    sk.add_constraint(Perpendicular(l1, l2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(l1);
    let p2 = sk.entity_params(l2);
    let dot = (p1[2] - p1[0]) * (p2[2] - p2[0]) + (p1[3] - p1[1]) * (p2[3] - p2[1]);
    assert_near(dot, 0.0, "perpendicular dot product");
}

/// Equal length constraint.
#[test]
fn equal_length() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 4.0, 0.0); // length 4
    let l2 = sk.add_line(0.0, 1.0, 2.0, 1.0); // length 2 → should become 4
    sk.fix_entity(l1);
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 0.0, 1.0));
    sk.add_constraint(Horizontal(l2));
    sk.add_constraint(EqualLength(l1, l2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p2 = sk.entity_params(l2);
    let len2 = ((p2[2] - p2[0]).powi(2) + (p2[3] - p2[1]).powi(2)).sqrt();
    assert_near(len2, 4.0, "equal length");
}

// ─── circle constraints ───────────────────────────────────────────────────────

/// Radius constraint.
#[test]
fn circle_radius() {
    let mut sk = Sketch::new();
    let c = sk.add_circle(0.0, 0.0, 1.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Center(c), 0.0, 0.0));
    sk.add_constraint(Radius {
        circle: c,
        radius: 3.0,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(c);
    assert_near(p[2], 3.0, "radius");
}

/// Point on circle.
#[test]
fn point_on_circle() {
    let mut sk = Sketch::new();
    let c = sk.add_circle(0.0, 0.0, 2.0);
    sk.fix_entity(c);
    let p = sk.add_point(3.0, 0.0); // outside circle
    sk.add_constraint(Constraint::fix_point(p, 0.0, 3.0)); // fix to (0,3)
                                                           // The point at (0,3) is at distance 3 from origin, but circle radius is 2 — can't both hold.
                                                           // Instead: just fix y and let PointOnCircle pull x to satisfy radius.
                                                           // Re-do: fix only y of p at 0, constrain PointOnCircle.
    let mut sk2 = Sketch::new();
    let c2 = sk2.add_circle(0.0, 0.0, 2.0);
    sk2.fix_entity(c2);
    let p2 = sk2.add_point(1.5, 1.0);
    // fix y at 0 so x can adjust to land on circle
    let p2_y = sk2.entities[p2].param_start + 1;
    sk2.fix_param(p2_y);
    sk2.add_constraint(PointOnCircle {
        point: PointRef::Point(p2),
        circle: c2,
    });
    let res2 = sk2.solve();
    assert!(res2.converged, "residual={}", res2.residual);
    let coords = sk2.point_coords(PointRef::Point(p2));
    let r = (coords.x * coords.x + coords.y * coords.y).sqrt();
    assert_near(r, 2.0, "point on circle radius");
}

// ─── DOF analysis ─────────────────────────────────────────────────────────────

#[test]
fn dof_fully_constrained_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 1.0, 0.0);
    // 4 params; add: fix_start(2) + horizontal(1) + line_length(1) = 4 equations
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    sk.add_constraint(Horizontal(l));
    sk.add_constraint(LineLength {
        line: l,
        length: 5.0,
    });
    assert_eq!(sk.dof(), 0, "fully constrained");
}

#[test]
fn dof_under_constrained() {
    let mut sk = Sketch::new();
    let _p = sk.add_point(1.0, 2.0); // 2 params, 0 constraints → DOF=2
    assert_eq!(sk.dof(), 2);
}

// ─── wire BRep ───────────────────────────────────────────────────────────────

/// Angle constraint: two lines meet at a 90° angle.
#[test]
fn angle_between_lines() {
    use std::f64::consts::FRAC_PI_2;
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 2.0, 0.0);
    let l2 = sk.add_line(0.0, 0.0, 1.5, 1.5);
    sk.fix_entity(l1);
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 0.0, 0.0));
    sk.add_constraint(Angle {
        l1,
        l2,
        angle_rad: FRAC_PI_2,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(l1);
    let p2 = sk.entity_params(l2);
    let dot = (p1[2] - p1[0]) * (p2[2] - p2[0]) + (p1[3] - p1[1]) * (p2[3] - p2[1]);
    assert_near(dot, 0.0, "90° angle => dot product = 0");
}

/// EqualRadius: two circles end up with the same radius.
#[test]
fn equal_radius_circles() {
    let mut sk = Sketch::new();
    let c1 = sk.add_circle(0.0, 0.0, 3.0);
    let c2 = sk.add_circle(5.0, 0.0, 1.0);
    sk.fix_entity(c1);
    sk.add_constraint(Constraint::fix_point(PointRef::Center(c2), 5.0, 0.0));
    sk.add_constraint(EqualRadius(c1, c2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(c1);
    let p2 = sk.entity_params(c2);
    assert_near(p1[2], p2[2], "equal radius");
}

/// PointOnLine: a free point is pulled onto a fixed line.
#[test]
fn point_on_line_constraint() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 4.0, 0.0);
    sk.fix_entity(l);
    let p = sk.add_point(2.0, 3.0);
    let px_idx = sk.entities[p].param_start;
    sk.fix_param(px_idx);
    sk.add_constraint(PointOnLine {
        point: PointRef::Point(p),
        line: l,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let coords = sk.point_coords(PointRef::Point(p));
    assert_near(coords.y, 0.0, "point y should be 0 (on horizontal line)");
}

/// Tangent (circle-line): circle with fixed radius moves so it is tangent to
/// a horizontal line at y=0.  The center y-coordinate should converge to r.
#[test]
fn circle_tangent_to_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(-5.0, 0.0, 5.0, 0.0);
    sk.fix_entity(l);
    // Circle at (0, 2), radius 1.  Fix x and radius; let y move to satisfy tangency.
    let c = sk.add_circle(0.0, 2.0, 1.0);
    let c_params = sk.entities[c].param_start;
    sk.fix_param(c_params); // fix cx
    sk.fix_param(c_params + 2); // fix r
    sk.add_constraint(Tangent { circle: c, line: l });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(c);
    // dist(center, line y=0) = |cy| = r = 1
    assert_near(
        p[1].abs(),
        p[2],
        "|center y| should equal radius for tangency",
    );
}

/// CircleCircleTangent (external): c2 slides along X until externally tangent
/// to c1 (r1=2, r2=1 → expected dist = 3).
#[test]
fn circle_circle_external_tangent() {
    use rcad_gcs::constraint::Constraint::CircleCircleTangent;
    let mut sk = Sketch::new();
    let c1 = sk.add_circle(0.0, 0.0, 2.0);
    sk.fix_entity(c1);
    // c2 starts at (6, 0) with r=1.  Fix y and r; let x move.
    let c2 = sk.add_circle(6.0, 0.0, 1.0);
    let c2_params = sk.entities[c2].param_start;
    sk.fix_param(c2_params + 1); // fix cy
    sk.fix_param(c2_params + 2); // fix r
    sk.add_constraint(CircleCircleTangent {
        c1,
        c2,
        external: true,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(c1);
    let p2 = sk.entity_params(c2);
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let dist = (dx * dx + dy * dy).sqrt();
    assert_near(dist, p1[2] + p2[2], "external tangency: dist = r1 + r2");
}

/// ArcArcTangent (external): two arcs become externally tangent.
#[test]
fn arc_arc_external_tangent() {
    use rcad_gcs::constraint::Constraint::ArcArcTangent;
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 2.0, 0.0, std::f64::consts::PI);
    sk.fix_entity(a1);
    // a2 starts at (8, 0) with r=1. Fix cy and r; let cx move.
    let a2 = sk.add_arc(8.0, 0.0, 1.0, 0.0, std::f64::consts::PI);
    let a2_params = sk.entities[a2].param_start;
    sk.fix_param(a2_params + 1); // fix cy
    sk.fix_param(a2_params + 2); // fix r
    sk.add_constraint(ArcArcTangent {
        a1,
        a2,
        external: true,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(a1);
    let p2 = sk.entity_params(a2);
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let dist = (dx * dx + dy * dy).sqrt();
    assert_near(dist, p1[2] + p2[2], "arc external tangency: dist = r1 + r2");
}

/// Symmetric: two points are mirror images across a horizontal line.
#[test]
fn symmetric_across_line() {
    use rcad_gcs::constraint::Constraint::Symmetric;
    let mut sk = Sketch::new();
    // Mirror line: y=2 (horizontal)
    let mirror = sk.add_line(-5.0, 2.0, 5.0, 2.0);
    sk.fix_entity(mirror);
    // p1 fixed at (3, 0); p2 free — should end up at (3, 4)
    let p1 = sk.add_point(3.0, 0.0);
    sk.add_constraint(Constraint::fix_point(p1, 3.0, 0.0));
    let p2 = sk.add_point(3.0, 1.0); // wrong y, should move to 4
                                     // Fix x of p2 so only y is free
    let p2_x = sk.entities[p2].param_start;
    sk.fix_param(p2_x);
    sk.add_constraint(Symmetric {
        p1: PointRef::Point(p1),
        p2: PointRef::Point(p2),
        line: mirror,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c2 = sk.point_coords(PointRef::Point(p2));
    assert_near(c2.x, 3.0, "symmetric x");
    assert_near(c2.y, 4.0, "symmetric y (mirror of y=0 across y=2)");
}

// ─── new constraint tests ─────────────────────────────────────────────────────

/// Concentric: two circles share the same center.
#[test]
fn concentric_circles() {
    let mut sk = Sketch::new();
    let c1 = sk.add_circle(0.0, 0.0, 3.0);
    sk.fix_entity(c1);
    let c2 = sk.add_circle(5.0, 4.0, 1.0); // wrong center
    sk.add_constraint(Constraint::Concentric(c1, c2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(c1);
    let p2 = sk.entity_params(c2);
    assert_near(p1[0], p2[0], "concentric cx");
    assert_near(p1[1], p2[1], "concentric cy");
}

/// Midpoint: a point is pulled to the midpoint of a line.
#[test]
fn midpoint_on_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 4.0, 2.0);
    sk.fix_entity(l);
    let p = sk.add_point(0.0, 0.0); // wrong position
    sk.add_constraint(Constraint::Midpoint {
        point: PointRef::Point(p),
        line: l,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let coords = sk.point_coords(PointRef::Point(p));
    assert_near(coords.x, 2.0, "midpoint x");
    assert_near(coords.y, 1.0, "midpoint y");
}

/// Diameter: circle diameter constrained to a value.
#[test]
fn diameter_constraint() {
    let mut sk = Sketch::new();
    let c = sk.add_circle(0.0, 0.0, 1.0); // radius 1, want diameter 6 → radius 3
    sk.add_constraint(Constraint::fix_point(PointRef::Center(c), 0.0, 0.0));
    sk.add_constraint(Constraint::Diameter {
        circle: c,
        diameter: 6.0,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(c);
    assert_near(p[2], 3.0, "radius should be diameter/2");
}

/// Rectangle with diagonal: 4 lines forming a rectangle + diagonal line.
#[test]
fn rectangle_with_diagonal() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 3.0, 0.0); // bottom
    let l2 = sk.add_line(3.0, 0.0, 3.0, 2.0); // right
    let l3 = sk.add_line(3.0, 2.0, 0.0, 2.0); // top
    let l4 = sk.add_line(0.0, 2.0, 0.0, 0.0); // left
    let d = sk.add_line(0.0, 0.0, 3.0, 2.0); // diagonal
                                             // Fix the rectangle
    sk.fix_entity(l1);
    // Coincident constraints for corners
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l1),
        PointRef::LineStart(l2),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l2),
        PointRef::LineStart(l3),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l3),
        PointRef::LineStart(l4),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l4),
        PointRef::LineStart(l1),
    ));
    // Diagonal connects l1 start to l3 start
    sk.add_constraint(Constraint::coincident(
        PointRef::LineStart(d),
        PointRef::LineStart(l1),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(d),
        PointRef::LineStart(l3),
    ));
    // Opposite sides equal
    sk.add_constraint(EqualLength(l1, l3));
    sk.add_constraint(EqualLength(l2, l4));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
}

/// Four-bar linkage: 4 connected lines, fix one bar, verify DOF.
#[test]
fn four_bar_linkage() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 4.0, 0.0); // ground (fixed)
    let l2 = sk.add_line(0.0, 0.0, 1.5, 2.5); // left arm
    let l3 = sk.add_line(1.5, 2.5, 5.5, 2.5); // coupler
    let l4 = sk.add_line(5.5, 2.5, 4.0, 0.0); // right arm
    sk.fix_entity(l1);
    // Pin joints
    sk.add_constraint(Constraint::coincident(
        PointRef::LineStart(l2),
        PointRef::LineStart(l1),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l2),
        PointRef::LineStart(l3),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l3),
        PointRef::LineStart(l4),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l4),
        PointRef::LineEnd(l1),
    ));
    // Fix lengths
    sk.add_constraint(LineLength {
        line: l2,
        length: 3.0,
    });
    sk.add_constraint(LineLength {
        line: l3,
        length: 4.0,
    });
    sk.add_constraint(LineLength {
        line: l4,
        length: 3.0,
    });
    // Four-bar has 1 DOF (input angle free)
    assert_eq!(sk.dof(), 1, "four-bar should have 1 DOF");
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
}

// ============================================================================
// OCCT GCS Alignment Tests - Edge Cases
// ============================================================================

/// Over-constrained but consistent: extra constraints that are satisfied.
/// OCCT GCS coverage: redundant_constraints, over_defined_system.
#[test]
fn over_constrained_consistent() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 5.0, 0.0);
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    sk.add_constraint(Horizontal(l));
    sk.add_constraint(LineLength {
        line: l,
        length: 5.0,
    }); // length already 5
    let res = sk.solve();
    // Should still converge despite redundant constraint
    assert!(
        res.converged || res.residual < 1e-6,
        "redundant but consistent"
    );
}

/// Under-constrained: solve should still work with partial constraints.
/// OCCT GCS coverage: under_defined_system, partial_solution.
#[test]
fn under_constrained_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 3.0, 4.0);
    // Only fix start point - line can still rotate and scale
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    assert!(sk.dof() > 0, "should be under-constrained");
    let res = sk.solve();
    assert!(res.converged, "under-constrained should still solve");
}

/// Distance between two fixed points: should remain unchanged.
/// OCCT GCS coverage: fixed_entities, distance_preservation.
#[test]
fn distance_between_fixed_points() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(3.0, 4.0);
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(p2, 3.0, 4.0));
    sk.add_constraint(Constraint::point_distance(p1, p2, 5.0)); // Already satisfied
    let res = sk.solve();
    assert!(res.converged, "distance already satisfied");
    let c1 = sk.point_coords(PointRef::Point(p1));
    let c2 = sk.point_coords(PointRef::Point(p2));
    assert_near(c1.x, 0.0, "p1 x fixed");
    assert_near(c1.y, 0.0, "p1 y fixed");
    assert_near(c2.x, 3.0, "p2 x fixed");
    assert_near(c2.y, 4.0, "p2 y fixed");
}

/// Arc with fixed radius and center.
/// OCCT GCS coverage: arc_entity, circular_constraint.
#[test]
fn arc_fixed_radius() {
    let mut sk = Sketch::new();
    let a = sk.add_arc(0.0, 0.0, 2.0, 0.0, std::f64::consts::PI);
    sk.add_constraint(Constraint::fix_point(PointRef::Center(a), 0.0, 0.0));
    sk.add_constraint(Radius {
        circle: a,
        radius: 5.0,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(a);
    assert_near(p[2], 5.0, "arc radius");
}

/// Three points collinear.
/// OCCT GCS coverage: collinear_constraint, three_point_alignment.
#[test]
fn three_points_collinear() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(2.0, 2.0);
    let p3 = sk.add_point(4.0, 1.0); // Not collinear initially
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(p2, 2.0, 2.0));
    // Use point_on_line via a line through p1-p2
    let l = sk.add_line(0.0, 0.0, 2.0, 2.0);
    sk.add_constraint(Constraint::coincident(PointRef::LineStart(l), p1));
    sk.add_constraint(Constraint::coincident(PointRef::LineEnd(l), p2));
    sk.add_constraint(PointOnLine {
        point: PointRef::Point(p3),
        line: l,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    // p3 should be on the line y=x
    let c3 = sk.point_coords(PointRef::Point(p3));
    assert_near(c3.x, c3.y, "p3 on line y=x");
}

/// Circle inside circle (internal tangency).
/// OCCT GCS coverage: internal_tangency, nested_circles.
#[test]
fn circle_circle_internal_tangent() {
    use rcad_gcs::constraint::Constraint::CircleCircleTangent;
    let mut sk = Sketch::new();
    let c1 = sk.add_circle(0.0, 0.0, 5.0);
    sk.fix_entity(c1);
    let c2 = sk.add_circle(2.0, 0.0, 2.0);
    let c2_params = sk.entities[c2].param_start;
    sk.fix_param(c2_params + 1); // fix cy
    sk.fix_param(c2_params + 2); // fix r
    sk.add_constraint(CircleCircleTangent {
        c1,
        c2,
        external: false,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(c1);
    let p2 = sk.entity_params(c2);
    let dist = ((p2[0] - p1[0]).powi(2) + (p2[1] - p1[1]).powi(2)).sqrt();
    assert_near(dist, p1[2] - p2[2], "internal tangency: dist = r1 - r2");
}

/// Point distance constraint between line endpoint and point.
/// OCCT GCS coverage: point_distance, endpoint_constraints.
#[test]
fn point_line_endpoint_distance() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 4.0, 0.0); // horizontal line
    sk.fix_entity(l);
    let p = sk.add_point(2.0, 5.0);
    let px = sk.entities[p].param_start;
    sk.fix_param(px); // fix x, let y adjust
                      // Distance from point to line endpoint
    sk.add_constraint(Constraint::point_distance(
        PointRef::Point(p),
        PointRef::LineEnd(l),
        3.0,
    ));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
}

/// Multiple coincident points (chain).
/// OCCT GCS coverage: chain_constraints, multi_coincident.
#[test]
fn chain_coincident_points() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let p3 = sk.add_point(2.0, 2.0);
    let p4 = sk.add_point(3.0, 3.0);
    sk.add_constraint(Constraint::fix_point(p1, 5.0, 5.0));
    sk.add_constraint(Constraint::coincident(p1, p2));
    sk.add_constraint(Constraint::coincident(p2, p3));
    sk.add_constraint(Constraint::coincident(p3, p4));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c4 = sk.point_coords(PointRef::Point(p4));
    assert_near(c4.x, 5.0, "chain x");
    assert_near(c4.y, 5.0, "chain y");
}

/// Line tangent to circle at endpoint.
/// OCCT GCS coverage: endpoint_tangency, geometric_continuity.
#[test]
fn line_endpoint_tangent_to_circle() {
    let mut sk = Sketch::new();
    let c = sk.add_circle(0.0, 0.0, 2.0);
    sk.fix_entity(c);
    // Line starting on circle, tangent to it
    let l = sk.add_line(2.0, 0.0, 4.0, 1.0);
    sk.add_constraint(PointOnCircle {
        point: PointRef::LineStart(l),
        circle: c,
    });
    sk.add_constraint(Tangent { circle: c, line: l });
    let res = sk.solve();
    // Tangent constraint should pull the line to be tangent
    assert!(
        res.converged || res.residual < 0.1,
        "tangent line to circle"
    );
}

/// Angle constraint with 45 degrees.
/// OCCT GCS coverage: arbitrary_angle, diagonal_lines.
#[test]
fn angle_45_degrees() {
    use std::f64::consts::FRAC_PI_4;
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 4.0, 0.0); // horizontal
    let l2 = sk.add_line(0.0, 0.0, 3.0, 1.0);
    sk.fix_entity(l1);
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 0.0, 0.0));
    sk.add_constraint(Angle {
        l1,
        l2,
        angle_rad: FRAC_PI_4,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p2 = sk.entity_params(l2);
    let dx = p2[2] - p2[0];
    let dy = p2[3] - p2[1];
    let angle = (dy / dx).atan();
    assert_near(angle, FRAC_PI_4, "45° angle");
}

/// Symmetric points across vertical line.
/// OCCT GCS coverage: symmetry_vertical, mirror_constraint.
#[test]
fn symmetric_across_vertical_line() {
    use rcad_gcs::constraint::Constraint::Symmetric;
    let mut sk = Sketch::new();
    let mirror = sk.add_line(2.0, -5.0, 2.0, 5.0); // vertical line at x=2
    sk.fix_entity(mirror);
    let p1 = sk.add_point(0.0, 3.0);
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 3.0));
    let p2 = sk.add_point(1.0, 3.0);
    let p2_y = sk.entities[p2].param_start + 1;
    sk.fix_param(p2_y);
    sk.add_constraint(Symmetric {
        p1: PointRef::Point(p1),
        p2: PointRef::Point(p2),
        line: mirror,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c2 = sk.point_coords(PointRef::Point(p2));
    assert_near(c2.x, 4.0, "symmetric x (mirror of 0 across x=2)");
}

/// Concentric arcs.
/// OCCT GCS coverage: arc_concentricity, shared_center.
#[test]
fn concentric_arcs() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 3.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(5.0, 2.0, 1.0, 0.0, std::f64::consts::PI);
    sk.fix_entity(a1);
    sk.add_constraint(Constraint::Concentric(a1, a2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(a1);
    let p2 = sk.entity_params(a2);
    assert_near(p1[0], p2[0], "concentric arc cx");
    assert_near(p1[1], p2[1], "concentric arc cy");
}

/// Complex triangle with angle and length constraints.
/// OCCT GCS coverage: triangle_solver, multi_constraint_system.
#[test]
fn triangle_with_constraints() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 4.0, 0.0); // base
    let l2 = sk.add_line(4.0, 0.0, 2.0, 3.0); // right side
    let l3 = sk.add_line(2.0, 3.0, 0.0, 0.0); // left side
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l1), 0.0, 0.0));
    sk.add_constraint(Horizontal(l1));
    sk.add_constraint(LineLength {
        line: l1,
        length: 5.0,
    });
    sk.add_constraint(LineLength {
        line: l2,
        length: 4.0,
    });
    sk.add_constraint(LineLength {
        line: l3,
        length: 3.0,
    });
    // Connect vertices
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l1),
        PointRef::LineStart(l2),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l2),
        PointRef::LineStart(l3),
    ));
    sk.add_constraint(Constraint::coincident(
        PointRef::LineEnd(l3),
        PointRef::LineStart(l1),
    ));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    // Verify 3-4-5 triangle properties
    let p1 = sk.entity_params(l1);
    assert_near((p1[2] - p1[0]).abs(), 5.0, "base length");
}

/// Zero-length line handling.
/// OCCT GCS coverage: degenerate_geometry, zero_measure.
#[test]
fn zero_length_line_handling() {
    let mut sk = Sketch::new();
    let l = sk.add_line(2.0, 3.0, 2.0, 3.0); // Zero length
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 2.0, 3.0));
    let res = sk.solve();
    // Should handle degenerate case gracefully
    assert!(res.converged || res.residual < 1e-6, "zero-length line");
}

/// Point on arc (using point_on_circle for arc's circle).
/// OCCT GCS coverage: point_on_arc, circular_arc_constraint.
#[test]
fn point_on_arc() {
    let mut sk = Sketch::new();
    let a = sk.add_arc(0.0, 0.0, 3.0, 0.0, std::f64::consts::PI);
    sk.fix_entity(a);
    let p = sk.add_point(3.0, 0.0); // On the arc (radius 3, angle 0)
    let px = sk.entities[p].param_start;
    sk.fix_param(px);
    sk.add_constraint(PointOnCircle {
        point: PointRef::Point(p),
        circle: a,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c = sk.point_coords(PointRef::Point(p));
    let r = (c.x * c.x + c.y * c.y).sqrt();
    assert_near(r, 3.0, "point on arc radius");
}

// ─── DistanceX / DistanceY (directed axis-aligned distances) ─────────────────

/// Directed X distance: `x(p1) − x(p2) = d` with p1 fixed at the origin
/// drives p2 to x = −d... sign check: DistanceX(p2, p1, d) means
/// `x(p2) − x(p1) = d`, so p2 lands at x = +d.
#[test]
fn distance_x_directed() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(2.0, 3.0);
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
    let p2_y = sk.entities[p2].param_start + 1;
    sk.fix_param(p2_y); // pin y so the single equation is well-posed
    sk.add_constraint(Constraint::distance_x(p2, p1, 5.0));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c2 = sk.point_coords(PointRef::Point(p2));
    assert_near(c2.x, 5.0, "directed DistanceX(p2,p1,5) puts p2 at x=5");
    assert_near(c2.y, 3.0, "p2 y untouched");
}

/// Directed Y distance: `y(p1) − y(p2) = d` drives p2 below p1 for d > 0.
#[test]
fn distance_y_directed() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(1.0, 6.0);
    let p2 = sk.add_point(4.0, 1.0);
    sk.add_constraint(Constraint::fix_point(p1, 1.0, 6.0));
    let p2_x = sk.entities[p2].param_start;
    sk.fix_param(p2_x); // pin x so the single equation is well-posed
    sk.add_constraint(Constraint::distance_y(p1, p2, 4.0));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c2 = sk.point_coords(PointRef::Point(p2));
    assert_near(c2.y, 2.0, "DistanceY(p1,p2,4) puts p2 at y = 6−4");
    assert_near(c2.x, 4.0, "p2 x untouched");
}

/// The axis-aligned edge-margin pattern the mounting-plate hole sketch
/// uses: an anchored horizontal skeleton line whose length is driven by
/// DistanceX between its own endpoints, and a hole circle pinned by
/// DistanceX/DistanceY margins to the skeleton endpoints.
#[test]
fn distance_xy_edge_margin_pattern() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 7.5, 0.4); // drifted skeleton
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    sk.add_constraint(Horizontal(l));
    sk.add_constraint(Constraint::distance_x(
        PointRef::LineEnd(l),
        PointRef::LineStart(l),
        10.0,
    ));

    let c = sk.add_circle(6.0, 5.0, 1.0); // drifted hole
    sk.add_constraint(Radius {
        circle: c,
        radius: 1.5,
    });
    // Margin 2 from the right skeleton end (directed: end.x − center.x = 2),
    // margin 3 up from the anchored start (center.y − start.y = 3).
    sk.add_constraint(Constraint::distance_x(
        PointRef::LineEnd(l),
        PointRef::Center(c),
        2.0,
    ));
    sk.add_constraint(Constraint::distance_y(
        PointRef::Center(c),
        PointRef::LineStart(l),
        3.0,
    ));

    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let le = sk.point_coords(PointRef::LineEnd(l));
    assert_near(le.x, 10.0, "skeleton length driven to 10");
    assert_near(le.y, 0.0, "skeleton horizontal");
    let cc = sk.point_coords(PointRef::Center(c));
    assert_near(cc.x, 8.0, "hole center 2 left of the right end");
    assert_near(cc.y, 3.0, "hole center 3 above the anchor");
    assert_near(sk.entity_params(c)[2], 1.5, "radius driven");
}

/// Infeasible: both points fixed 1 apart while DistanceX demands 5 —
/// the system must report non-convergence (and never panic).
#[test]
fn distance_x_infeasible_conflict() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(p2, 1.0, 1.0));
    sk.add_constraint(Constraint::distance_x(p2, p1, 5.0));
    let res = sk.solve();
    assert!(!res.converged, "conflicting DistanceX must not converge");
}

/// Point-symmetric pair: p2 = 2·about − p1 (the stays-centered guarantee).
#[test]
fn symmetric_point_pair() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(3.0, 7.0);
    let p2 = sk.add_point(5.0, 1.0);
    let about = sk.add_point(1.0, 2.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Point(p1), 3.0, 7.0));
    sk.add_constraint(Constraint::fix_point(PointRef::Point(about), 1.0, 2.0));
    // p2 free: symmetric about (1,2) puts it at (−1, −3).
    sk.add_constraint(Constraint::SymmetricPoint {
        p1: PointRef::Point(p1),
        p2: PointRef::Point(p2),
        about: PointRef::Point(about),
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c = sk.point_coords(PointRef::Point(p2));
    assert_near(c.x, -1.0, "symmetric x = 2·1 − 3");
    assert_near(c.y, -3.0, "symmetric y = 2·2 − 7");
}
