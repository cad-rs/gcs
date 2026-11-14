//! 3D GCS corpus (solver 3D-0 acceptance).
//!
//! Convergence, conflict detection, redundancy diagnosis, determinism,
//! incremental re-solve, and plane parameterization hardening over the
//! [`SpaceSketch`] API. The corpus mirrors the 2D solver's discipline: every
//! case asserts on the solver/diagnosis answer, not on intermediate trace
//! output.

use std::time::Instant;

use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::entity::SpacePointRef;
use rcad_gcs::space3d::SpaceSketch;
use rcad_gcs::EntityStatus;

// ── helpers ──────────────────────────────────────────────────────────────────

fn dist(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    let (dx, dy, dz) = (a.0 - b.0, a.1 - b.1, a.2 - b.2);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn solve(sk: &mut SpaceSketch) {
    let r = sk.solve();
    assert!(
        r.converged,
        "solve did not converge: residual={}, iterations={}",
        r.residual, r.iterations
    );
}

fn assert_near(a: f64, b: f64, tol: f64, label: &str) {
    assert!((a - b).abs() < tol, "{label}: expected {b}, got {a}");
}

fn plane_normal(sk: &SpaceSketch, plane: usize) -> (f64, f64, f64) {
    let e = &sk.entities[plane];
    (
        sk.params[e.param(0)],
        sk.params[e.param(1)],
        sk.params[e.param(2)],
    )
}

fn point3(sk: &SpaceSketch, id: usize) -> (f64, f64, f64) {
    sk.point_coords(id).unwrap()
}

// ── convergence corpus ───────────────────────────────────────────────────────

#[test]
fn convergence_corpus() {
    // 1. Distance chain between three free points anchored at one fixed point.
    {
        let mut sk = SpaceSketch::new();
        let p0 = sk.add_point(0.0, 0.0, 0.0);
        let p1 = sk.add_point(2.0, 1.0, 0.0);
        let p2 = sk.add_point(4.0, 3.0, 1.0);
        sk.add_constraint(SpaceConstraint::fix_point(p0, 0.0, 0.0, 0.0));
        sk.add_constraint(SpaceConstraint::point_distance(p0, p1, 3.0));
        sk.add_constraint(SpaceConstraint::point_distance(p1, p2, 4.0));
        solve(&mut sk);
        assert_near(
            dist(point3(&sk, p0), point3(&sk, p1)),
            3.0,
            1e-8,
            "chain d01",
        );
        assert_near(
            dist(point3(&sk, p1), point3(&sk, p2)),
            4.0,
            1e-8,
            "chain d12",
        );
    }

    // 2. Driven line: fixed start, prescribed length and axial direction.
    {
        let mut sk = SpaceSketch::new();
        let l = sk.add_line(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        sk.add_constraint(SpaceConstraint::fix_point(
            SpacePointRef::LineStart(l),
            0.0,
            0.0,
            0.0,
        ));
        sk.add_constraint(SpaceConstraint::LineLength {
            line: l,
            length: 5.0,
        });
        sk.add_constraint(SpaceConstraint::LineParallelToDir {
            line: l,
            dx: 0.0,
            dy: 0.0,
            dz: 1.0,
        });
        solve(&mut sk);
        let (_, _, _, x2, y2, z2) = sk.line_endpoints(l).unwrap();
        assert_near(
            (x2 * x2 + y2 * y2 + z2 * z2).sqrt(),
            5.0,
            1e-8,
            "line length",
        );
        assert_near(x2, 0.0, 1e-8, "line dx");
        assert_near(y2, 0.0, 1e-8, "line dy");
    }

    // 3. Point pulled onto a sphere surface (non-equatorial start).
    {
        let mut sk = SpaceSketch::new();
        let p = sk.add_point(1.0, 2.0, 1.0);
        let s = sk.add_sphere(0.0, 0.0, 0.0, 3.0);
        sk.fix_entity(s);
        let start = sk.entities[p].param_start;
        sk.fix_param(start);
        sk.fix_param(start + 1);
        sk.add_constraint(SpaceConstraint::PointOnSphere {
            point: SpacePointRef::Point(p),
            sphere: s,
        });
        solve(&mut sk);
        assert_near(
            dist(point3(&sk, p), (0.0, 0.0, 0.0)),
            3.0,
            1e-8,
            "point on sphere",
        );
    }

    // 4. Point pulled onto a cylinder surface.
    {
        let mut sk = SpaceSketch::new();
        let p = sk.add_point(1.0, 1.0, 0.0);
        let c = sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0);
        sk.fix_entity(c);
        let start = sk.entities[p].param_start;
        sk.fix_param(start);
        sk.fix_param(start + 2);
        sk.add_constraint(SpaceConstraint::PointOnCylinder {
            point: SpacePointRef::Point(p),
            cylinder: c,
        });
        solve(&mut sk);
        let q = point3(&sk, p);
        let perp = (q.0 * q.0 + q.1 * q.1).sqrt();
        assert_near(perp, 2.0, 1e-8, "point on cylinder");
    }

    // 5. Line made perpendicular to a fixed plane.
    {
        let mut sk = SpaceSketch::new();
        let l = sk.add_line(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let plane = sk.add_plane(0.0, 0.0, 1.0, 0.0);
        sk.fix_entity(plane);
        sk.add_constraint(SpaceConstraint::fix_point(
            SpacePointRef::LineStart(l),
            0.0,
            0.0,
            0.0,
        ));
        sk.add_constraint(SpaceConstraint::LinePerpendicularPlane { line: l, plane });
        solve(&mut sk);
        let (_, _, _, x2, y2, z2) = sk.line_endpoints(l).unwrap();
        assert_near(x2, 0.0, 1e-8, "perp-plane dx");
        assert_near(y2, 0.0, 1e-8, "perp-plane dy");
        assert!(z2.abs() > 0.5, "line should keep a z direction, got {z2}");
    }

    // 6. Angle between two lines driven to 90 degrees.
    {
        let mut sk = SpaceSketch::new();
        let l1 = sk.add_line(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        sk.fix_entity(l1);
        let l2 = sk.add_line(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        sk.add_constraint(SpaceConstraint::fix_point(
            SpacePointRef::LineStart(l2),
            0.0,
            0.0,
            0.0,
        ));
        sk.add_constraint(SpaceConstraint::AngleLineLine {
            l1,
            l2,
            angle_rad: std::f64::consts::FRAC_PI_2,
        });
        solve(&mut sk);
        let (_, _, _, x2, y2, z2) = sk.line_endpoints(l2).unwrap();
        assert!(x2.abs() < 1e-8, "perpendicular x: {x2}");
        let _ = (y2, z2);
    }

    // 7. Free plane fully determined by a normal target and a fixed point.
    {
        let mut sk = SpaceSketch::new();
        let plane = sk.add_plane(0.3, -0.4, 0.5, 1.2);
        let p = sk.add_point(0.0, 0.0, 5.0);
        sk.fix_entity(p);
        sk.add_constraint(SpaceConstraint::PlaneNormal {
            plane,
            nx: 0.0,
            ny: 0.0,
            nz: 1.0,
        });
        sk.add_constraint(SpaceConstraint::PointOnPlane {
            point: SpacePointRef::Point(p),
            plane,
        });
        solve(&mut sk);
        let n = plane_normal(&sk, plane);
        assert_near(n.0, 0.0, 1e-8, "plane nx");
        assert_near(n.1, 0.0, 1e-8, "plane ny");
        assert_near(n.2, 1.0, 1e-8, "plane nz");
        let e = &sk.entities[plane];
        assert_near(sk.params[e.param(3)], 5.0, 1e-8, "plane d");
        assert_eq!(
            sk.diagnose().dof,
            0,
            "determined plane is fully constrained"
        );
    }

    // 8. Externally tangent spheres with one free center coordinate.
    {
        let mut sk = SpaceSketch::new();
        let s1 = sk.add_sphere(0.0, 0.0, 0.0, 2.0);
        sk.fix_entity(s1);
        let s2 = sk.add_sphere(5.0, 0.0, 0.0, 1.0);
        let start = sk.entities[s2].param_start;
        for k in [1usize, 2, 3] {
            sk.fix_param(start + k);
        }
        sk.add_constraint(SpaceConstraint::SphereTangent {
            s1,
            s2,
            external: true,
        });
        solve(&mut sk);
        let e = &sk.entities[s2];
        let cx = sk.params[e.param(0)];
        assert_near(cx, 3.0, 1e-8, "sphere tangent center x");
    }

    // 9. Tangent cylinders with one free axis-distance coordinate.
    {
        let mut sk = SpaceSketch::new();
        let c1 = sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
        sk.fix_entity(c1);
        let c2 = sk.add_cylinder(5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.5);
        let start = sk.entities[c2].param_start;
        for k in [1usize, 2, 3, 4, 5, 6] {
            sk.fix_param(start + k);
        }
        sk.add_constraint(SpaceConstraint::CylinderTangent {
            c1,
            c2,
            external: true,
        });
        solve(&mut sk);
        let e = &sk.entities[c2];
        assert_near(
            sk.params[e.param(0)],
            2.5,
            1e-8,
            "cylinder tangent center x",
        );
    }

    // 10. Plane tangent to a sphere with a free offset.
    {
        let mut sk = SpaceSketch::new();
        let sphere = sk.add_sphere(0.0, 0.0, 0.0, 3.0);
        sk.fix_entity(sphere);
        let plane = sk.add_plane(0.0, 0.0, 1.0, 1.0);
        let start = sk.entities[plane].param_start;
        for k in [0usize, 1, 2] {
            sk.fix_param(start + k);
        }
        sk.add_constraint(SpaceConstraint::PlaneTangentToSphere { plane, sphere });
        solve(&mut sk);
        let e = &sk.entities[plane];
        assert_near(sk.params[e.param(3)], 3.0, 1e-8, "tangent plane offset");
    }

    // 11. Line pulled into a fixed plane.
    {
        let mut sk = SpaceSketch::new();
        let plane = sk.add_plane(0.0, 0.0, 1.0, 2.0);
        sk.fix_entity(plane);
        let l = sk.add_line(0.0, 0.0, 0.0, 1.0, 1.0, 0.0);
        let start = sk.entities[l].param_start;
        for k in [0usize, 1, 3, 4] {
            sk.fix_param(start + k);
        }
        sk.add_constraint(SpaceConstraint::LineInPlane { line: l, plane });
        solve(&mut sk);
        assert_near(sk.params[start + 2], 2.0, 1e-8, "line start z");
        assert_near(sk.params[start + 5], 2.0, 1e-8, "line end z");
    }

    // 12. Concentric spheres.
    {
        let mut sk = SpaceSketch::new();
        let s1 = sk.add_sphere(1.0, 2.0, 3.0, 1.0);
        sk.fix_entity(s1);
        let s2 = sk.add_sphere(5.0, 5.0, 5.0, 2.0);
        let start = sk.entities[s2].param_start;
        sk.fix_param(start + 3);
        sk.add_constraint(SpaceConstraint::Concentric3D(s1, s2));
        solve(&mut sk);
        let e = &sk.entities[s2];
        let c = (
            sk.params[e.param(0)],
            sk.params[e.param(1)],
            sk.params[e.param(2)],
        );
        assert_near(dist(c, (1.0, 2.0, 3.0)), 0.0, 1e-8, "concentric centers");
    }

    // 13. Parallel lines.
    {
        let mut sk = SpaceSketch::new();
        let l1 = sk.add_line(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        sk.fix_entity(l1);
        let l2 = sk.add_line(2.0, 0.0, 0.0, 3.0, 1.0, 1.0);
        sk.add_constraint(SpaceConstraint::fix_point(
            SpacePointRef::LineStart(l2),
            2.0,
            0.0,
            0.0,
        ));
        sk.add_constraint(SpaceConstraint::LineParallelLine(l1, l2));
        solve(&mut sk);
        let (_, _, _, x2, y2, z2) = sk.line_endpoints(l2).unwrap();
        assert_near(x2 - 2.0, 0.0, 1e-8, "parallel dx");
        assert_near(z2, 0.0, 1e-8, "parallel dz");
        let _ = y2;
    }
}

// ── conflict corpus ──────────────────────────────────────────────────────────

#[test]
fn conflict_double_fixed_point() {
    let mut sk = SpaceSketch::new();
    let p = sk.add_point(0.0, 0.0, 0.0);
    sk.add_constraint(SpaceConstraint::fix_point(p, 0.0, 0.0, 0.0));
    sk.add_constraint(SpaceConstraint::fix_point(p, 1.0, 1.0, 1.0));
    let d = sk.diagnose();
    assert!(!d.conflicting.is_empty(), "{d:?}");
    assert_every_entity(&sk, EntityStatus::Conflicting);
}

#[test]
fn conflict_distance_beyond_pinned_span() {
    let mut sk = SpaceSketch::new();
    let p0 = sk.add_point(0.0, 0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0, 0.0);
    sk.add_constraint(SpaceConstraint::fix_point(p0, 0.0, 0.0, 0.0));
    sk.add_constraint(SpaceConstraint::fix_point(p1, 1.0, 0.0, 0.0));
    sk.add_constraint(SpaceConstraint::point_distance(p0, p1, 5.0));
    let d = sk.diagnose();
    assert!(!d.conflicting.is_empty(), "{d:?}");
    assert_every_entity(&sk, EntityStatus::Conflicting);
}

#[test]
fn conflict_line_length_vs_pinned_ends() {
    let mut sk = SpaceSketch::new();
    let l = sk.add_line(0.0, 0.0, 0.0, 3.0, 0.0, 0.0);
    sk.add_constraint(SpaceConstraint::fix_point(
        SpacePointRef::LineStart(l),
        0.0,
        0.0,
        0.0,
    ));
    sk.add_constraint(SpaceConstraint::fix_point(
        SpacePointRef::LineEnd(l),
        3.0,
        0.0,
        0.0,
    ));
    sk.add_constraint(SpaceConstraint::LineLength {
        line: l,
        length: 7.0,
    });
    let d = sk.diagnose();
    assert!(!d.conflicting.is_empty(), "{d:?}");
    assert_every_entity(&sk, EntityStatus::Conflicting);
}

#[test]
fn conflict_parallel_and_perpendicular() {
    let mut sk = SpaceSketch::new();
    let l1 = sk.add_line(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    let l2 = sk.add_line(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    sk.fix_entity(l1);
    sk.fix_entity(l2);
    sk.add_constraint(SpaceConstraint::LineParallelLine(l1, l2));
    sk.add_constraint(SpaceConstraint::LinePerpendicularLine(l1, l2));
    let d = sk.diagnose();
    assert!(!d.conflicting.is_empty(), "{d:?}");
    assert_every_entity(&sk, EntityStatus::Conflicting);
}

#[test]
fn conflict_plane_normals() {
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(0.3, -0.4, 0.5, 1.2);
    sk.add_constraint(SpaceConstraint::PlaneNormal {
        plane,
        nx: 0.0,
        ny: 0.0,
        nz: 1.0,
    });
    sk.add_constraint(SpaceConstraint::PlaneNormal {
        plane,
        nx: 1.0,
        ny: 0.0,
        nz: 0.0,
    });
    let d = sk.diagnose();
    assert!(!d.conflicting.is_empty(), "{d:?}");
    assert_every_entity(&sk, EntityStatus::Conflicting);
}

#[test]
fn conflict_sphere_radius() {
    let mut sk = SpaceSketch::new();
    let s = sk.add_sphere(0.0, 0.0, 0.0, 2.0);
    sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere: s,
        radius: 2.0,
    });
    sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere: s,
        radius: 3.0,
    });
    let d = sk.diagnose();
    assert!(!d.conflicting.is_empty(), "{d:?}");
    assert_every_entity(&sk, EntityStatus::Conflicting);
}

#[test]
fn conflict_cylinder_radius() {
    let mut sk = SpaceSketch::new();
    let c = sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    sk.add_constraint(SpaceConstraint::CylinderRadius {
        cylinder: c,
        radius: 1.0,
    });
    sk.add_constraint(SpaceConstraint::CylinderRadius {
        cylinder: c,
        radius: 2.0,
    });
    let d = sk.diagnose();
    assert!(!d.conflicting.is_empty(), "{d:?}");
    assert_every_entity(&sk, EntityStatus::Conflicting);
}

#[test]
fn conflict_point_pinned_off_plane() {
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(0.0, 0.0, 1.0, 0.0);
    sk.fix_entity(plane);
    let p = sk.add_point(0.0, 0.0, 5.0);
    sk.fix_entity(p);
    sk.add_constraint(SpaceConstraint::PointOnPlane {
        point: SpacePointRef::Point(p),
        plane,
    });
    let d = sk.diagnose();
    assert!(!d.conflicting.is_empty(), "{d:?}");
    assert_every_entity(&sk, EntityStatus::Conflicting);
}

// ── redundancy corpus ────────────────────────────────────────────────────────

fn assert_redundant(sk: &SpaceSketch, label: &str) {
    let d = sk.diagnose();
    assert!(
        !d.redundant.is_empty(),
        "{label}: expected a redundant row, got {d:?}"
    );
    assert!(
        d.conflicting.is_empty(),
        "{label}: misclassified as conflicting: {d:?}"
    );
    // Every entity of these cases carries the duplicated constraint;
    // solved, that is all there is to say about it (as authored it is
    // also unsolved, which outranks over-defined).
    let mut solved = sk.clone();
    assert!(solved.solve().converged, "{label}");
    assert_every_entity(&solved, EntityStatus::OverDefined);
}

/// Every entity of `sk` reports `state` (P4.1).
fn assert_every_entity(sk: &SpaceSketch, state: EntityStatus) {
    let status = sk.status();
    for (e, got) in status.entities.iter().enumerate() {
        assert_eq!(*got, Some(state), "entity {e}: {status:?}");
    }
}

#[test]
fn redundancy_corpus() {
    // 1. Duplicate fixed.
    {
        let mut sk = SpaceSketch::new();
        let p = sk.add_point(0.0, 0.0, 0.0);
        sk.add_constraint(SpaceConstraint::fix_point(p, 3.0, 4.0, 5.0));
        sk.add_constraint(SpaceConstraint::fix_point(p, 3.0, 4.0, 5.0));
        assert_redundant(&sk, "duplicate fixed");
    }

    // 2. Duplicate coincident.
    {
        let mut sk = SpaceSketch::new();
        let p1 = sk.add_point(0.0, 0.0, 0.0);
        let p2 = sk.add_point(5.0, 5.0, 5.0);
        sk.add_constraint(SpaceConstraint::coincident(p1, p2));
        sk.add_constraint(SpaceConstraint::coincident(p1, p2));
        assert_redundant(&sk, "duplicate coincident");
    }

    // 3. Duplicate point-on-plane.
    {
        let mut sk = SpaceSketch::new();
        let p = sk.add_point(0.0, 0.0, 0.0);
        let plane = sk.add_plane(0.0, 0.0, 1.0, 1.0);
        sk.add_constraint(SpaceConstraint::PointOnPlane {
            point: SpacePointRef::Point(p),
            plane,
        });
        sk.add_constraint(SpaceConstraint::PointOnPlane {
            point: SpacePointRef::Point(p),
            plane,
        });
        assert_redundant(&sk, "duplicate point-on-plane");
    }

    // 4. Duplicate sphere radius.
    {
        let mut sk = SpaceSketch::new();
        let s = sk.add_sphere(0.0, 0.0, 0.0, 2.0);
        sk.add_constraint(SpaceConstraint::SphereRadius {
            sphere: s,
            radius: 2.0,
        });
        sk.add_constraint(SpaceConstraint::SphereRadius {
            sphere: s,
            radius: 2.0,
        });
        assert_redundant(&sk, "duplicate sphere radius");
    }

    // 5. Duplicate cylinder radius.
    {
        let mut sk = SpaceSketch::new();
        let c = sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
        sk.add_constraint(SpaceConstraint::CylinderRadius {
            cylinder: c,
            radius: 2.0,
        });
        sk.add_constraint(SpaceConstraint::CylinderRadius {
            cylinder: c,
            radius: 2.0,
        });
        assert_redundant(&sk, "duplicate cylinder radius");
    }

    // 6. Duplicate line length.
    {
        let mut sk = SpaceSketch::new();
        let l = sk.add_line(0.0, 0.0, 0.0, 3.0, 4.0, 0.0);
        sk.add_constraint(SpaceConstraint::LineLength {
            line: l,
            length: 5.0,
        });
        sk.add_constraint(SpaceConstraint::LineLength {
            line: l,
            length: 5.0,
        });
        assert_redundant(&sk, "duplicate line length");
    }

    // 7. Duplicate perpendicular.
    {
        let mut sk = SpaceSketch::new();
        let l1 = sk.add_line(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        let l2 = sk.add_line(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        sk.add_constraint(SpaceConstraint::LinePerpendicularLine(l1, l2));
        sk.add_constraint(SpaceConstraint::LinePerpendicularLine(l1, l2));
        assert_redundant(&sk, "duplicate perpendicular");
    }

    // 8. Duplicate point distance.
    {
        let mut sk = SpaceSketch::new();
        let p1 = sk.add_point(0.0, 0.0, 0.0);
        let p2 = sk.add_point(3.0, 4.0, 0.0);
        sk.add_constraint(SpaceConstraint::point_distance(p1, p2, 5.0));
        sk.add_constraint(SpaceConstraint::point_distance(p1, p2, 5.0));
        assert_redundant(&sk, "duplicate point distance");
    }
}

#[test]
fn dependency_group_names_the_whole_relation() {
    let mut sk = SpaceSketch::new();
    let s = sk.add_sphere(0.0, 0.0, 0.0, 2.0);
    sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere: s,
        radius: 2.0,
    });
    sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere: s,
        radius: 2.0,
    });
    let d = sk.diagnose();
    assert!(
        d.dependency_groups.iter().any(|g| g == &vec![0, 1]),
        "duplicate rows form one dependency group: {d:?}"
    );
}

// ── DOF (rank-based) corpus ──────────────────────────────────────────────────

#[test]
fn rank_based_dof_corpus() {
    let mut sk = SpaceSketch::new();
    let p = sk.add_point(1.0, 2.0, 3.0);
    assert_eq!(sk.diagnose().dof, 3, "free point has three DOF");
    assert_every_entity(&sk, EntityStatus::UnderDefined);
    sk.add_constraint(SpaceConstraint::fix_point(p, 1.0, 2.0, 3.0));
    assert_eq!(sk.diagnose().dof, 0, "fixed point has zero DOF");
    assert_eq!(sk.diagnose().rank, 3, "fixed point has rank three");
    assert_every_entity(&sk, EntityStatus::WellDefined);

    let mut sk = SpaceSketch::new();
    let a = sk.add_point(0.0, 0.0, 0.0);
    let b = sk.add_point(1.0, 0.0, 0.0);
    sk.add_constraint(SpaceConstraint::point_distance(a, b, 2.0));
    let d = sk.diagnose();
    assert_eq!(
        (d.rank, d.dof),
        (1, 5),
        "one distance removes one of six DOF: {d:?}"
    );

    let mut sk = SpaceSketch::new();
    sk.add_line(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    assert_eq!(sk.diagnose().dof, 6, "free line has six DOF");
    assert_every_entity(&sk, EntityStatus::UnderDefined);

    let mut sk = SpaceSketch::new();
    sk.add_sphere(0.0, 0.0, 0.0, 1.0);
    assert_eq!(sk.diagnose().dof, 4, "free sphere has four DOF");
    assert_every_entity(&sk, EntityStatus::UnderDefined);

    // A free plane has four storage coordinates but only three physical
    // DOF: the (normal, offset) scale direction is gauge.
    let mut sk = SpaceSketch::new();
    sk.add_plane(0.3, -0.4, 0.5, 1.2);
    assert_eq!(sk.diagnose().dof, 3, "free plane has three DOF");
    assert_every_entity(&sk, EntityStatus::UnderDefined);

    // Plane with a fixed offset: the gauge orbit is unavailable, the normal's
    // three coordinates are all physical.
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(0.3, -0.4, 0.5, 1.2);
    let start = sk.entities[plane].param_start;
    sk.fix_param(start + 3);
    assert_eq!(
        sk.diagnose().dof,
        3,
        "offset-fixed plane keeps three normal DOF"
    );
    assert_every_entity(&sk, EntityStatus::UnderDefined);
}

// ── determinism ──────────────────────────────────────────────────────────────

fn deterministic_system() -> SpaceSketch {
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(0.3, -0.4, 0.5, 1.2);
    let p = sk.add_point(0.0, 0.0, 4.0);
    sk.fix_entity(p);
    sk.add_constraint(SpaceConstraint::PlaneNormal {
        plane,
        nx: 0.0,
        ny: 0.0,
        nz: 1.0,
    });
    sk.add_constraint(SpaceConstraint::PointOnPlane {
        point: SpacePointRef::Point(p),
        plane,
    });
    sk
}

#[test]
fn repeated_solve_is_bitwise_deterministic() {
    let mut a = deterministic_system();
    let mut b = deterministic_system();
    solve(&mut a);
    solve(&mut b);
    assert_eq!(
        a.params_snapshot(),
        b.params_snapshot(),
        "same input must give identical bits"
    );
}

#[test]
fn re_solve_at_the_solution_is_bitwise_stable() {
    let mut sk = deterministic_system();
    solve(&mut sk);
    let first = sk.params_snapshot();
    solve(&mut sk);
    assert_eq!(
        sk.params_snapshot(),
        first,
        "warm start at the solution must not move"
    );
}

#[test]
fn per_component_solve_matches_all_at_once() {
    fn build() -> SpaceSketch {
        let mut sk = SpaceSketch::new();
        let pa = sk.add_point(0.0, 0.0, 0.0);
        let pb = sk.add_point(4.0, 0.0, 0.0);
        sk.add_constraint(SpaceConstraint::fix_point(pa, 0.0, 0.0, 0.0));
        sk.add_constraint(SpaceConstraint::point_distance(pa, pb, 5.0));
        let s = sk.add_sphere(9.0, 9.0, 9.0, 2.0);
        sk.add_constraint(SpaceConstraint::SphereRadius {
            sphere: s,
            radius: 2.5,
        });
        sk
    }

    let mut all = build();
    solve(&mut all);

    let mut piecewise = build();
    let comps = piecewise.components();
    assert_eq!(comps.len(), 2, "two independent groups: {comps:?}");
    for c in 0..comps.len() {
        let r = piecewise.solve_selected(Some(&[c]));
        assert!(r.converged, "component {c}: residual {}", r.residual);
    }
    for (a, b) in all.params.iter().zip(piecewise.params.iter()) {
        assert!(
            (a - b).abs() < 1e-12,
            "per-component params diverge: {a} vs {b}"
        );
    }
}

#[test]
fn warm_start_from_drift_converges_to_the_same_fixed_point() {
    let mut clean = deterministic_system();
    solve(&mut clean);
    let target = clean.params_snapshot();

    let mut drifted = deterministic_system();
    drifted.drift_free_params(0.35, 0x5EED);
    solve(&mut drifted);
    for (i, (a, b)) in target.iter().zip(drifted.params.iter()).enumerate() {
        assert!(
            (a - b).abs() < 1e-8,
            "warm start changed the fixed point at param {i}: {a} vs {b}"
        );
    }
}

// ── incremental (component-local) solve ─────────────────────────────────────

fn incremental_system() -> SpaceSketch {
    let mut sk = SpaceSketch::new();
    let pa = sk.add_point(0.0, 0.0, 0.0);
    let pb = sk.add_point(4.0, 0.0, 0.0);
    sk.add_constraint(SpaceConstraint::fix_point(pa, 0.0, 0.0, 0.0));
    sk.add_constraint(SpaceConstraint::point_distance(pa, pb, 5.0));
    let s = sk.add_sphere(9.0, 9.0, 9.0, 2.0);
    sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere: s,
        radius: 2.5,
    });
    sk
}

#[test]
fn incremental_resolve_leaves_other_components_bitwise_untouched() {
    let mut sk = incremental_system();
    solve(&mut sk);
    let sphere_start = sk.entities[2].param_start;
    let before_sphere: Vec<f64> = sk.params[sphere_start..sphere_start + 4].to_vec();

    // Drift group A's point and retarget its driving distance (constraint 1).
    let pb_start = sk.entities[1].param_start;
    sk.params[pb_start] += 0.8;
    sk.params[pb_start + 1] += 0.5;
    sk.constraints[1] = SpaceConstraint::PointDistance {
        p1: SpacePointRef::Point(0),
        p2: SpacePointRef::Point(1),
        distance: 3.0,
    };

    let r = sk.solve_incremental(&[1]);
    assert!(r.converged, "incremental residual {}", r.residual);

    // Group B (sphere) parameters are bitwise unchanged.
    for (k, before) in before_sphere.iter().enumerate() {
        assert_eq!(
            sk.params[sphere_start + k],
            *before,
            "untouched component param {} moved",
            sphere_start + k
        );
    }
    // Group A actually reached the new distance.
    let p0 = point3(&sk, 0);
    let p1 = point3(&sk, 1);
    assert_near(dist(p0, p1), 3.0, 1e-8, "incremental group A distance");
}

#[test]
fn incremental_selection_matches_a_full_solve() {
    let mut incremental = incremental_system();
    solve(&mut incremental);
    let mut full = incremental.clone();

    let pb_start = incremental.entities[1].param_start;
    for sk in [&mut incremental, &mut full] {
        sk.params[pb_start] += 0.8;
        sk.params[pb_start + 1] += 0.5;
        sk.constraints[1] = SpaceConstraint::PointDistance {
            p1: SpacePointRef::Point(0),
            p2: SpacePointRef::Point(1),
            distance: 3.0,
        };
    }
    let r_inc = incremental.solve_incremental(&[1]);
    assert!(r_inc.converged, "incremental residual {}", r_inc.residual);
    let r_full = full.solve();
    assert!(r_full.converged, "full residual {}", r_full.residual);

    for (a, b) in incremental.params.iter().zip(full.params.iter()) {
        assert!((a - b).abs() < 1e-12, "incremental vs full: {a} vs {b}");
    }
}

#[test]
fn empty_incremental_selection_does_not_move() {
    let mut sk = incremental_system();
    solve(&mut sk);
    let before = sk.params_snapshot();
    let r = sk.solve_incremental(&[]);
    assert!(r.converged, "settled system should report converged: {r:?}");
    assert_eq!(
        sk.params_snapshot(),
        before,
        "empty selection must not move params"
    );
}

#[test]
fn solve_selected_leaves_unselected_components_bitwise_untouched() {
    let mut sk = incremental_system();
    solve(&mut sk);
    let sphere_start = sk.entities[2].param_start;
    let before: Vec<f64> = sk.params[sphere_start..sphere_start + 4].to_vec();

    let pb_start = sk.entities[1].param_start;
    sk.params[pb_start] += 0.7;
    let comps = sk.components();
    let a_ord = comps.iter().position(|g| g.contains(&1)).expect("group A");
    let r = sk.solve_selected(Some(&[a_ord]));
    assert!(r.converged, "residual {}", r.residual);
    for (k, v) in before.iter().enumerate() {
        assert_eq!(
            sk.params[sphere_start + k],
            *v,
            "sphere param {} moved",
            sphere_start + k
        );
    }
}

// ── plane parameterization hardening ────────────────────────────────────────

#[test]
fn plane_insertion_pins_the_unit_normal_representation() {
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(3.0, 4.0, 0.0, 10.0);
    let e = &sk.entities[plane];
    assert_near(sk.params[e.param(0)], 0.6, 1e-15, "nx");
    assert_near(sk.params[e.param(1)], 0.8, 1e-15, "ny");
    assert_near(sk.params[e.param(2)], 0.0, 1e-15, "nz");
    assert_near(sk.params[e.param(3)], 2.0, 1e-15, "d");
}

#[test]
fn equivalent_plane_scalings_produce_identical_systems() {
    let mut a = SpaceSketch::new();
    let pa = a.add_plane(0.0, 0.0, 1.0, 10.0);
    let mut b = SpaceSketch::new();
    let pb = b.add_plane(0.0, 0.0, 5.0, 50.0);
    assert_eq!(
        a.params, b.params,
        "scaled representations must normalize identically"
    );
    let _ = (pa, pb);
    assert_eq!(
        a.diagnose(),
        b.diagnose(),
        "diagnosis must be scale invariant"
    );
}

#[test]
fn fixed_plane_parameters_are_not_renormalized() {
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(0.0, 0.0, 1.0, 2.0);
    let start = sk.entities[plane].param_start;
    // Deliberately de-normalize, then mask one coordinate: the gauge fix
    // must leave the whole plane alone (normalization would move a masked
    // value).
    sk.params[start] = 2.0;
    sk.params[start + 3] = 8.0;
    sk.fix_param(start);
    let before = sk.params_snapshot();
    sk.normalize_planes();
    assert_eq!(
        sk.params_snapshot(),
        before,
        "masked plane was renormalized"
    );
}

#[test]
fn solve_keeps_free_planes_unit_normal() {
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(6.0, -3.0, 0.0, 2.0);
    let p = sk.add_point(0.0, 0.0, 0.0);
    sk.fix_entity(p);
    sk.add_constraint(SpaceConstraint::PointOnPlane {
        point: SpacePointRef::Point(p),
        plane,
    });
    solve(&mut sk);
    let n = plane_normal(&sk, plane);
    let len = (n.0 * n.0 + n.1 * n.1 + n.2 * n.2).sqrt();
    assert_near(len, 1.0, 1e-12, "plane normal unit length after solve");
    let d = sk.params[sk.entities[plane].param(3)];
    assert_near(d, 0.0, 1e-8, "point at origin lies on the plane");
}

// ── performance evidence (ignored by default) ───────────────────────────────

/// System-level timing on a 150-point chain (~450 parameters): diagnosis,
/// full solve, and the no-op incremental path. The fair analytic-vs-finite-
/// difference assembly comparison lives in the library test
/// `space3d::solver::tests::analytic_jacobian_assembly_bench` (it needs the
/// crate-private assembly entry point). Run this one with:
/// `cargo test -p rcad-gcs --test space3d_corpus -- --ignored --nocapture space3d_system_bench`
#[test]
#[ignore]
fn space3d_system_bench() {
    let n = 150usize;
    let mut sk = SpaceSketch::new();
    for i in 0..n {
        sk.add_point(i as f64, (i % 7) as f64 * 0.1, 0.0);
    }
    sk.add_constraint(SpaceConstraint::fix_point(0, 0.0, 0.0, 0.0));
    for i in 0..n - 1 {
        sk.add_constraint(SpaceConstraint::point_distance(i, i + 1, 1.0));
    }

    let t0 = Instant::now();
    let d = sk.diagnose();
    let diagnose_us = t0.elapsed().as_secs_f64() * 1e6;
    println!(
        "diagnose (analytic assembly + QR): {diagnose_us:.0} us, rank {}, dof {}",
        d.rank, d.dof
    );

    let t1 = Instant::now();
    let r = sk.solve();
    let full_us = t1.elapsed().as_secs_f64() * 1e6;
    println!(
        "full solve: {full_us:.0} us (converged {}, iterations {})",
        r.converged, r.iterations
    );

    let t2 = Instant::now();
    let _ = sk.solve_incremental(&[]);
    let incremental_us = t2.elapsed().as_secs_f64() * 1e6;
    println!("empty incremental selection: {incremental_us:.0} us");
}
