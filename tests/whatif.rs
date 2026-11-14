//! What-if probes (P4.5) against the full diagnosis. The corpus gates
//! (`would_add_matches_the_diagnosis` in `tests/corpus.rs` and
//! `tests/assembly_corpus.rs`) compare every answer on the 2D and
//! assembly corpora; these cover the 3D sketch (its gauge: plane normal
//! length, cylinder axis slide) and the paths the corpora seldom take.

use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::space3d::{SpaceConstraint, SpaceEntityKind, SpacePointRef, SpaceSketch};
use rcad_gcs::{Addition, Constraint, Diagnosis, Sketch, Verdict};

/// What the diagnosis says about adding a constraint of `k` equations,
/// from the diagnoses before and after.
fn diagnosed(before: &Diagnosis, after: &Diagnosis, k: usize) -> Addition {
    let grew = |now: &[usize], then: &[usize]| now.iter().any(|x| !then.contains(x));
    let verdict = if after.rank as i64 - before.rank as i64 == k as i64 {
        Verdict::Independent
    } else if grew(&after.conflicting, &before.conflicting) {
        Verdict::Conflicting
    } else if grew(&after.partially_redundant, &before.partially_redundant) {
        Verdict::PartiallyRedundant
    } else {
        Verdict::Redundant
    };
    Addition {
        verdict,
        removes: (before.dof - after.dof).max(0) as usize,
    }
}

fn diagnosed_2d(sk: &Sketch, c: &Constraint) -> Addition {
    let mut with = sk.clone();
    with.add_constraint(c.clone());
    diagnosed(&sk.diagnose(), &with.diagnose(), c.equation_count())
}

fn diagnosed_3d(sk: &SpaceSketch, c: &SpaceConstraint) -> Addition {
    let mut with = sk.clone();
    with.add_constraint(c.clone());
    diagnosed(&sk.diagnose(), &with.diagnose(), c.equation_count())
}

// ── 2D paths ────────────────────────────────────────────────────────────────

#[test]
fn a_candidate_that_joins_two_components_or_touches_only_held_geometry() {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(4.0, 0.0);
    let held = sk.add_point(9.0, 9.0);
    sk.fix_entity(held);
    sk.add_constraint(Constraint::fix_point(a, 0.0, 0.0));
    sk.add_constraint(DistanceY {
        p1: PointRef::Point(b),
        p2: PointRef::Point(a),
        distance: 0.0,
    });
    assert!(sk.solve().converged);
    let probe = sk.what_if();
    let candidates = [
        // Joins the two: b's x is free.
        Constraint::point_distance(a, b, 4.0),
        Constraint::point_distance(a, b, 5.0),
        // b's y again, agreeing and not.
        Fixed {
            point: PointRef::Point(b),
            x: 4.0,
            y: 0.0,
        },
        Fixed {
            point: PointRef::Point(b),
            x: 4.0,
            y: 1.0,
        },
        // Only held geometry: met or not.
        Fixed {
            point: PointRef::Point(held),
            x: 9.0,
            y: 9.0,
        },
        Fixed {
            point: PointRef::Point(held),
            x: 9.5,
            y: 9.0,
        },
        // Reaches the held point from b.
        Constraint::point_distance(b, held, 9.85),
    ];
    for c in &candidates {
        assert_eq!(probe.would_add(c), diagnosed_2d(&sk, c), "{c:?}");
    }
    // A Fixed on b where its y is already held takes one freedom, not two.
    assert_eq!(
        probe.would_add(&candidates[2]),
        Addition {
            verdict: Verdict::Redundant,
            removes: 1
        }
    );
    assert_eq!(
        probe.would_add(&candidates[3]).verdict,
        Verdict::Conflicting
    );
    // The probe leaves the sketch as it was.
    drop(probe);
    assert_eq!(sk.constraints_len(), 2);
}

// ── 3D ──────────────────────────────────────────────────────────────────────

/// A plane stored with a normal of length ≠ 1, points on it, a line, a
/// sphere and a cylinder, partly constrained, solved.
fn space_fixture() -> SpaceSketch {
    let mut sk = SpaceSketch::new();
    let plane = sk.add_plane(0.1, 0.2, 2.0, 0.6);
    let p0 = sk.add_point(0.0, 0.0, 0.0);
    let p1 = sk.add_point(3.0, 0.2, 0.1);
    let line = sk.add_line(0.0, 0.0, 0.0, 0.1, 0.2, 2.0);
    let sphere = sk.add_sphere(5.0, 5.0, 5.0, 1.0);
    let cylinder = sk.add_cylinder(0.0, 4.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    let p2 = sk.add_point(1.0, 4.0, 0.5);
    sk.add_constraint(SpaceConstraint::fix_point(p0, 0.0, 0.0, 0.0));
    sk.add_constraint(SpaceConstraint::PlaneNormal {
        plane,
        nx: 0.0,
        ny: 0.0,
        nz: 1.0,
    });
    sk.add_constraint(SpaceConstraint::PointOnPlane {
        point: SpacePointRef::Point(p0),
        plane,
    });
    sk.add_constraint(SpaceConstraint::point_distance(p0, p1, 3.0));
    sk.add_constraint(SpaceConstraint::PointOnPlane {
        point: SpacePointRef::Point(p1),
        plane,
    });
    sk.add_constraint(SpaceConstraint::coincident(
        SpacePointRef::LineStart(line),
        p0,
    ));
    sk.add_constraint(SpaceConstraint::LineLength { line, length: 2.0 });
    sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere,
        radius: 1.0,
    });
    sk.add_constraint(SpaceConstraint::CylinderRadius {
        cylinder,
        radius: 1.0,
    });
    sk.add_constraint(SpaceConstraint::PointOnCylinder {
        point: SpacePointRef::Point(p2),
        cylinder,
    });
    assert!(sk.solve().converged);
    sk
}

/// Candidates of every kind on the fixture's entities, dimensions at what
/// they measure and off it.
fn space_candidates(sk: &SpaceSketch) -> Vec<SpaceConstraint> {
    let mut points = Vec::new();
    let (mut planes, mut lines, mut spheres, mut cylinders) = (vec![], vec![], vec![], vec![]);
    for (e, entity) in sk.entities.iter().enumerate() {
        match entity.kind {
            SpaceEntityKind::SpacePoint => points.push(SpacePointRef::Point(e)),
            SpaceEntityKind::SpaceLine => {
                points.push(SpacePointRef::LineStart(e));
                points.push(SpacePointRef::LineEnd(e));
                lines.push(e);
            }
            SpaceEntityKind::Plane => planes.push(e),
            SpaceEntityKind::Sphere => {
                points.push(SpacePointRef::SphereCenter(e));
                spheres.push(e);
            }
            SpaceEntityKind::Cylinder => cylinders.push(e),
            _ => {}
        }
    }
    let measured = |c: &SpaceConstraint, off: f64| {
        let mut c = c.clone();
        let m = c.measure(&sk.params, &sk.entities).expect("a dimension");
        c.set_value(m + off);
        c
    };
    let mut out = Vec::new();
    for &p in &points {
        let (x, y, z) = p.point_coords(&sk.params, &sk.entities);
        out.push(SpaceConstraint::fix_point(p, x, y, z));
        out.push(SpaceConstraint::fix_point(p, x + 0.5, y, z));
        for &plane in &planes {
            out.push(SpaceConstraint::PointOnPlane { point: p, plane });
        }
        for &sphere in &spheres {
            out.push(SpaceConstraint::PointOnSphere { point: p, sphere });
        }
        for &cylinder in &cylinders {
            out.push(SpaceConstraint::PointOnCylinder { point: p, cylinder });
        }
    }
    for (i, &a) in points.iter().enumerate() {
        for &b in &points[i + 1..] {
            out.push(SpaceConstraint::coincident(a, b));
            for off in [0.0, 0.7] {
                out.push(measured(&SpaceConstraint::point_distance(a, b, 1.0), off));
            }
        }
    }
    for &plane in &planes {
        let start = sk.entities[plane].param_start;
        let n = &sk.params[start..start + 3];
        out.push(SpaceConstraint::PlaneNormal {
            plane,
            nx: n[0],
            ny: n[1],
            nz: n[2],
        });
        out.push(SpaceConstraint::PlaneNormal {
            plane,
            nx: 1.0,
            ny: 0.0,
            nz: 0.0,
        });
    }
    for (i, &l) in lines.iter().enumerate() {
        for off in [0.0, 1.0] {
            out.push(measured(
                &SpaceConstraint::LineLength {
                    line: l,
                    length: 1.0,
                },
                off,
            ));
        }
        for &m in &lines[i + 1..] {
            out.push(SpaceConstraint::LineParallelLine(l, m));
            out.push(SpaceConstraint::LinePerpendicularLine(l, m));
        }
    }
    for &s in &spheres {
        for off in [0.0, 0.5] {
            out.push(measured(
                &SpaceConstraint::SphereRadius {
                    sphere: s,
                    radius: 1.0,
                },
                off,
            ));
        }
    }
    for &c in &cylinders {
        for off in [0.0, 0.3] {
            out.push(measured(
                &SpaceConstraint::CylinderRadius {
                    cylinder: c,
                    radius: 1.0,
                },
                off,
            ));
        }
    }
    out
}

#[test]
fn space_would_add_matches_the_diagnosis() {
    let mut free = SpaceSketch::new();
    free.add_point(1.0, 2.0, 3.0);
    free.add_line(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    free.add_plane(0.3, -0.4, 0.5, 1.2);
    free.add_cylinder(1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.5);
    free.add_sphere(0.0, 3.0, 0.0, 1.0);
    let (mut total, mut right) = (0, 0);
    for sk in [space_fixture(), free] {
        let probe = sk.what_if();
        for c in space_candidates(&sk) {
            total += 1;
            let got = probe.would_add(&c);
            let want = diagnosed_3d(&sk, &c);
            if got == want {
                right += 1;
            } else {
                println!("{c:?}: would_add {got:?}, diagnosis {want:?}");
            }
        }
    }
    println!("3D would_add: {right}/{total} as the diagnosis");
    assert_eq!(right, total);
    assert!(total > 100, "{total}");
}

#[test]
fn pinning_a_cylinder_axis_point_takes_three_ranks_and_two_freedoms() {
    // A cylinder through a held point: its axis point sliding along the
    // axis is how it is stored, and nothing pins it. Making it concentric
    // with a held sphere pins all three coordinates of that point (three
    // ranks) but takes only the two freedoms across the axis: the slide
    // was never one. (The cylinder is in a component already, so this is
    // the incremental test.)
    let mut sk = SpaceSketch::new();
    let cylinder = sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    let sphere = sk.add_sphere(0.0, 0.0, 0.0, 2.0);
    sk.fix_entity(sphere);
    let p = sk.add_point(1.0, 0.0, 0.0);
    sk.fix_entity(p);
    sk.add_constraint(SpaceConstraint::PointOnCylinder {
        point: SpacePointRef::Point(p),
        cylinder,
    });
    let probe = sk.what_if();
    let pin = SpaceConstraint::Concentric3D(cylinder, sphere);
    let got = probe.would_add(&pin);
    assert_eq!(
        got,
        Addition {
            verdict: Verdict::Independent,
            removes: 2
        }
    );
    assert_eq!(got, diagnosed_3d(&sk, &pin));
}
