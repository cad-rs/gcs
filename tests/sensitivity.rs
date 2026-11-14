//! Driven dimensions and their sensitivities (P4.7). The corpus gate
//! (`sensitivities_match_finite_differences` in `tests/corpus.rs`) checks
//! every derivative on the 2D corpus against finite differences, and that
//! a driven dimension has none exactly when adding it would be
//! independent; these cover the cases named in the docs, the 3D sketch and
//! assemblies.

use glam::{DQuat, DVec3};
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::rigid3d::{ConstraintSolver, ConstraintSystem, PoseSolveOptions};
use rcad_gcs::space3d::{SpaceConstraint, SpaceSketch};
use rcad_gcs::{Constraint, Sketch};

/// A 4 × 3 rectangle pinned at its origin, its width and height driven
/// (when `sized`); returns the sketch, the origin, the far corner and the
/// indices of the width and height.
fn rectangle(sized: bool) -> (Sketch, PointRef, PointRef, usize, usize) {
    let mut sk = Sketch::new();
    let o = sk.add_point(0.0, 0.0);
    let bottom = sk.add_line(0.0, 0.0, 4.0, 0.0);
    let right = sk.add_line(4.0, 0.0, 4.0, 3.0);
    let top = sk.add_line(4.0, 3.0, 0.0, 3.0);
    let left = sk.add_line(0.0, 3.0, 0.0, 0.0);
    sk.add_constraint(Constraint::fix_point(o, 0.0, 0.0));
    sk.add_constraint(Coincident(PointRef::LineStart(bottom), PointRef::Point(o)));
    for (a, b) in [(bottom, right), (right, top), (top, left), (left, bottom)] {
        sk.add_constraint(Coincident(PointRef::LineEnd(a), PointRef::LineStart(b)));
    }
    sk.add_constraint(Horizontal(bottom));
    sk.add_constraint(Horizontal(top));
    sk.add_constraint(Vertical(right));
    sk.add_constraint(Vertical(left));
    let width = sk.add_constraint(LineLength {
        line: bottom,
        length: 4.0,
    });
    let height = if sized {
        sk.add_constraint(LineLength {
            line: right,
            length: 3.0,
        })
    } else {
        usize::MAX
    };
    assert!(sk.solve().converged);
    (
        sk,
        PointRef::Point(o),
        PointRef::LineEnd(right),
        width,
        height,
    )
}

#[test]
fn a_rectangle_diagonal_follows_its_width_and_height() {
    let (sk, o, corner, width, height) = rectangle(true);
    let diagonal = Constraint::point_distance(o, corner, 0.0);
    assert!((sk.measure_constraint(&diagonal).unwrap() - 5.0).abs() < 1e-12);
    let s = sk.sensitivities(&diagonal).expect("determined");
    assert_eq!(s.len(), 2, "{s:?}");
    assert_eq!((s[0].0, s[1].0), (width, height));
    assert!(
        (s[0].1 - 0.8).abs() < 1e-9 && (s[1].1 - 0.6).abs() < 1e-9,
        "{s:?}"
    );
}

#[test]
fn an_undetermined_driven_dimension_has_no_sensitivity() {
    // Without a height the diagonal can change on its own.
    let (sk, o, corner, _, _) = rectangle(false);
    let diagonal = Constraint::point_distance(o, corner, 0.0);
    assert_eq!(sk.sensitivities(&diagonal), None);
    // The width it does follow: one to one.
    let along = DistanceX {
        p1: corner,
        p2: o,
        distance: 0.0,
    };
    let s = sk.sensitivities(&along).expect("determined");
    assert_eq!(s.len(), 1);
    assert!((s[0].1 - 1.0).abs() < 1e-9, "{s:?}");
    // Not a single value: no sensitivity.
    assert_eq!(sk.sensitivities(&Horizontal(1)), None);
}

#[test]
fn redundant_driving_dimensions_share_the_derivative() {
    // The width given twice: they cannot change alone; the least-squares
    // split of smallest norm gives each half.
    let (mut sk, o, corner, width, _) = rectangle(true);
    let again = sk.add_constraint(LineLength {
        line: 1,
        length: 4.0,
    });
    let along = DistanceX {
        p1: corner,
        p2: o,
        distance: 0.0,
    };
    let s = sk.sensitivities(&along).expect("determined");
    assert_eq!(
        s.iter().map(|&(c, _)| c).collect::<Vec<_>>(),
        vec![width, again]
    );
    assert!(s.iter().all(|&(_, d)| (d - 0.5).abs() < 1e-9), "{s:?}");
}

#[test]
fn a_dimension_between_held_points_depends_on_nothing() {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(3.0, 4.0);
    sk.fix_entity(a);
    sk.fix_entity(b);
    let d = Constraint::point_distance(a, b, 0.0);
    assert_eq!(sk.measure_constraint(&d), Some(5.0));
    assert_eq!(sk.sensitivities(&d), Some(Vec::new()));
}

// ── 3D ──────────────────────────────────────────────────────────────────────

#[test]
fn space_sensitivities_match_finite_differences() {
    // A fixed point, a point at a driven distance on a plane with a driven
    // normal, a line of driven length from the first, a sphere of driven
    // radius at a driven distance.
    let mut sk = SpaceSketch::new();
    let p0 = sk.add_point(0.0, 0.0, 0.0);
    let p1 = sk.add_point(3.0, 1.0, 0.0);
    let plane = sk.add_plane(0.0, 0.0, 1.0, 0.0);
    sk.fix_entity(plane);
    let s = sk.add_sphere(2.0, 2.0, 2.0, 0.5);
    sk.add_constraint(SpaceConstraint::fix_point(p0, 0.0, 0.0, 0.0));
    sk.add_constraint(SpaceConstraint::point_distance(p0, p1, 3.2));
    sk.add_constraint(SpaceConstraint::PointOnPlane {
        point: p1.into(),
        plane,
    });
    sk.add_constraint(SpaceConstraint::PointDistanceAlongDir {
        p1: p1.into(),
        p2: p0.into(),
        dx: 0.0,
        dy: 1.0,
        dz: 0.0,
        distance: 1.0,
    });
    sk.add_constraint(SpaceConstraint::SphereRadius {
        sphere: s,
        radius: 0.5,
    });
    sk.add_constraint(SpaceConstraint::fix_point(
        rcad_gcs::space3d::SpacePointRef::SphereCenter(s),
        2.0,
        2.0,
        2.0,
    ));
    assert!(sk.solve().converged);
    let driven = [
        SpaceConstraint::point_distance(p1, rcad_gcs::space3d::SpacePointRef::SphereCenter(s), 0.0),
        SpaceConstraint::PointDistanceAlongDir {
            p1: p1.into(),
            p2: p0.into(),
            dx: 1.0,
            dy: 0.0,
            dz: 0.0,
            distance: 0.0,
        },
    ];
    for d in &driven {
        let got = sk.sensitivities(d).expect("determined");
        for j in 0..sk.constraints.len() {
            let Some(v) = sk.constraints[j].value() else {
                continue;
            };
            let h = 1e-5;
            let moved = |dv: f64| {
                let mut t = sk.clone();
                t.set_constraint_value(j, v + dv);
                assert!(t.solve().converged);
                t.measure_constraint(d).unwrap()
            };
            let fd = (moved(h) - moved(-h)) / (2.0 * h);
            let want = got.iter().find(|(c, _)| *c == j).map_or(0.0, |&(_, x)| x);
            assert!((fd - want).abs() < 1e-5, "{d:?} by {j}: {want} vs {fd}");
        }
    }
    // The sphere's own distance from p1 depends on the p0–p1 distance.
    let got = sk.sensitivities(&driven[0]).unwrap();
    assert!(got.iter().any(|&(c, _)| c == 1), "{got:?}");
}

// ── assemblies ──────────────────────────────────────────────────────────────

#[test]
fn assembly_sensitivities_match_finite_differences() {
    // A grounded body and a second body on a ball joint at (2, 0, 0), held
    // by two driven distances and a driven angle.
    use rcad_gcs::rigid3d::Constraint as Mate;
    let mut system = ConstraintSystem::with_components(2);
    system.add_constraint(Mate::fixed_with_orientation(
        0,
        0,
        DVec3::ZERO,
        DQuat::IDENTITY,
    ));
    system.add_constraint(Mate::coincident(
        1,
        0,
        DVec3::new(2.0, 0.0, 0.0),
        1,
        DVec3::ZERO,
    ));
    system.add_constraint(Mate::distance(
        2,
        0,
        DVec3::new(0.0, 2.0, 0.0),
        1,
        DVec3::new(1.0, 0.0, 0.0),
        2.9,
    ));
    system.add_constraint(Mate::distance(
        3,
        0,
        DVec3::new(0.0, 0.0, 2.0),
        1,
        DVec3::new(0.0, 1.0, 0.0),
        3.1,
    ));
    system.add_constraint(Mate::angle(4, 0, DVec3::Z, 1, DVec3::Z, 0.3));
    let solver = ConstraintSolver::with_options(PoseSolveOptions {
        tolerance: 1e-13,
        ..PoseSolveOptions::default()
    });
    // The values are what an authored pose measures, so it is a solution.
    let pos = vec![DVec3::ZERO, DVec3::new(2.0, 0.0, 0.0)];
    let rot = vec![
        DQuat::IDENTITY,
        DQuat::from_euler(glam::EulerRot::XYZ, 0.2, -0.3, 0.5),
    ];
    for j in 2..5 {
        let m = system.constraints[j].kind.measure(&pos, &rot).unwrap();
        system.set_constraint_value(j, m);
    }
    let r = solver.solve_from(&system, pos.clone(), rot.clone());
    assert!(r.success);
    let (pos, rot) = (r.positions, r.rotations);
    let driven = Mate::distance(
        9,
        0,
        DVec3::new(-1.0, 0.0, 0.0),
        1,
        DVec3::new(0.5, 0.5, 0.5),
        0.0,
    );
    let got = rcad_gcs::rigid3d::assembly_sensitivities(&system, &pos, &rot, &driven)
        .expect("determined");
    assert!(!got.is_empty());
    for j in 2..5 {
        let v = system.constraints[j].kind.value().unwrap();
        let h = 1e-5;
        let moved = |dv: f64| {
            let mut s = system.clone();
            s.set_constraint_value(j, v + dv);
            let r = solver.solve_from(&s, pos.clone(), rot.clone());
            assert!(r.success);
            driven.kind.measure(&r.positions, &r.rotations).unwrap()
        };
        let fd = (moved(h) - moved(-h)) / (2.0 * h);
        let want = got.iter().find(|(c, _)| *c == j).map_or(0.0, |&(_, x)| x);
        assert!((fd - want).abs() < 1e-5, "by {j}: {want} vs {fd}");
    }
    // Without the angle the body can still spin: not determined.
    let mut loose = system.clone();
    loose.set_constraint_enabled(4, false);
    assert_eq!(
        rcad_gcs::rigid3d::assembly_sensitivities(&loose, &pos, &rot, &driven),
        None
    );
}
