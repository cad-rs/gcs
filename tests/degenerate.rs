//! Degenerate geometry: zero-length lines, zero radii, coincident points,
//! and numbers that are not numbers.
//!
//! `report_degenerate` prints what a solve and a diagnosis do on each case;
//! it is a measurement and ignored by default:
//!
//! ```text
//! cargo test --test degenerate -- --ignored --nocapture
//! ```

use rcad_gcs::constraint::Constraint;
use rcad_gcs::entity::PointRef;
use rcad_gcs::Sketch;
use std::panic::{catch_unwind, AssertUnwindSafe};

/// A named way to build a degenerate sketch.
type Case = (&'static str, fn() -> Sketch);

/// A line from the origin to `(x, y)` whose start is held.
fn line_from_origin(sk: &mut Sketch, x: f64, y: f64) -> usize {
    let line = sk.add_line(0.0, 0.0, x, y);
    let start = sk.entities[line].param_start;
    sk.fix_param(start);
    sk.fix_param(start + 1);
    line
}

fn cases() -> Vec<Case> {
    vec![
        ("length of a zero-length line", || {
            let mut sk = Sketch::new();
            let l = line_from_origin(&mut sk, 0.0, 0.0);
            sk.add_constraint(Constraint::LineLength {
                line: l,
                length: 2.0,
            });
            sk
        }),
        ("point on a zero-length line", || {
            let mut sk = Sketch::new();
            let l = line_from_origin(&mut sk, 0.0, 0.0);
            let p = sk.add_point(1.0, 1.0);
            sk.add_constraint(Constraint::PointOnLine {
                point: PointRef::Point(p),
                line: l,
            });
            sk
        }),
        ("circle tangent to a zero-length line", || {
            let mut sk = Sketch::new();
            let l = line_from_origin(&mut sk, 0.0, 0.0);
            let c = sk.add_circle(1.0, 1.0, 0.5);
            sk.add_constraint(Constraint::Tangent { circle: c, line: l });
            sk
        }),
        ("parallel to a zero-length line", || {
            let mut sk = Sketch::new();
            let l = line_from_origin(&mut sk, 0.0, 0.0);
            let m = sk.add_line(1.0, 0.0, 2.0, 1.0);
            sk.add_constraint(Constraint::Parallel(m, l));
            sk
        }),
        ("angle with a zero-length line", || {
            let mut sk = Sketch::new();
            let l = line_from_origin(&mut sk, 0.0, 0.0);
            let m = sk.add_line(1.0, 0.0, 2.0, 1.0);
            sk.add_constraint(Constraint::Angle {
                l1: m,
                l2: l,
                angle_rad: 0.5,
            });
            sk
        }),
        ("equal length with a zero-length line", || {
            let mut sk = Sketch::new();
            let l = line_from_origin(&mut sk, 0.0, 0.0);
            let m = sk.add_line(1.0, 0.0, 2.0, 1.0);
            sk.fix_entity(m);
            sk.add_constraint(Constraint::EqualLength(m, l));
            sk
        }),
        ("symmetric about a zero-length line", || {
            let mut sk = Sketch::new();
            let l = line_from_origin(&mut sk, 0.0, 0.0);
            let a = sk.add_point(1.0, 0.0);
            let b = sk.add_point(-1.0, 0.5);
            sk.add_constraint(Constraint::Symmetric {
                p1: PointRef::Point(a),
                p2: PointRef::Point(b),
                line: l,
            });
            sk
        }),
        ("distance between coincident points", || {
            let mut sk = Sketch::new();
            let a = sk.add_point(0.0, 0.0);
            sk.fix_entity(a);
            let b = sk.add_point(0.0, 0.0);
            sk.add_constraint(Constraint::point_distance(
                PointRef::Point(a),
                PointRef::Point(b),
                1.0,
            ));
            sk
        }),
        ("point on a circle, at its centre", || {
            let mut sk = Sketch::new();
            let c = sk.add_circle(0.0, 0.0, 1.0);
            sk.fix_entity(c);
            let p = sk.add_point(0.0, 0.0);
            sk.add_constraint(Constraint::PointOnCircle {
                point: PointRef::Point(p),
                circle: c,
            });
            sk
        }),
        ("tangent circles, concentric", || {
            let mut sk = Sketch::new();
            let a = sk.add_circle(0.0, 0.0, 1.0);
            sk.fix_entity(a);
            let b = sk.add_circle(0.0, 0.0, 0.5);
            sk.add_constraint(Constraint::CircleCircleTangent {
                c1: a,
                c2: b,
                external: true,
            });
            sk
        }),
        ("point on a circle of zero radius", || {
            let mut sk = Sketch::new();
            let c = sk.add_circle(0.0, 0.0, 0.0);
            let p = sk.add_point(1.0, 0.0);
            sk.fix_entity(p);
            sk.add_constraint(Constraint::PointOnCircle {
                point: PointRef::Point(p),
                circle: c,
            });
            sk.add_constraint(Constraint::fix_point(PointRef::Center(c), 0.0, 0.0));
            sk
        }),
        ("arc of zero radius, its end held", || {
            let mut sk = Sketch::new();
            let a = sk.add_arc(0.0, 0.0, 0.0, 0.0, 1.0);
            sk.add_constraint(Constraint::fix_point(PointRef::Center(a), 0.0, 0.0));
            sk.add_constraint(Constraint::fix_point(PointRef::ArcEnd(a), 0.0, 2.0));
            sk
        }),
        ("a coordinate that is NaN", || {
            let mut sk = Sketch::new();
            let p = sk.add_point(f64::NAN, 0.0);
            sk.add_constraint(Constraint::fix_point(PointRef::Point(p), 1.0, 1.0));
            sk
        }),
        ("a coordinate that is infinite", || {
            let mut sk = Sketch::new();
            let p = sk.add_point(f64::INFINITY, 0.0);
            sk.add_constraint(Constraint::fix_point(PointRef::Point(p), 1.0, 1.0));
            sk
        }),
        ("a dimension that is NaN", || {
            let mut sk = Sketch::new();
            let l = line_from_origin(&mut sk, 1.0, 0.0);
            sk.add_constraint(Constraint::LineLength {
                line: l,
                length: f64::NAN,
            });
            sk
        }),
        ("a held coordinate that is NaN", || {
            let mut sk = Sketch::new();
            let a = sk.add_point(f64::NAN, 0.0);
            sk.fix_entity(a);
            let b = sk.add_point(1.0, 0.0);
            sk.add_constraint(Constraint::point_distance(
                PointRef::Point(a),
                PointRef::Point(b),
                1.0,
            ));
            sk
        }),
    ]
}

fn finite(sk: &Sketch) -> bool {
    sk.params_snapshot().iter().all(|v| v.is_finite())
}

#[test]
#[ignore]
fn report_degenerate() {
    for (name, build) in cases() {
        let sk = build();
        let start = sk.params_snapshot();
        let started_finite = finite(&sk);
        let solved = catch_unwind(AssertUnwindSafe(|| {
            let mut sk = build();
            let r = sk.solve();
            let after = sk.params_snapshot();
            (r, finite(&sk), after)
        }));
        let diagnosed = catch_unwind(AssertUnwindSafe(|| build().diagnose()));
        match solved {
            Ok((r, finite_after, after)) => println!(
                "{name}: {:?} converged {} residual {:e} it {} | params finite {} -> {} | moved {}",
                r.status,
                r.converged,
                r.residual,
                r.iterations,
                started_finite,
                finite_after,
                after
                    .iter()
                    .zip(&start)
                    .any(|(a, b)| a.to_bits() != b.to_bits())
            ),
            Err(_) => println!("{name}: solve PANICKED"),
        }
        match diagnosed {
            Ok(d) => println!(
                "    diagnose: dof {} rank {} conflicting {:?} redundant {:?} partially {:?}",
                d.dof, d.rank, d.conflicting, d.redundant, d.partially_redundant
            ),
            Err(_) => println!("    diagnose PANICKED"),
        }
    }
}

// ── gates ───────────────────────────────────────────────────────────────────

/// No case panics, and none leaves a parameter that is not a number: a
/// solve that cannot proceed says why and puts the geometry back.
#[test]
fn no_degenerate_case_panics_or_leaves_nan() {
    for (name, build) in cases() {
        let mut sk = build();
        let start = sk.params_snapshot();
        let r = sk.solve();
        if !r.converged {
            assert_eq!(
                sk.params_snapshot()
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                start.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                "{name}: a failed solve moved the geometry"
            );
        } else {
            assert!(finite(&sk), "{name}: converged to a non-finite state");
        }
        let _ = sk.diagnose();
    }
}

/// A line of zero length has no direction, so its length has no gradient
/// there: the solve cannot start, and says so as it does for a distance
/// between coincident points. It used to divide zero by zero and report
/// `Diverged`, which is for numbers that are not numbers.
#[test]
fn a_zero_length_line_stalls_like_coincident_points() {
    let mut line = cases()[0].1();
    let mut points = cases()
        .into_iter()
        .find(|(name, _)| *name == "distance between coincident points")
        .unwrap()
        .1();
    let (a, b) = (line.solve(), points.solve());
    assert_eq!(b.status, rcad_gcs::SolveStatus::Stalled);
    assert_eq!(a.status, b.status, "{a:?}");
}

/// The other two families on the same footing: a coordinate that is not a
/// number, and a 3D line of zero length given a length, end in a status
/// and leave the geometry as it was.
#[test]
fn the_3d_and_rigid_families_survive_degenerate_input() {
    use glam::{DQuat, DVec3};
    use rcad_gcs::rigid3d::{Constraint as Mate, ConstraintSolver, ConstraintSystem};
    use rcad_gcs::space3d::constraint::SpaceConstraint;
    use rcad_gcs::space3d::SpaceSketch;
    use rcad_gcs::SolveStatus;

    // 3D: a NaN coordinate, held to a point.
    let mut sk = SpaceSketch::new();
    let p = sk.add_point(f64::NAN, 0.0, 0.0);
    sk.add_constraint(SpaceConstraint::Fixed {
        point: rcad_gcs::space3d::entity::SpacePointRef::Point(p),
        x: 1.0,
        y: 1.0,
        z: 1.0,
    });
    let r = sk.solve();
    assert_eq!(r.status, SolveStatus::Diverged, "{r:?}");
    assert!(sk.params_snapshot()[0].is_nan());
    let _ = sk.diagnose();

    // 3D: a line of zero length given a length.
    let mut sk = SpaceSketch::new();
    let l = sk.add_line(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for k in 0..3 {
        sk.fix_param(sk.entities[l].param_start + k);
    }
    sk.add_constraint(SpaceConstraint::LineLength {
        line: l,
        length: 2.0,
    });
    let before = sk.params_snapshot();
    let r = sk.solve();
    assert!(!r.converged, "{r:?}");
    assert_ne!(r.status, SolveStatus::Diverged, "{r:?}");
    assert_eq!(sk.params_snapshot(), before);
    let _ = sk.diagnose();

    // Rigid: a part whose start position is NaN.
    let mut system = ConstraintSystem::with_components(2);
    system.add_constraint(Mate::fixed_with_orientation(
        0,
        0,
        DVec3::ZERO,
        DQuat::IDENTITY,
    ));
    system.add_constraint(Mate::coincident(1, 0, DVec3::X, 1, DVec3::ZERO));
    let result = ConstraintSolver::new().solve_from(
        &system,
        vec![DVec3::ZERO, DVec3::new(f64::NAN, 0.0, 0.0)],
        vec![DQuat::IDENTITY, DQuat::IDENTITY],
    );
    assert!(!result.success);
    assert_eq!(result.status, SolveStatus::Diverged, "{:?}", result.status);
}
