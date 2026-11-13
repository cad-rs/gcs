//! The same shape at different sizes and places.
//!
//! A rectangle with a circle on one corner, fully dimensioned, built at a
//! size `s` and an offset `o` from the origin. Geometrically it is one
//! model; only the numbers differ. Whether it solves should not depend on
//! them.
//!
//! `report_scale` prints what the solver does on each; it is a measurement
//! and ignored by default:
//!
//! ```text
//! cargo test --release --test scale -- --ignored --nocapture
//! ```

use rcad_gcs::constraint::Constraint;
use rcad_gcs::entity::PointRef;
use rcad_gcs::{Sketch, SolveResult, SolveStatus};

/// A fully constrained `4s × 3s` rectangle whose lower-left corner is pinned
/// at `(o, o)`, with a circle of radius `s` centred on its upper-right
/// corner and a point held on that circle on the diagonal. Every entity
/// starts a relative `off` away from where the constraints want it.
fn model(s: f64, o: f64, off: f64) -> Sketch {
    let (w, h) = (4.0 * s, 3.0 * s);
    let d = off * s;
    let mut sk = Sketch::new();
    let bottom = sk.add_line(o + d, o - d, o + w - d, o + d);
    let right = sk.add_line(o + w + d, o + d, o + w - d, o + h + d);
    let top = sk.add_line(o + w + d, o + h - d, o - d, o + h + d);
    let left = sk.add_line(o + d, o + h + d, o - d, o + d);
    let circle = sk.add_circle(o + w - d, o + h + d, s * (1.0 + off));
    let point = sk.add_point(o + w + 0.8 * s + d, o + h + 0.6 * s - d);

    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(bottom), o, o));
    for (end, start) in [(bottom, right), (right, top), (top, left), (left, bottom)] {
        sk.add_constraint(Constraint::coincident(
            PointRef::LineEnd(end),
            PointRef::LineStart(start),
        ));
    }
    sk.add_constraint(Constraint::Horizontal(bottom));
    sk.add_constraint(Constraint::Vertical(right));
    sk.add_constraint(Constraint::Horizontal(top));
    sk.add_constraint(Constraint::Vertical(left));
    sk.add_constraint(Constraint::LineLength {
        line: bottom,
        length: w,
    });
    sk.add_constraint(Constraint::LineLength {
        line: right,
        length: h,
    });
    sk.add_constraint(Constraint::coincident(
        PointRef::Center(circle),
        PointRef::LineEnd(right),
    ));
    sk.add_constraint(Constraint::Radius { circle, radius: s });
    sk.add_constraint(Constraint::PointOnCircle {
        point: PointRef::Point(point),
        circle,
    });
    // The point sits on the 3-4-5 diagonal through the centre.
    sk.add_constraint(Constraint::distance_x(
        PointRef::Point(point),
        PointRef::Center(circle),
        0.8 * s,
    ));
    sk
}

/// The largest error, relative to the size of the model, of the solved
/// corner positions against where they belong.
fn worst_relative_error(sk: &Sketch, s: f64, o: f64) -> f64 {
    let expect = [
        (o, o),
        (o + 4.0 * s, o),
        (o + 4.0 * s, o + 3.0 * s),
        (o, o + 3.0 * s),
    ];
    (0..4)
        .map(|line| {
            let p = sk.point_coords(PointRef::LineStart(line));
            let (x, y) = expect[line];
            (p.x - x).abs().max((p.y - y).abs()) / s
        })
        .fold(0.0_f64, f64::max)
}

fn solve(s: f64, o: f64) -> (SolveResult, f64, i64) {
    let mut sk = model(s, o, 0.02);
    let dof = sk.dof();
    let result = sk.solve();
    (result, worst_relative_error(&sk, s, o), dof)
}

#[test]
#[ignore]
fn report_scale() {
    println!("size      offset    converged  status         iters  residual    rel.error   dof");
    for s in [1.0e-6, 1.0e-3, 1.0, 1.0e3, 1.0e6] {
        for o in [0.0, 1.0e3 * s, 1.0e6 * s] {
            let (r, error, dof) = solve(s, o);
            println!(
                "{s:<8.0e}  {o:<8.0e}  {:<9}  {:<13}  {:>5}  {:<10.3e}  {:<10.3e}  {dof}",
                r.converged,
                format!("{:?}", r.status),
                r.iterations,
                r.residual,
                error
            );
        }
    }
}

/// A 3-4-5 triangle turned by an angle that is no round number, hung on a
/// corner at `(o, 0.7·o)`: nothing in it lands on values a binary float
/// holds exactly. `ab` is its first side, held at 0.3 rad to a fixed
/// reference line; the third corner is placed by its two distances.
fn triangle(s: f64, o: f64, off: f64) -> (Sketch, [f64; 2]) {
    let (ax, ay) = (o, 0.7 * o);
    let turn: f64 = 0.3;
    let b = (ax + 3.0 * s * turn.cos(), ay + 3.0 * s * turn.sin());
    // C is at right angles to AB from B, 4s away.
    let c = (b.0 - 4.0 * s * turn.sin(), b.1 + 4.0 * s * turn.cos());
    let d = off * s;

    let mut sk = Sketch::new();
    let reference = sk.add_line(ax, ay, ax + s, ay);
    sk.fix_entity(reference);
    let ab = sk.add_line(ax + d, ay - d, b.0 - d, b.1 + d);
    let corner = sk.add_point(c.0 + d, c.1 + d);
    sk.add_constraint(Constraint::coincident(
        PointRef::LineStart(ab),
        PointRef::LineStart(reference),
    ));
    sk.add_constraint(Constraint::LineLength {
        line: ab,
        length: 3.0 * s,
    });
    sk.add_constraint(Constraint::Angle {
        l1: reference,
        l2: ab,
        angle_rad: turn,
    });
    sk.add_constraint(Constraint::point_distance(
        PointRef::Point(corner),
        PointRef::LineEnd(ab),
        4.0 * s,
    ));
    sk.add_constraint(Constraint::point_distance(
        PointRef::Point(corner),
        PointRef::LineStart(ab),
        5.0 * s,
    ));
    (sk, [c.0, c.1])
}

#[test]
#[ignore]
fn report_scale_triangle() {
    println!("size      offset    converged  status         iters  residual    rel.error   dof");
    for s in [1.0e-3, 1.0, 1.0e3] {
        for factor in [0.0, 1.0e3, 1.0e6, 1.0e9, 1.0e12] {
            let o = factor * s * std::f64::consts::FRAC_PI_3;
            let (mut sk, c) = triangle(s, o, 0.02);
            let dof = sk.dof();
            let r = sk.solve();
            let p = sk.point_coords(PointRef::Point(2));
            let error = (p.x - c[0]).abs().max((p.y - c[1]).abs()) / s;
            println!(
                "{s:<8.0e}  {o:<8.1e}  {:<9}  {:<13}  {:>5}  {:<10.3e}  {:<10.3e}  {dof}",
                r.converged,
                format!("{:?}", r.status),
                r.iterations,
                r.residual,
                error
            );
        }
    }
}

// ── gates ───────────────────────────────────────────────────────────────────

/// What the absolute tolerance leaves at ordinary sizes, plus what the
/// coordinates themselves resolve at magnitude `x`.
fn error_allowed(x: f64) -> f64 {
    1.0e-9 + 64.0 * f64::EPSILON * x
}

/// The same two models, from a thousandth of a unit to a million units
/// across and from the origin to a billion times their own size away:
/// every one solves, to the accuracy its coordinates can hold. (Smaller
/// and larger at the origin: `the_triangle_solves_exactly_at_any_size`.)
#[test]
fn the_same_model_solves_at_any_size_and_place() {
    for s in [1.0e-3, 1.0, 1.0e3, 1.0e6] {
        for factor in [0.0, 1.0e3, 1.0e6, 1.0e9] {
            // The rectangle: round numbers.
            let o = factor * s;
            let (r, relative, dof) = solve(s, o);
            assert_eq!(dof, 0);
            assert!(r.converged, "rectangle {s:e} at {o:e}: {r:?}");
            assert!(
                relative * s <= error_allowed(o + 5.0 * s),
                "rectangle {s:e} at {o:e}: off by {:e}",
                relative * s
            );

            // The triangle: nothing representable exactly.
            let o = factor * s * std::f64::consts::FRAC_PI_3;
            let (mut sk, c) = triangle(s, o, 0.02);
            let r = sk.solve();
            assert!(r.converged, "triangle {s:e} at {o:e}: {r:?}");
            let p = sk.point_coords(PointRef::Point(2));
            let error = (p.x - c[0]).abs().max((p.y - c[1]).abs());
            assert!(
                error <= error_allowed(2.0 * o + 5.0 * s),
                "triangle {s:e} at {o:e}: off by {error:e}"
            );
            // Near the origin the tolerance is met outright; the precision
            // exit is for where it cannot be.
            if factor <= 1.0e3 && s <= 1.0 {
                assert_eq!(r.status, SolveStatus::Converged, "{s:e} at {o:e}");
            }
            // A billion units out, a float is spaced 1e-7 apart: the
            // tolerance cannot be met there, and the solve says so.
            if o >= 1.0e9 {
                assert_eq!(
                    r.status,
                    SolveStatus::ConvergedToPrecision,
                    "{s:e} at {o:e}"
                );
                assert!(r.residual > rcad_gcs::numerics::RESIDUAL_TOL);
            }
        }
    }
}

/// Being far from the origin excuses rounding, not a contradiction: a
/// second, different length for the same side is still a failed solve, and
/// the geometry is left where it was. (Up to a point: at 1e12 a float is
/// spaced 1e-4 apart, and a contradiction of 0.01 is a few dozen of those.
/// There it reads as rounding, and no arithmetic in doubles could tell.)
#[test]
fn a_conflict_far_from_the_origin_is_still_a_conflict() {
    for factor in [0.0, 1.0e6, 1.0e9] {
        let o = factor * std::f64::consts::FRAC_PI_3;
        let (mut sk, _) = triangle(1.0, o, 0.02);
        // The corner is already 5 from the first end of `ab`.
        sk.add_constraint(Constraint::point_distance(
            PointRef::Point(2),
            PointRef::LineStart(1),
            5.01,
        ));
        let before = sk.params_snapshot();
        let r = sk.solve();
        assert!(!r.converged, "at {o:e}: {r:?}");
        assert_eq!(r.status, SolveStatus::Stalled, "at {o:e}");
        assert_eq!(sk.params_snapshot(), before);
        assert!(!sk.diagnose().conflicting.is_empty(), "at {o:e}");
    }
}

/// Asking for more than double precision holds is not a failure either: a
/// tolerance of 1e-18 at unit scale ends at the resolution of the numbers.
#[test]
fn an_unreachable_tolerance_ends_at_the_precision_of_the_numbers() {
    let (mut sk, c) = triangle(1.0, 10.0, 0.02);
    let r = sk.solve_with(&rcad_gcs::SolveOptions {
        tolerance: 1.0e-18,
        ..Default::default()
    });
    assert!(r.converged, "{r:?}");
    assert_eq!(r.status, SolveStatus::ConvergedToPrecision);
    let p = sk.point_coords(PointRef::Point(2));
    assert!((p.x - c[0]).abs().max((p.y - c[1]).abs()) < 1e-13);
}

/// The 3D-geometric family: a tetrahedron of distances on a fixed corner,
/// far from the origin.
#[test]
fn a_3d_sketch_solves_far_from_the_origin() {
    use rcad_gcs::space3d::constraint::SpaceConstraint;
    use rcad_gcs::space3d::SpaceSketch;

    for o in [0.0, 1.0e6, 1.0e9, 1.0e12] {
        let o = o * std::f64::consts::FRAC_PI_3;
        let corners = [
            (o, 0.5 * o, -0.3 * o),
            (o + 3.0, 0.5 * o, -0.3 * o),
            (o, 0.5 * o + 4.0, -0.3 * o),
            (o, 0.5 * o, -0.3 * o + 2.0),
        ];
        let mut sk = SpaceSketch::new();
        let ids: Vec<usize> = corners
            .iter()
            .enumerate()
            .map(|(i, c)| {
                // Start a little off, except the corner that is pinned.
                let d = if i == 0 { 0.0 } else { 0.03 };
                sk.add_point(c.0 + d, c.1 - d, c.2 + d)
            })
            .collect();
        sk.fix_entity(ids[0]);
        // Hold the frame: the second corner along X, the third in the XY
        // plane, so that the six distances leave nothing free.
        sk.fix_param(sk.entities[ids[1]].param_start + 1);
        sk.fix_param(sk.entities[ids[1]].param_start + 2);
        sk.fix_param(sk.entities[ids[2]].param_start + 2);
        let mut params = sk.params_snapshot();
        params[sk.entities[ids[1]].param_start + 1] = corners[1].1;
        params[sk.entities[ids[1]].param_start + 2] = corners[1].2;
        params[sk.entities[ids[2]].param_start + 2] = corners[2].2;
        sk.set_params(params);
        for a in 0..4 {
            for b in a + 1..4 {
                let (p, q) = (corners[a], corners[b]);
                let distance =
                    ((p.0 - q.0).powi(2) + (p.1 - q.1).powi(2) + (p.2 - q.2).powi(2)).sqrt();
                sk.add_constraint(SpaceConstraint::point_distance(ids[a], ids[b], distance));
            }
        }
        assert_eq!(sk.dof(), 0, "at {o:e}");
        let r = sk.solve();
        assert!(r.converged, "at {o:e}: {r:?}");
        for (id, c) in ids.iter().zip(corners) {
            let p = sk.point_coords(*id).unwrap();
            let error = (p.0 - c.0)
                .abs()
                .max((p.1 - c.1).abs())
                .max((p.2 - c.2).abs());
            assert!(
                error <= error_allowed(2.0 * o + 10.0),
                "at {o:e}: off by {error:e}"
            );
        }
    }
}

/// The rigid-pose family: a part locked to a ground that stands far from
/// the origin, asked for with a tolerance the translation cannot resolve
/// there.
#[test]
fn an_assembly_solves_far_from_the_origin() {
    use glam::{DQuat, DVec3};
    use rcad_gcs::rigid3d::{
        Constraint as Mate, ConstraintSolver, ConstraintSystem, PoseSolveOptions,
    };

    let solver = ConstraintSolver::with_options(PoseSolveOptions {
        tolerance: 1.0e-12,
        ..PoseSolveOptions::default()
    });
    let turn = DQuat::from_euler(glam::EulerRot::XYZ, 0.4, -0.3, 0.7);
    let relative = DVec3::new(1.5, -2.0, 0.5);
    for o in [0.0, 1.0e6, 1.0e9, 1.0e12] {
        let ground = DVec3::new(o, -0.6 * o, 0.25 * o) * std::f64::consts::FRAC_PI_3;
        let ground_turn = DQuat::from_rotation_z(0.3);
        let mut system = ConstraintSystem::with_components(3);
        system.add_constraint(Mate::fixed_with_orientation(0, 0, ground, ground_turn));
        system.add_constraint(Mate::lock(1, 0, 1, relative, turn));
        // A second part hung on the first by three points.
        for (k, local) in [DVec3::ZERO, DVec3::X * 2.0, DVec3::new(0.0, 3.0, 1.0)]
            .into_iter()
            .enumerate()
        {
            system.add_constraint(Mate::coincident(
                10 + k,
                1,
                local + DVec3::new(0.0, 0.0, 4.0),
                2,
                local,
            ));
        }
        // Start the parts near the ground, not at the origin.
        let start = vec![
            ground,
            ground + DVec3::splat(0.3),
            ground - DVec3::splat(0.2),
        ];
        let result = solver.solve_from(
            &system,
            start,
            vec![ground_turn, DQuat::IDENTITY, DQuat::from_rotation_x(0.2)],
        );
        assert!(
            result.success,
            "at {o:e}: {:?} error {:e}",
            result.status, result.error
        );
        let expected = ground + ground_turn * relative;
        let error = (result.positions[1] - expected).length();
        assert!(
            error <= error_allowed(4.0 * ground.length() + 10.0),
            "at {o:e}: off by {error:e}"
        );
        let hung = result.positions[1] + result.rotations[1] * DVec3::new(0.0, 0.0, 4.0);
        assert!(
            (result.positions[2] - hung).length() <= error_allowed(4.0 * ground.length() + 10.0),
            "at {o:e}"
        );
    }
}

/// The triangle at the origin, from a millionth of a unit to a million
/// units across, at the default tolerance and at a tighter one: its third
/// corner lands within a rounding-level fraction of its size every time.
///
/// Every 2D equation is a length, the angle included (it is the sine of the
/// angle error times the lengths of its lines). Before that, the angle was
/// stated as a product of two lengths: at 1e-5 across it fell under the
/// tolerance with the corner 0.12% of the size out of place, and at 1e-6
/// with a tolerance of 1e-16 the solve never converged.
#[test]
fn the_triangle_solves_exactly_at_any_size() {
    for tolerance in [1.0e-10, 1.0e-16] {
        for s in [1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0, 1.0e3, 1.0e6] {
            let (mut sk, c) = triangle(s, 0.0, 0.02);
            let r = sk.solve_with(&rcad_gcs::SolveOptions {
                tolerance,
                ..Default::default()
            });
            assert!(r.converged, "size {s:e}, tolerance {tolerance:e}: {r:?}");
            let p = sk.point_coords(PointRef::Point(2));
            let error = (p.x - c[0]).abs().max((p.y - c[1]).abs()) / s;
            assert!(
                error <= 1.0e-12,
                "size {s:e}, tolerance {tolerance:e}: off by {error:e} of the size"
            );
        }
    }
}

/// The triangle at the origin, small and large, at the default tolerance
/// and at tighter ones (the measurement behind the gate above).
#[test]
#[ignore]
fn report_small_geometry() {
    for tolerance in [1.0e-10, 1.0e-16, 1.0e-22] {
        for s in [1.0e6, 1.0e3, 1.0, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6] {
            let (mut sk, c) = triangle(s, 0.0, 0.02);
            let r = sk.solve_with(&rcad_gcs::SolveOptions {
                tolerance,
                ..Default::default()
            });
            let p = sk.point_coords(PointRef::Point(2));
            let error = (p.x - c[0]).abs().max((p.y - c[1]).abs()) / s;
            println!(
                "tol {tolerance:e} size {s:e}: {:?} it {} residual {:e} relative error {error:e}",
                r.status, r.iterations, r.residual
            );
        }
    }
}

/// The 3D-geometric family's equal length, on small lines: a held line of
/// length `3s` and a second one that may only stretch along X. Its residual
/// was the difference of the squared lengths, an area: at 1e-8 across it
/// starts under the default tolerance, and the solve reported converged
/// with the line still 10% short.
#[test]
fn a_3d_equal_length_holds_on_small_lines() {
    use rcad_gcs::space3d::constraint::SpaceConstraint;
    use rcad_gcs::space3d::SpaceSketch;

    for (s, tolerance) in [1.0e-8, 1.0e-7, 1.0e-6, 1.0e-3, 1.0, 1.0e3]
        .into_iter()
        .flat_map(|s| [(s, 1.0e-10), (s, 1.0e-16)])
    {
        let mut sk = SpaceSketch::new();
        let held = sk.add_line(0.0, 0.0, 0.0, 0.0, 3.0 * s, 0.0);
        sk.fix_entity(held);
        let line = sk.add_line(0.0, 0.0, s, 2.7 * s, 0.0, s);
        let start = sk.entities[line].param_start;
        for k in [0, 1, 2, 4, 5] {
            sk.fix_param(start + k);
        }
        sk.add_constraint(SpaceConstraint::EqualLineLength(held, line));
        let r = sk.solve_with(&rcad_gcs::SolveOptions {
            tolerance,
            ..Default::default()
        });
        assert!(r.converged, "size {s:e}, tolerance {tolerance:e}: {r:?}");
        let x = sk.params_snapshot()[start + 3];
        assert!(
            (x - 3.0 * s).abs() <= 1.0e-12 * s,
            "size {s:e}, tolerance {tolerance:e}: the end is at {x:e}, not {:e}",
            3.0 * s
        );
    }
}

/// The triangle without its angle, so free to turn about its pinned corner:
/// how far its three lengths are from holding after a solve, relative to
/// its size. A component with freedom left gets no convergence polish, so
/// this is what the tolerance alone buys.
#[test]
#[ignore]
fn report_small_geometry_free() {
    for s in [1.0, 1.0e-2, 1.0e-4, 1.0e-6] {
        let (mut sk, _) = triangle(s, 0.0, 0.02);
        assert!(sk.remove_constraint(2));
        let r = sk.solve();
        let a = sk.point_coords(PointRef::LineStart(1));
        let b = sk.point_coords(PointRef::LineEnd(1));
        let c = sk.point_coords(PointRef::Point(2));
        let d = |p: glam::DVec2, q: glam::DVec2| (p - q).length();
        let worst = [(d(a, b), 3.0), (d(b, c), 4.0), (d(a, c), 5.0)]
            .iter()
            .map(|(l, k)| (l - k * s).abs() / s)
            .fold(0.0_f64, f64::max);
        println!(
            "size {s:e}: {:?} it {} residual {:e} worst relative error {worst:e}",
            r.status, r.iterations, r.residual
        );
    }
}
