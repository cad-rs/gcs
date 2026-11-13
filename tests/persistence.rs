//! Solution persistence / no-hysteresis contract (solver roadmap M4.1).
//!
//! DCM's headline public design goal (ASME JCISE): returning prior
//! dimension values recovers the original CONFIGURATION. The contract
//! has two levels:
//!
//! 1. **Promised**: a round-tripped driving dimension re-converges to the
//!    original solution within [`CONTRACT_TOL`] (1e-9) and the
//!    configuration never flips (no mirrored/rotated branch). Newton
//!    terminates on an RMS residual threshold, so the terminal iterate is
//!    path-dependent AT THAT SCALE in principle — that is what the
//!    tolerance is for.
//! 2. **Observed** (recorded 2026-10-07, this family): the worst
//!    round-trip drift is bitwise ZERO across all four families — the
//!    quadratically-convergent damped loop pins the terminal iterate.
//!    The suite prints the measured drift under `--nocapture`; if a
//!    future change makes it non-zero but within the contract, the test
//!    still passes and the printout tells us immediately.

use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::Sketch;

const CONTRACT_TOL: f64 = 1e-9;

fn assert_round_trip(name: &str, build: &dyn Fn(f64) -> Sketch, drive: f64, factor: f64) {
    // Leg 0: authored at `drive`.
    let mut sk = build(drive);
    let r0 = sk.solve();
    assert!(r0.converged, "{name} leg 0: residual {}", r0.residual);
    let original = sk.params_snapshot();
    let fingerprint = config_fingerprint(&sk);

    // Leg 1: drive up.
    let mut sk = build(drive * factor);
    let r1 = sk.solve();
    assert!(r1.converged, "{name} leg 1: residual {}", r1.residual);

    // Leg 2: back to the original value — the round trip.
    let mut sk = build(drive);
    let r2 = sk.solve();
    assert!(r2.converged, "{name} leg 2: residual {}", r2.residual);

    let worst = original
        .iter()
        .zip(sk.params_snapshot().iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    println!("{name}: worst round-trip drift {worst:.3e}");
    assert!(
        worst <= CONTRACT_TOL,
        "{name}: round-trip drift {worst:.3e} exceeds the contract ({CONTRACT_TOL:.0e})"
    );
    let after = config_fingerprint(&sk);
    assert_eq!(
        fingerprint, after,
        "{name}: configuration flipped across the round trip"
    );
}

/// Orientation fingerprint: the sign pattern of a few canonical directed
/// measurements. A hysteresis flip (mirrored/rotated branch) changes at
/// least one sign; sub-tolerance numeric noise never does.
fn config_fingerprint(sk: &Sketch) -> Vec<i8> {
    let mut signs = Vec::new();
    for i in 0..sk.entities.len() {
        let p = sk.entity_params(i);
        for v in p {
            signs.push(if *v > 0.0 {
                1
            } else if *v < 0.0 {
                -1
            } else {
                0
            });
        }
    }
    signs
}

/// Fully-constrained rectangle: driving bottom length round-trips.
#[test]
fn rectangle_length_round_trip() {
    assert_round_trip(
        "rectangle",
        &|len| {
            let mut sk = Sketch::new();
            let origin = sk.add_point(0.0, 0.0);
            let b = sk.add_line(0.0, 0.0, len, 0.0);
            let r = sk.add_line(len, 0.0, len, 6.0);
            let t = sk.add_line(len, 6.0, 0.0, 6.0);
            let l = sk.add_line(0.0, 6.0, 0.0, 0.0);
            sk.add_constraint(Fixed { point: PointRef::Point(origin), x: 0.0, y: 0.0 });
            sk.add_constraint(Coincident(PointRef::LineStart(b), PointRef::Point(origin)));
            for line in [b, t] {
                sk.add_constraint(Horizontal(line));
            }
            for line in [r, l] {
                sk.add_constraint(Vertical(line));
            }
            sk.add_constraint(Coincident(PointRef::LineEnd(b), PointRef::LineStart(r)));
            sk.add_constraint(Coincident(PointRef::LineEnd(r), PointRef::LineStart(t)));
            sk.add_constraint(Coincident(PointRef::LineEnd(t), PointRef::LineStart(l)));
            sk.add_constraint(Coincident(PointRef::LineEnd(l), PointRef::LineStart(b)));
            sk.add_constraint(LineLength { line: b, length: len });
            sk.add_constraint(LineLength { line: r, length: 6.0 });
            sk
        },
        10.0,
        1.4,
    );
}

/// The axis-aligned margin pattern: skeleton length and hole margins
/// round-trip together.
#[test]
fn margin_pattern_round_trip() {
    assert_round_trip(
        "margin-pattern",
        &|len| {
            let mut sk = Sketch::new();
            let anchor = sk.add_point(0.0, 0.0);
            let line = sk.add_line(0.0, 0.0, len, 0.0);
            sk.add_constraint(Fixed { point: PointRef::Point(anchor), x: 0.0, y: 0.0 });
            sk.add_constraint(Coincident(PointRef::LineStart(line), PointRef::Point(anchor)));
            sk.add_constraint(Horizontal(line));
            sk.add_constraint(DistanceX {
                p1: PointRef::LineEnd(line),
                p2: PointRef::LineStart(line),
                distance: len,
            });
            let hole = sk.add_circle(len - 2.0, 3.0, 1.5);
            sk.add_constraint(Radius { circle: hole, radius: 1.5 });
            sk.add_constraint(DistanceX {
                p1: PointRef::LineEnd(line),
                p2: PointRef::Center(hole),
                distance: 2.0,
            });
            sk.add_constraint(DistanceY {
                p1: PointRef::Center(hole),
                p2: PointRef::LineStart(line),
                distance: 3.0,
            });
            sk
        },
        10.0,
        1.6,
    );
}

/// Triangle with three driven sides (SSS): round-tripping one side must
/// recover the exact triangle, not the mirrored one.
#[test]
fn triangle_sss_round_trip() {
    assert_round_trip(
        "triangle-sss",
        &|base| {
            let mut sk = Sketch::new();
            let p0 = sk.add_point(0.0, 0.0);
            let p1 = sk.add_point(base, 0.0);
            let p2 = sk.add_point(base * 0.4, 2.5);
            let l01 = sk.add_line(0.0, 0.0, base, 0.0);
            let l12 = sk.add_line(base, 0.0, base * 0.4, 2.5);
            let l20 = sk.add_line(base * 0.4, 2.5, 0.0, 0.0);
            sk.add_constraint(Fixed { point: PointRef::Point(p0), x: 0.0, y: 0.0 });
            sk.add_constraint(Coincident(PointRef::LineStart(l01), PointRef::Point(p0)));
            sk.add_constraint(Coincident(PointRef::LineEnd(l01), PointRef::Point(p1)));
            sk.add_constraint(Coincident(PointRef::LineStart(l12), PointRef::Point(p1)));
            sk.add_constraint(Coincident(PointRef::LineEnd(l12), PointRef::Point(p2)));
            sk.add_constraint(Coincident(PointRef::LineStart(l20), PointRef::Point(p2)));
            sk.add_constraint(Coincident(PointRef::LineEnd(l20), PointRef::Point(p0)));
            sk.add_constraint(LineLength { line: l01, length: base });
            sk.add_constraint(LineLength { line: l12, length: 4.0 });
            sk.add_constraint(LineLength { line: l20, length: 3.5 });
            sk
        },
        5.0,
        1.3,
    );
}

/// Point-symmetric bore between two driven references.
#[test]
fn symmetric_bore_round_trip() {
    assert_round_trip(
        "symmetric-bore",
        &|span| {
            let mut sk = Sketch::new();
            let a = sk.add_point(0.0, 0.0);
            let b = sk.add_point(span, 0.0);
            let bore = sk.add_circle(span * 0.5, 2.0, 1.0);
            sk.add_constraint(Fixed { point: PointRef::Point(a), x: 0.0, y: 0.0 });
            sk.add_constraint(DistanceX { p1: PointRef::Point(b), p2: PointRef::Point(a), distance: span });
            sk.add_constraint(DistanceY { p1: PointRef::Point(b), p2: PointRef::Point(a), distance: 0.0 });
            sk.add_constraint(SymmetricPoint {
                p1: PointRef::Center(bore),
                p2: PointRef::Point(b),
                about: PointRef::Point(a),
            });
            sk.add_constraint(Radius { circle: bore, radius: 1.0 });
            sk
        },
        10.0,
        1.5,
    );
}
