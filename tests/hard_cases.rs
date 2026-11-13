//! Hard cases for the step control: systems that are well posed but badly
//! conditioned, where a plain damped Gauss–Newton step struggles.
//!
//! `report_serial_truss` and `report_needle_triangle` are measurements,
//! not gates: they print how the solver does on each case and are ignored
//! by default. Run it optimised:
//!
//! ```text
//! cargo test --release --test hard_cases -- --ignored --nocapture
//! ```

use glam::DVec3;
use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::SpaceSketch;
use rcad_gcs::SolveOptions;

/// A serial 3D truss: every point is tied by distances to the `ties.len()`
/// points before it (three make it rigid), the first few points are fixed.
/// Errors made early amplify along the chain, which is what makes it hard.
/// Authored at the solution.
fn truss(points: usize, ties: &[usize]) -> SpaceSketch {
    // A helix keeps consecutive points well apart and never coplanar.
    let position = |i: usize| {
        let t = i as f64 * 0.9;
        DVec3::new(2.0 * t.cos(), 2.0 * t.sin(), 0.7 * i as f64)
    };
    let mut sk = SpaceSketch::new();
    let ids: Vec<usize> = (0..points)
        .map(|i| {
            let p = position(i);
            sk.add_point(p.x, p.y, p.z)
        })
        .collect();
    let anchored = ties.iter().copied().max().unwrap_or(0);
    for &id in &ids[..anchored] {
        sk.fix_entity(id);
    }
    for i in anchored..points {
        for &back in ties {
            let distance = (position(i) - position(i - back)).length();
            sk.add_constraint(SpaceConstraint::point_distance(
                ids[i],
                ids[i - back],
                distance,
            ));
        }
    }
    sk
}

/// Largest distance between a point of the solved truss and where it was
/// authored: zero when the solve came back to the same configuration.
fn displacement(solved: &SpaceSketch, authored: &SpaceSketch) -> f64 {
    solved
        .params_snapshot()
        .iter()
        .zip(authored.params_snapshot())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max)
}

struct Outcome {
    converged: bool,
    iterations: usize,
    residual: f64,
    displacement: f64,
    millis: f64,
}

fn run(points: usize, ties: &[usize], amplitude: f64, max_iterations: usize) -> Outcome {
    let authored = truss(points, ties);
    let mut sk = authored.clone();
    sk.drift_free_params(amplitude, 7);
    let options = SolveOptions {
        max_iterations,
        ..SolveOptions::default()
    };
    let started = std::time::Instant::now();
    let result = sk.solve_with(&options);
    Outcome {
        converged: result.converged,
        iterations: result.iterations,
        residual: result.residual,
        displacement: displacement(&sk, &authored),
        millis: started.elapsed().as_secs_f64() * 1e3,
    }
}

#[test]
#[ignore]
fn report_serial_truss() {
    println!("serial truss (ties 1,2,3), default options, iteration limit 100");
    println!("points  amp   converged  iterations  residual    moved       ms");
    for amplitude in [0.05, 0.2] {
        for points in [20usize, 80, 120, 160, 200, 260] {
            let o = run(points, &[1, 2, 3], amplitude, 100);
            println!(
                "{points:>6}  {amplitude:<4}  {:<9}  {:>10}  {:<10.3e}  {:<10.3e}  {:.0}",
                o.converged, o.iterations, o.residual, o.displacement, o.millis
            );
        }
    }
    println!("braced truss (ties 1,2,3,5), for comparison");
    for points in [20usize, 80, 260] {
        let o = run(points, &[1, 2, 3, 5], 0.05, 100);
        println!(
            "{points:>6}  0.05  {:<9}  {:>10}  {:<10.3e}  {:<10.3e}  {:.0}",
            o.converged, o.iterations, o.residual, o.displacement, o.millis
        );
    }
}

/// From a small perturbation the serial truss comes back to where it was
/// authored, at every length tried here. Before the step control became a
/// trust region the 160-point truss did not converge within 100 iterations
/// (and neither did the 260-point one, which the report above covers; it is
/// too slow for an unoptimised test run).
#[test]
fn serial_trusses_converge_back() {
    for points in [20usize, 80, 160] {
        let o = run(points, &[1, 2, 3], 0.05, 100);
        assert!(o.converged, "{points} points: residual {:.3e}", o.residual);
        assert!(
            o.displacement < 1e-6,
            "{points} points moved {:.3e}",
            o.displacement
        );
    }
}

/// A needle triangle: `A` and `B` held 1 apart, `C` 2 from `A` and `1 + δ`
/// from `B`. `C` then sits `≈ 2√δ` off the line `AB`, the two distance
/// rows are nearly parallel (`κ(J) ≈ 1/(2√δ)`, squared in `JᵀJ`), and at
/// `δ = 0` the circles touch: a near-double root. Prints, per `δ`, start
/// and tolerance, what the solve does and how far `C` lands from the exact
/// point, relative to its height off the line.
///
/// Measured (P2.2, see `docs/CATCHUP_PLAN.md` §10): at `δ ≥ 1e-8` all is
/// well. At the default tolerance every `δ ≤ 1e-10` reports `Converged`
/// with the height wrong by up to its own size or more (the residual does
/// not bound the position near tangency). With a tight tolerance `δ = 1e-12`
/// and `1e-14` end as `MaxIterations`, while Newton with an exact 2 × 2
/// solve converges in 12–21 iterations; with the regularisation `λ`
/// lowered from 1e-10 to 1e-20 the solver does too at `1e-12`, and at
/// `1e-14` one start still stalls.
#[test]
#[ignore]
fn report_needle_triangle() {
    use rcad_gcs::constraint::Constraint;
    use rcad_gcs::entity::PointRef;
    use rcad_gcs::Sketch;

    for delta in [1e-6, 1e-8, 1e-10, 1e-12, 1e-14] {
        let d: f64 = delta;
        let x = 2.0 - d - 0.5 * d * d;
        let y = (4.0 - x * x).sqrt();
        for (off, tolerance) in [(1e-3, 1e-10), (1e-3, 1e-16), (0.1, 1e-16)] {
            let mut sk = Sketch::new();
            let a = sk.add_point(0.0, 0.0);
            let b = sk.add_point(1.0, 0.0);
            sk.fix_entity(a);
            sk.fix_entity(b);
            let c = sk.add_point(x + off, y + off);
            sk.add_constraint(Constraint::point_distance(
                PointRef::Point(a),
                PointRef::Point(c),
                2.0,
            ));
            sk.add_constraint(Constraint::point_distance(
                PointRef::Point(b),
                PointRef::Point(c),
                1.0 + d,
            ));
            let r = sk.solve_with(&SolveOptions {
                tolerance,
                ..Default::default()
            });
            let p = sk.point_coords(PointRef::Point(c));
            println!(
                "δ {d:e} start {off:e} tolerance {tolerance:e}: {:?} it {} residual {:.1e} | C off by {:.1e} of its height",
                r.status,
                r.iterations,
                r.residual,
                ((p.y - y) / y).abs().max(((p.x - x) / y).abs())
            );
        }
    }
}
