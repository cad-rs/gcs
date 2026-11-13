//! Solver benchmarks. Run with `cargo bench`; see `benches/BASELINE.md` for
//! recorded numbers.
//!
//! The 2D fixtures are rectangles (an origin point plus four lines, closed
//! corner to corner, with horizontal/vertical sides and driven width and
//! height: 18 parameters and 18 equations each):
//!
//! - **chain** — every rectangle hangs off the previous one's corner, so the
//!   whole sketch is ONE connected component. This is the scaling test of
//!   the linear algebra.
//! - **field** — independent, individually anchored rectangles. This is the
//!   scaling test of the per-component decomposition.
//!
//! The 3D-geometric fixture is a **truss** of points tied together by
//! distances, again one connected component.
//!
//! The assembly fixtures are a ground body plus rigid bodies: a **chain**
//! (each body mated to the previous one) and a **star** (each body pinned
//! to the ground), both a single connected cluster.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use glam::{DQuat, DVec3};
use rcad_gcs::constraint::Constraint;
use rcad_gcs::entity::{EntityId, PointRef};
use rcad_gcs::rigid3d::{
    diagnose_assembly, Constraint as Mate, ConstraintSolver, ConstraintSystem,
};
use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::SpaceSketch;
use rcad_gcs::Sketch;

const WIDTH: f64 = 4.0;
const HEIGHT: f64 = 3.0;

/// Where a rectangle's origin corner is held.
enum Anchor {
    /// Pinned at its authored position.
    Fixed,
    /// Coincident with a point of an earlier rectangle.
    At(PointRef),
}

/// Add one rectangle with its origin corner at `(ox, oy)`; returns the id
/// of its bottom line.
fn push_rectangle(sk: &mut Sketch, ox: f64, oy: f64, anchor: Anchor) -> EntityId {
    let origin = sk.add_point(ox, oy);
    let bottom = sk.add_line(ox, oy, ox + WIDTH, oy);
    let right = sk.add_line(ox + WIDTH, oy, ox + WIDTH, oy + HEIGHT);
    let top = sk.add_line(ox + WIDTH, oy + HEIGHT, ox, oy + HEIGHT);
    let left = sk.add_line(ox, oy + HEIGHT, ox, oy);

    match anchor {
        Anchor::Fixed => sk.add_constraint(Constraint::fix_point(origin, ox, oy)),
        Anchor::At(point) => {
            sk.add_constraint(Constraint::coincident(PointRef::Point(origin), point))
        }
    }
    sk.add_constraint(Constraint::coincident(
        PointRef::LineStart(bottom),
        PointRef::Point(origin),
    ));
    for (end, start) in [(bottom, right), (right, top), (top, left), (left, bottom)] {
        sk.add_constraint(Constraint::coincident(
            PointRef::LineEnd(end),
            PointRef::LineStart(start),
        ));
    }
    sk.add_constraint(Constraint::Horizontal(bottom));
    sk.add_constraint(Constraint::Horizontal(top));
    sk.add_constraint(Constraint::Vertical(right));
    sk.add_constraint(Constraint::Vertical(left));
    sk.add_constraint(Constraint::LineLength {
        line: bottom,
        length: WIDTH,
    });
    sk.add_constraint(Constraint::LineLength {
        line: right,
        length: HEIGHT,
    });
    bottom
}

/// `count` rectangles in one connected component, authored at the solution.
fn rectangle_chain(count: usize) -> Sketch {
    let mut sk = Sketch::new();
    let mut previous: Option<EntityId> = None;
    for i in 0..count {
        let anchor = match previous {
            None => Anchor::Fixed,
            Some(bottom) => Anchor::At(PointRef::LineEnd(bottom)),
        };
        previous = Some(push_rectangle(&mut sk, i as f64 * WIDTH, 0.0, anchor));
    }
    sk
}

/// `count` independent rectangles, authored at the solution.
fn rectangle_field(count: usize) -> Sketch {
    let mut sk = Sketch::new();
    for i in 0..count {
        let (col, row) = ((i % 20) as f64, (i / 20) as f64);
        push_rectangle(&mut sk, col * 10.0, row * 10.0, Anchor::Fixed);
    }
    sk
}

/// The fixture knocked off its solution by `amp` (deterministic).
fn drifted(mut sk: Sketch, amp: f64) -> Sketch {
    sk.drift_free_params(amp, 7);
    sk
}

fn parameter_count(sk: &Sketch) -> usize {
    sk.params_snapshot().len()
}

fn bench_sketch(c: &mut Criterion) {
    // Cold solve of one connected component: every iteration factors the
    // whole system.
    let mut group = c.benchmark_group("sketch/solve_chain");
    group.sample_size(10);
    for count in [3usize, 11, 44] {
        let start = drifted(rectangle_chain(count), 0.05);
        assert!(start.clone().solve().converged, "chain {count} must solve");
        group.bench_with_input(
            BenchmarkId::from_parameter(parameter_count(&start)),
            &start,
            |b, start| {
                b.iter_batched(
                    || start.clone(),
                    |mut sk| black_box(sk.solve()),
                    BatchSize::LargeInput,
                )
            },
        );
    }
    group.finish();

    // Warm re-solve: the same component a small step away from its solution
    // (what a dimension edit or a drag frame looks like).
    let mut group = c.benchmark_group("sketch/resolve_chain");
    group.sample_size(10);
    for count in [3usize, 11, 44] {
        let mut solved = drifted(rectangle_chain(count), 0.05);
        assert!(solved.solve().converged);
        let start = drifted(solved, 1.0e-4);
        group.bench_with_input(
            BenchmarkId::from_parameter(parameter_count(&start)),
            &start,
            |b, start| {
                b.iter_batched(
                    || start.clone(),
                    |mut sk| black_box(sk.solve()),
                    BatchSize::LargeInput,
                )
            },
        );
    }
    group.finish();

    // Many small components: the decomposition carries the scaling.
    let mut group = c.benchmark_group("sketch/solve_field");
    group.sample_size(10);
    for count in [20usize, 200] {
        let start = drifted(rectangle_field(count), 0.05);
        assert!(start.clone().solve().converged, "field {count} must solve");
        group.bench_with_input(
            BenchmarkId::from_parameter(parameter_count(&start)),
            &start,
            |b, start| {
                b.iter_batched(
                    || start.clone(),
                    |mut sk| black_box(sk.solve()),
                    BatchSize::LargeInput,
                )
            },
        );
    }
    group.finish();

    // Rank-revealing diagnosis of one connected component.
    let mut group = c.benchmark_group("sketch/diagnose_chain");
    group.sample_size(10);
    for count in [3usize, 11, 44] {
        let sk = rectangle_chain(count);
        assert_eq!(sk.diagnose().dof, 0);
        group.bench_with_input(
            BenchmarkId::from_parameter(parameter_count(&sk)),
            &sk,
            |b, sk| b.iter(|| black_box(sk.diagnose())),
        );
    }
    group.finish();
}

/// Earlier points each truss point is tied to by a distance: its three
/// predecessors, which make it rigid, plus one brace further back.
///
/// The brace is a redundant (consistent) equation and it is there on
/// purpose. Without it the truss is a serial chain in which every point
/// hangs off the three before it, and errors amplify along the chain: from
/// a 0.05 perturbation the solver needs 5 iterations at 20 points, 17 at
/// 120, and does not converge within 100 at 160 or 260. That is a
/// conditioning case for the step control to win later, not a stable
/// timing fixture. With the brace every size converges in 4–5 iterations.
const TRUSS_TIES: [usize; 4] = [1, 2, 3, 5];

/// A braced 3D truss of points tied together by distances (see
/// [`TRUSS_TIES`]), its first five points fixed. One connected component of
/// `3 * points` parameters, authored at the solution.
fn truss(points: usize) -> SpaceSketch {
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
    let anchored = TRUSS_TIES.into_iter().max().unwrap_or(0);
    for &id in &ids[..anchored] {
        sk.fix_entity(id);
    }
    for i in anchored..points {
        for back in TRUSS_TIES {
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

fn bench_space(c: &mut Criterion) {
    // Cold solve of one connected 3D component.
    let mut group = c.benchmark_group("space/solve_truss");
    group.sample_size(10);
    for points in [20usize, 80, 260] {
        let mut start = truss(points);
        start.drift_free_params(0.05, 7);
        assert!(start.clone().solve().converged, "truss {points} must solve");
        let parameters = start.params_snapshot().len();
        group.bench_with_input(
            BenchmarkId::from_parameter(parameters),
            &start,
            |b, start| {
                b.iter_batched(
                    || start.clone(),
                    |mut sk| black_box(sk.solve()),
                    BatchSize::LargeInput,
                )
            },
        );
    }
    group.finish();

    // Rank-revealing diagnosis of the same component. Every brace is a
    // redundant row, so this is the diagnosis with dependencies to resolve
    // (the sketch chains above have none).
    let mut group = c.benchmark_group("space/diagnose_truss");
    group.sample_size(10);
    for points in [20usize, 80, 260] {
        let sk = truss(points);
        assert_eq!(sk.diagnose().dof, 0);
        let parameters = sk.params_snapshot().len();
        group.bench_with_input(BenchmarkId::from_parameter(parameters), &sk, |b, sk| {
            b.iter(|| black_box(sk.diagnose()))
        });
    }
    group.finish();
}

/// A system holding only the grounded body 0.
fn ground() -> ConstraintSystem {
    let mut system = ConstraintSystem::with_components(1);
    system.add_constraint(Mate::fixed_with_orientation(
        0,
        0,
        DVec3::ZERO,
        DQuat::IDENTITY,
    ));
    system
}

/// Ground plus `bodies` rigid bodies, each mated to the previous one by two
/// coincident pairs and one angle (six equations for its six freedoms).
fn assembly_chain(bodies: usize) -> ConstraintSystem {
    let mut system = ground();
    for i in 0..bodies {
        let base = system.num_components - 1;
        let body = system.num_components;
        system.num_components += 1;
        let id = 10 + i * 10;
        system.add_constraint(Mate::coincident(
            id,
            base,
            DVec3::new(3.0, 0.0, 0.0),
            body,
            DVec3::ZERO,
        ));
        system.add_constraint(Mate::coincident(
            id + 1,
            base,
            DVec3::new(5.0, 0.0, 0.0),
            body,
            DVec3::new(2.0, 0.0, 0.0),
        ));
        system.add_constraint(Mate::angle(
            id + 2,
            base,
            DVec3::new(0.0, 1.0, 0.0),
            body,
            DVec3::new(0.0, 0.0, 1.0),
            std::f64::consts::FRAC_PI_2,
        ));
    }
    system
}

/// Ground plus `bodies` rigid bodies, each pinned to the ground by three
/// coincident pairs (nine consistent equations for its six freedoms).
fn assembly_star(bodies: usize) -> ConstraintSystem {
    let mut system = ground();
    for i in 0..bodies {
        let body = system.num_components;
        system.num_components += 1;
        let t = DVec3::new((i % 5) as f64 * 8.0, (i / 5) as f64 * 6.0, 1.0);
        let id = 10 + i * 10;
        let pairs = [
            (t, DVec3::ZERO),
            (t + DVec3::new(2.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)),
            (t + DVec3::new(0.0, 3.0, 1.0), DVec3::new(0.0, 3.0, 1.0)),
        ];
        for (k, (world, local)) in pairs.into_iter().enumerate() {
            system.add_constraint(Mate::coincident(id + k, 0, world, body, local));
        }
    }
    system
}

/// A named assembly builder, parameterised by the number of bodies.
type AssemblyFixture = (&'static str, fn(usize) -> ConstraintSystem);

fn bench_assembly(c: &mut Criterion) {
    let solver = ConstraintSolver::new();
    let fixtures: [AssemblyFixture; 2] = [("chain", assembly_chain), ("star", assembly_star)];

    // Cold solve from the seed pose (every free body at the origin).
    for (name, build) in fixtures {
        let mut group = c.benchmark_group(format!("assembly/solve_{name}"));
        group.sample_size(10);
        for bodies in [5usize, 20, 40] {
            let system = build(bodies);
            assert!(solver.solve(&system).success, "{name} {bodies} must solve");
            group.bench_with_input(BenchmarkId::from_parameter(bodies), &system, |b, system| {
                b.iter(|| black_box(solver.solve(system)))
            });
        }
        group.finish();
    }

    // Rank-revealing diagnosis at the solved pose.
    let mut group = c.benchmark_group("assembly/diagnose_star");
    group.sample_size(10);
    for bodies in [5usize, 20, 40] {
        let system = assembly_star(bodies);
        let solved = solver.solve(&system);
        group.bench_with_input(BenchmarkId::from_parameter(bodies), &system, |b, system| {
            b.iter(|| {
                black_box(diagnose_assembly(
                    system,
                    &solved.positions,
                    &solved.rotations,
                ))
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_sketch, bench_space, bench_assembly);
criterion_main!(benches);
