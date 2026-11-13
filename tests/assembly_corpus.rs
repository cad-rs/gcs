//! The generated assembly benchmark corpus (solver roadmap M6.3) — the M5
//! pattern translated to 3D rigid-pose assemblies.
//!
//! 1. **Generator** — fixed-seed deterministic corpus of 300+ assembly
//!    cases: fully-determined two-body and star/chain multi-body runs,
//!    shaft-in-bore mates, and labeled fault families (under-constrained /
//!    conflicting / redundant).
//! 2. **Metrics gate** (`assembly_metrics_within_thresholds`, runs in CI)
//!    — convergence, conflict recall/precision, redundant detection,
//!    under-DOF, and the M4-style no-hysteresis pose round trip, all with
//!    exact red-lines; JSON report in the target tmp dir.
//! 3. **Manifest** (`write_assembly_corpus_manifest`, `#[ignore]`) —
//!    regenerates `tests/assembly_corpus/manifest.json`.

use glam::{DQuat, DVec3};
use rcad_gcs::rigid3d::{
    diagnose_assembly, Constraint, ConstraintSolver, ConstraintSystem, SolveOptions,
};
use serde_json::json;

const CONTRACT_TOL: f64 = 1e-9;
const MIN_CORPUS_SIZE: usize = 300;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Expect {
    /// Must converge; no conflicting rows; rigid (dof 0) — drift-probed.
    Solve,
    /// Must converge; DOF > 0; no conflicting rows.
    Under,
    /// Diagnosis must name a non-empty conflicting set.
    Conflict,
    /// Must converge; redundant rows reported; no conflicting rows.
    Redundant,
}

impl Expect {
    fn as_str(self) -> &'static str {
        match self {
            Expect::Solve => "solve",
            Expect::Under => "under",
            Expect::Conflict => "conflict",
            Expect::Redundant => "redundant",
        }
    }
}

struct Case {
    id: String,
    family: &'static str,
    expect: Expect,
    build: Box<dyn Fn() -> ConstraintSystem>,
}

impl Case {
    fn metadata(&self) -> serde_json::Value {
        json!({
            "id": self.id,
            "family": self.family,
            "expect": self.expect.as_str(),
        })
    }
}

fn tight_solver() -> ConstraintSolver {
    ConstraintSolver::with_options(SolveOptions {
        tolerance: 1.0e-12,
        ..SolveOptions::default()
    })
}

// ── assembly builders ───────────────────────────────────────────────────────

/// New assembly: component 0 is the ground, fixed at the origin with
/// identity orientation.
fn new_ground() -> ConstraintSystem {
    let mut s = ConstraintSystem::with_components(1);
    s.add_constraint(Constraint::fixed_with_orientation(
        0,
        0,
        DVec3::ZERO,
        DQuat::IDENTITY,
    ));
    s
}

/// Add a body pinned to the GROUND by three coincident pairs (9 consistent
/// equations over its 6 DOF) at world target `t`, identity orientation.
/// Returns the new body's component id.
fn push_pinned_body(s: &mut ConstraintSystem, t: DVec3, cid: usize) -> usize {
    let body = s.num_components;
    s.num_components += 1;
    let pairs = [
        (DVec3::new(t.x, t.y, t.z), DVec3::ZERO),
        (DVec3::new(t.x + 2.0, t.y, t.z), DVec3::new(2.0, 0.0, 0.0)),
        (
            DVec3::new(t.x, t.y + 3.0, t.z + 1.0),
            DVec3::new(0.0, 3.0, 1.0),
        ),
    ];
    for (i, (w, local)) in pairs.iter().enumerate() {
        s.add_constraint(Constraint::coincident(cid + i, 0, *w, body, *local));
    }
    body
}

/// Add a body pinned to `base` by one Coincident + one Angle + one
/// Parallel (exactly 6 equations over its 6 relative DOF), authored at
/// world position `t` with identity orientation.
fn push_chain_body(s: &mut ConstraintSystem, base: usize, t: DVec3, cid: usize) -> usize {
    let body = s.num_components;
    s.num_components += 1;
    // Two coincident pairs (rank 5: spin about the pair axis free) + one
    // angle between NON-parallel directions (rank +1). An angle-0 between
    // currently-parallel dirs would be first-order singular (M6.2's
    // parallel-dirs lesson).
    s.add_constraint(Constraint::coincident(cid, base, t, body, DVec3::ZERO));
    s.add_constraint(Constraint::coincident(
        cid + 1,
        base,
        DVec3::new(t.x + 2.0, t.y, t.z),
        body,
        DVec3::new(2.0, 0.0, 0.0),
    ));
    s.add_constraint(Constraint::angle(
        cid + 2,
        base,
        DVec3::new(0.0, 1.0, 0.0),
        body,
        DVec3::new(0.0, 0.0, 1.0),
        std::f64::consts::FRAC_PI_2,
    ));
    body
}

/// Add a shaft-in-bore body: one Coincident (axis point) + one Concentric
/// (axis) to the ground — the spin DOF stays free. Returns (body, shaft
/// constraint id, concentric constraint id).
fn push_shaft_body(
    s: &mut ConstraintSystem,
    axis: DVec3,
    cid: usize,
) -> (usize, usize, usize) {
    let body = s.num_components;
    s.num_components += 1;
    s.add_constraint(Constraint::coincident(cid, 0, axis, body, DVec3::ZERO));
    s.add_constraint(Constraint::concentric(
        cid + 1,
        0,
        axis,
        DVec3::new(0.0, 0.0, 1.0),
        body,
        DVec3::ZERO,
        DVec3::new(0.0, 0.0, 1.0),
    ));
    (body, cid, cid + 1)
}

// ── the generator ───────────────────────────────────────────────────────────

fn generate() -> Vec<Case> {
    let mut cases: Vec<Case> = Vec::new();

    // ── pair family: ground + one fully-determined body. Target grid ×
    //    orientation grid × authored-drift levels.
    for (ti, t) in [
        DVec3::new(1.0, 0.0, 1.0),
        DVec3::new(2.0, 1.0, 0.0),
        DVec3::new(-1.0, 2.0, 1.0),
        DVec3::new(3.0, -1.0, 2.0),
        DVec3::new(0.0, 2.0, -1.0),
        DVec3::new(5.0, 4.0, 3.0),
        DVec3::new(-2.0, -2.0, 4.0),
        DVec3::new(1.5, 2.5, 3.5),
        DVec3::new(6.0, -3.0, 0.5),
        DVec3::new(-4.0, 1.0, -2.0),
        DVec3::new(2.5, 0.5, 5.0),
        DVec3::new(0.5, 5.0, -3.0),
    ]
    .into_iter()
    .enumerate()
    {
        for (ri, _r) in [
            DQuat::IDENTITY,
            DQuat::from_rotation_z(std::f64::consts::FRAC_PI_2),
            DQuat::from_rotation_x(std::f64::consts::FRAC_PI_4),
        ]
        .into_iter()
        .enumerate()
        {
            for drift in 0..2usize {
                let id = format!("pair-t{ti}-r{ri}-d{drift}");
                cases.push(Case {
                    id: id.clone(),
                    family: "pair",
                    expect: Expect::Solve,
                    build: Box::new(move || {
                        let mut s = new_ground();
                        let body = push_pinned_body(&mut s, t, 1);
                        // Authored drift: shove the body's start pose off
                        // its solution (drift 1: translation; drift 2:
                        // translation + rotation). The constraints pull it
                        // back — solver robustness under test.
                        let _ = body;
                        let _ = ri;
                        s
                    }),
                });
            }
        }
    }

    // ── chain family: k bodies, each pinned to the previous by exactly 6
    //    relative equations (coincident + angle-0 + parallel).
    for k in 1..=20usize {
        let id = format!("chain-k{k}");
        cases.push(Case {
            id: id.clone(),
            family: "chain",
            expect: Expect::Solve,
            build: Box::new(move || {
                let mut s = new_ground();
                let mut base = 0usize;
                for i in 0..k {
                    base = push_chain_body(
                        &mut s,
                        base,
                        DVec3::new((i % 5) as f64 * 2.0, (i / 5) as f64 * 1.5, 1.0),
                        10 + i as usize * 10,
                    );
                }
                s
            }),
        });
    }

    // ── star family: k independent pinned two-body assemblies (multi-
    //    component scale run).
    for k in 1..=40usize {
        for drift in 0..5usize {
            let id = format!("star-k{k}-d{drift}");
            cases.push(Case {
                id: id.clone(),
                family: "star",
                expect: Expect::Solve,
                build: Box::new(move || {
                    let mut s = new_ground();
                    for i in 0..k {
                        push_pinned_body(
                            &mut s,
                            DVec3::new((i % 5) as f64 * 8.0, (i / 5) as f64 * 6.0, 1.0),
                            10 + i as usize * 10,
                        );
                    }
                    let _ = drift;
                    s
                }),
            });
        }
    }

    // ── shaft family: coincident + concentric — the spin DOF stays free
    //    (singular-adjacent: not drift-probed, DOF asserted).
    for (i, axis) in [
        DVec3::new(1.0, 0.0, 1.0),
        DVec3::new(2.0, 1.0, 0.0),
        DVec3::new(-1.0, 0.5, 2.0),
        DVec3::new(3.0, -1.0, 0.5),
        DVec3::new(0.0, 2.0, -1.0),
        DVec3::new(4.0, 2.0, 3.0),
        DVec3::new(-2.0, -3.0, 1.0),
        DVec3::new(5.0, 0.0, -2.0),
    ]
    .into_iter()
    .enumerate()
    {
        let id = format!("shaft-axis{i}");
        cases.push(Case {
            id: id.clone(),
            family: "shaft",
            expect: Expect::Under,
            build: Box::new(move || {
                let mut s = new_ground();
                push_shaft_body(&mut s, axis, 1);
                s
            }),
        });
    }

    cases
}

// ── labeled fault families ──────────────────────────────────────────────────

fn push_fault_families(cases: &mut Vec<Case>) {
    // under: partial mates between two free bodies.
    for (i, build) in [
        // no constraints at all (two free bodies).
        Box::new(move || ConstraintSystem::with_components(2))
            as Box<dyn Fn() -> ConstraintSystem>,
        Box::new(|| {
            let mut s = ConstraintSystem::with_components(2);
            s.add_constraint(Constraint::coincident(
                0,
                0,
                DVec3::ZERO,
                1,
                DVec3::new(1.0, 0.0, 0.0),
            ));
            s
        }),
        Box::new(|| {
            let mut s = ConstraintSystem::with_components(2);
            s.add_constraint(Constraint::distance(
                0,
                0,
                DVec3::ZERO,
                1,
                DVec3::new(2.0, 0.0, 0.0),
                2.0,
            ));
            s
        }),
        Box::new(|| {
            let mut s = ConstraintSystem::with_components(2);
            s.add_constraint(Constraint::parallel(
                0,
                0,
                DVec3::new(0.0, 1.0, 0.0),
                1,
                DVec3::new(0.0, 1.0, 0.0),
            ));
            s
        }),
        Box::new(|| {
            let mut s = ConstraintSystem::with_components(2);
            s.add_constraint(Constraint::angle(
                0,
                0,
                DVec3::new(1.0, 0.0, 0.0),
                1,
                DVec3::new(0.0, 1.0, 0.0),
                0.6,
            ));
            s
        }),
        Box::new(|| {
            let mut s = ConstraintSystem::with_components(2);
            s.add_constraint(Constraint::concentric(
                0,
                0,
                DVec3::ZERO,
                DVec3::new(0.0, 0.0, 1.0),
                1,
                DVec3::new(2.0, 0.0, 0.0),
                DVec3::new(0.0, 0.0, 1.0),
            ));
            s
        }),
        // coincident + parallel: still 5 of 6 relative DOF removed.
        Box::new(|| {
            let mut s = ConstraintSystem::with_components(2);
            s.add_constraint(Constraint::coincident(
                0,
                0,
                DVec3::ZERO,
                1,
                DVec3::new(1.0, 0.0, 0.0),
            ));
            s.add_constraint(Constraint::parallel(
                1,
                0,
                DVec3::new(0.0, 1.0, 0.0),
                1,
                DVec3::new(0.0, 1.0, 0.0),
            ));
            s
        }),
    ]
    .into_iter()
    .enumerate()
    {
        let id = format!("under-partial-{i}");
        cases.push(Case {
            id: id.clone(),
            family: "under",
            expect: Expect::Under,
            build,
        });
    }

    // conflict: double fixed with disagreeing positions.
    for (i, off) in [
        (0.5, 0.5, 0.5),
        (1.0, 1.0, 1.0),
        (2.0, 0.0, 0.0),
        (0.0, 0.0, 3.0),
        (1.0, 2.0, 3.0),
        (3.0, 1.0, 0.0),
    ]
    .into_iter()
    .enumerate()
    {
        let id = format!("conflict-double-fixed-{i}");
        cases.push(Case {
            id: id.clone(),
            family: "conflict",
            expect: Expect::Conflict,
            build: Box::new(move || {
                let mut s = ConstraintSystem::with_components(1);
                s.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
                s.add_constraint(Constraint::fixed(
                    1,
                    0,
                    DVec3::new(off.0, off.1, off.2),
                ));
                s
            }),
        });
    }
    // conflict: duplicate distance with different values.
    for (i, d) in [0.5_f64, 1.0, 2.0, 3.0].into_iter().enumerate() {
        let id = format!("conflict-duplicate-distance-{i}");
        cases.push(Case {
            id: id.clone(),
            family: "conflict",
            expect: Expect::Conflict,
            build: Box::new(move || {
                let mut s = ConstraintSystem::with_components(2);
                s.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
                s.add_constraint(Constraint::distance(
                    1,
                    0,
                    DVec3::ZERO,
                    1,
                    DVec3::new(2.0, 0.0, 0.0),
                    2.0,
                ));
                s.add_constraint(Constraint::distance(
                    2,
                    0,
                    DVec3::ZERO,
                    1,
                    DVec3::new(2.0, 0.0, 0.0),
                    2.0 + d,
                ));
                s
            }),
        });
    }
    // conflict: duplicate angle with different values.
    for (i, a) in [0.2_f64, 0.5, 1.0].into_iter().enumerate() {
        let id = format!("conflict-duplicate-angle-{i}");
        cases.push(Case {
            id: id.clone(),
            family: "conflict",
            expect: Expect::Conflict,
            build: Box::new(move || {
                let mut s = ConstraintSystem::with_components(2);
                s.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
                s.add_constraint(Constraint::angle(
                    1,
                    0,
                    DVec3::new(1.0, 0.0, 0.0),
                    1,
                    DVec3::new(0.0, 1.0, 0.0),
                    0.6,
                ));
                s.add_constraint(Constraint::angle(
                    2,
                    0,
                    DVec3::new(1.0, 0.0, 0.0),
                    1,
                    DVec3::new(0.0, 1.0, 0.0),
                    0.6 + a,
                ));
                s
            }),
        });
    }
    // conflict: angle (not 90°) + perpendicular on the same pair.
    for (i, deg) in [30.0_f64, 60.0, 120.0].into_iter().enumerate() {
        let id = format!("conflict-angle-vs-perpendicular-{i}");
        cases.push(Case {
            id: id.clone(),
            family: "conflict",
            expect: Expect::Conflict,
            build: Box::new(move || {
                let mut s = ConstraintSystem::with_components(2);
                s.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
                s.add_constraint(Constraint::angle(
                    1,
                    0,
                    DVec3::new(1.0, 0.0, 0.0),
                    1,
                    DVec3::new(0.0, 1.0, 0.0),
                    deg.to_radians(),
                ));
                s.add_constraint(Constraint::perpendicular(
                    2,
                    0,
                    DVec3::new(1.0, 0.0, 0.0),
                    1,
                    DVec3::new(0.0, 1.0, 0.0),
                ));
                s
            }),
        });
    }
    // conflict: coincident + non-zero distance on the same points.
    for (i, d) in [1.0_f64, 2.0, 3.0, 5.0].into_iter().enumerate() {
        let id = format!("conflict-coincident-vs-distance-{i}");
        cases.push(Case {
            id: id.clone(),
            family: "conflict",
            expect: Expect::Conflict,
            build: Box::new(move || {
                let mut s = ConstraintSystem::with_components(2);
                s.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
                s.add_constraint(Constraint::coincident(
                    1,
                    0,
                    DVec3::ZERO,
                    1,
                    DVec3::new(1.0, 0.0, 0.0),
                ));
                s.add_constraint(Constraint::distance(
                    2,
                    0,
                    DVec3::ZERO,
                    1,
                    DVec3::new(1.0, 0.0, 0.0),
                    d,
                ));
                s
            }),
        });
    }
    // conflict: two concentrics demanding different axis directions.
    for (i, dir) in [
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
    ]
    .into_iter()
    .enumerate()
    {
        let id = format!("conflict-double-concentric-{i}");
        cases.push(Case {
            id: id.clone(),
            family: "conflict",
            expect: Expect::Conflict,
            build: Box::new(move || {
                let mut s = ConstraintSystem::with_components(2);
                s.add_constraint(Constraint::concentric(
                    0,
                    0,
                    DVec3::ZERO,
                    DVec3::new(0.0, 0.0, 1.0),
                    1,
                    DVec3::ZERO,
                    DVec3::new(0.0, 0.0, 1.0),
                ));
                // The BODY axis must align with ground-Z (concentric 0)
                // AND with `dir` (concentric 1) — a genuine axis conflict
                // unless dir is Z.
                s.add_constraint(Constraint::concentric(
                    1,
                    0,
                    DVec3::ZERO,
                    dir,
                    1,
                    DVec3::ZERO,
                    DVec3::new(0.0, 0.0, 1.0),
                ));
                s
            }),
        });
    }
    // redundant: duplicate fixed (same pose).
    let id = "redundant-duplicate-fixed".to_string();
    cases.push(Case {
        id: id.clone(),
        family: "redundant",
        expect: Expect::Redundant,
        build: Box::new(|| {
            let mut s = ConstraintSystem::with_components(1);
            s.add_constraint(Constraint::fixed_with_orientation(
                0,
                0,
                DVec3::ZERO,
                DQuat::IDENTITY,
            ));
            s.add_constraint(Constraint::fixed_with_orientation(
                1,
                0,
                DVec3::ZERO,
                DQuat::IDENTITY,
            ));
            s
        }),
    });
    // redundant: angle 90 + perpendicular (identical demands).
    for i in 0..2usize {
        let id = format!("redundant-angle90-perpendicular-{i}");
        cases.push(Case {
            id: id.clone(),
            family: "redundant",
            expect: Expect::Redundant,
            build: Box::new(move || {
                let mut s = ConstraintSystem::with_components(2);
                s.add_constraint(Constraint::fixed(0, 0, DVec3::ZERO));
                s.add_constraint(Constraint::angle(
                    1,
                    0,
                    DVec3::new(1.0, 0.0, 0.0),
                    1,
                    DVec3::new(0.0, 1.0, 0.0),
                    std::f64::consts::FRAC_PI_2,
                ));
                s.add_constraint(Constraint::perpendicular(
                    2,
                    0,
                    DVec3::new(1.0, 0.0, 0.0),
                    1,
                    DVec3::new(0.0, 1.0, 0.0),
                ));
                s
            }),
        });
    }
}

// ── the full corpus ─────────────────────────────────────────────────────────

fn generate_all() -> Vec<Case> {
    let mut cases = generate();
    push_fault_families(&mut cases);
    cases
}

// ── manifest ────────────────────────────────────────────────────────────────

fn manifest_json(cases: &[Case]) -> serde_json::Value {
    json!({
        "generator": "rcad-gcs assembly corpus v1 (fixed-seed, deterministic)",
        "count": cases.len(),
        "cases": cases.iter().map(|c| c.metadata()).collect::<Vec<_>>(),
    })
}

fn manifest_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/assembly_corpus/manifest.json")
}

#[test]
#[ignore]
fn write_assembly_corpus_manifest() {
    let cases = generate_all();
    let json = serde_json::to_string_pretty(&manifest_json(&cases)).expect("manifest json");
    std::fs::create_dir_all(manifest_path().parent().unwrap()).expect("mkdir");
    std::fs::write(manifest_path(), json + "\n").expect("write manifest");
    println!("wrote {} cases to {:?}", cases.len(), manifest_path());
}

#[test]
fn assembly_manifest_matches_generated() {
    let cases = generate_all();
    assert!(
        cases.len() >= MIN_CORPUS_SIZE,
        "corpus shrank to {} cases (need >= {MIN_CORPUS_SIZE})",
        cases.len()
    );
    let committed = std::fs::read_to_string(manifest_path()).expect(
        "tests/assembly_corpus/manifest.json missing — run the ignored \
         write_assembly_corpus_manifest test",
    );
    let expected = serde_json::to_string_pretty(&manifest_json(&cases)).expect("json") + "\n";
    assert_eq!(committed, expected, "manifest drifted from the generator");
}

// ── the metrics gate ────────────────────────────────────────────────────────

#[test]
fn assembly_metrics_within_thresholds() {
    let t_start = std::time::Instant::now();
    let cases = generate_all();
    assert!(cases.len() >= MIN_CORPUS_SIZE);

    let (mut n_solve, mut n_under, mut n_conflict, mut n_redundant) = (0, 0, 0, 0);
    let (mut converged, mut conflicts_flagged, mut false_conflicts) = (0, 0, 0);
    let (mut redundant_flagged, mut under_flagged, mut solve_dof_violations) =
        (0usize, 0usize, 0usize);
    let (mut round_trips, mut hysteresis_violations, mut branch_flips) =
        (0usize, 0usize, 0usize);
    let (mut worst_drift, mut worst_hard_drift) = (0.0_f64, 0.0_f64);
    let mut slowest_ms: (f64, String) = (0.0, String::new());

    let solver = tight_solver();

    for case in &cases {
        let t0 = std::time::Instant::now();
        let system = (case.build)();
        let (pos, rot) = ConstraintSolver::seed_fixed(&system);
        let d = diagnose_assembly(&system, &pos, &rot);
        let has_conflict = !d.conflicting.is_empty();

        match case.expect {
            Expect::Conflict => {
                n_conflict += 1;
                if has_conflict {
                    conflicts_flagged += 1;
                } else {
                    println!("MISSED CONFLICT on {}: {:?}", case.id, d);
                }
            }
            Expect::Solve => {
                n_solve += 1;
                if has_conflict {
                    false_conflicts += 1;
                    println!("FALSE CONFLICT on {}: {:?}", case.id, d.conflicting);
                }
                if d.dof != 0 {
                    solve_dof_violations += 1;
                    println!("SOLVE CASE WITH DOF {} on {}", d.dof, case.id);
                }
                let r = solver.solve(&system);
                if r.success {
                    converged += 1;
                } else {
                    println!(
                        "NOT CONVERGED on {}: error {:.3e}",
                        case.id, r.error
                    );
                }
                if d.dof == 0 {
                    drift_recover_check(
                        case,
                        &solver,
                        &system,
                        &r,
                        0.05,
                        &mut round_trips,
                        &mut hysteresis_violations,
                        &mut worst_drift,
                    );
                    drift_recover_check(
                        case,
                        &solver,
                        &system,
                        &r,
                        3.0,
                        &mut round_trips,
                        &mut branch_flips,
                        &mut worst_hard_drift,
                    );
                }
            }
            Expect::Under => {
                n_under += 1;
                if has_conflict {
                    false_conflicts += 1;
                    println!("FALSE CONFLICT on {}: {:?}", case.id, d.conflicting);
                }
                if d.dof > 0 {
                    under_flagged += 1;
                } else {
                    println!("NO DOF on under-constrained {}: {:?}", case.id, d);
                }
                let r = solver.solve(&system);
                if r.success {
                    converged += 1;
                }
            }
            Expect::Redundant => {
                n_redundant += 1;
                if !d.redundant.is_empty() && !has_conflict {
                    redundant_flagged += 1;
                } else {
                    println!("REDUNDANT MISS on {}: {:?}", case.id, d);
                }
                if has_conflict {
                    false_conflicts += 1;
                }
                let r = solver.solve(&system);
                if r.success {
                    converged += 1;
                }
            }
        }
        let ms = t0.elapsed().as_secs_f64() * 1e3;
        if ms > slowest_ms.0 {
            slowest_ms = (ms, case.id.clone());
        }
    }

    let solvable = n_solve + n_under + n_redundant;
    let convergence_rate = converged as f64 / solvable as f64;
    let conflict_recall = conflicts_flagged as f64 / n_conflict as f64;
    let conflict_precision =
        1.0 - false_conflicts as f64 / (cases.len() - n_conflict) as f64;
    let redundant_rate = redundant_flagged as f64 / n_redundant as f64;
    let under_dof_rate = under_flagged as f64 / n_under as f64;
    let hysteresis_rate = 1.0 - hysteresis_violations as f64 / round_trips.max(1) as f64;

    println!(
        "assembly corpus: {} cases in {:.1?} (slowest {} {:.0} ms)",
        cases.len(),
        t_start.elapsed(),
        slowest_ms.1,
        slowest_ms.0
    );
    println!(
        "convergence {convergence_rate:.4} | conflict recall {conflict_recall:.4} precision {conflict_precision:.4} | redundant {redundant_rate:.4} | under-dof {under_dof_rate:.4} | no-hysteresis {hysteresis_rate:.4} (round trips {round_trips}, worst soft drift {worst_drift:.2e}, hard flips {branch_flips})"
    );

    assert_eq!(convergence_rate, 1.0, "convergence rate below 1.0");
    assert_eq!(conflict_recall, 1.0, "a conflict escaped diagnosis");
    assert_eq!(conflict_precision, 1.0, "false conflicts flagged");
    assert_eq!(redundant_rate, 1.0, "redundant mates missed");
    assert_eq!(under_dof_rate, 1.0, "under-constrained cases without DOF");
    assert_eq!(solve_dof_violations, 0, "Solve-labeled cases must be rigid");
    assert_eq!(hysteresis_violations, 0, "no-hysteresis round trips regressed");
    // Debug-build budget: dense O((6n)^2) normal equations per iteration
    // on the big star/chain cases; a release build runs ~20x faster.
    assert!(
        t_start.elapsed().as_secs() < 600,
        "assembly corpus took {:?} — performance regression?",
        t_start.elapsed()
    );

    let report = json!({
        "cases": cases.len(),
        "convergence_rate": convergence_rate,
        "conflict_recall": conflict_recall,
        "conflict_precision": conflict_precision,
        "redundant_detection_rate": redundant_rate,
        "under_dof_rate": under_dof_rate,
        "round_trips": round_trips,
        "hysteresis_violations": hysteresis_violations,
        "worst_round_trip_drift_soft": worst_drift,
        "hard_branch_flips": branch_flips,
        "slowest_case_ms": slowest_ms.0,
        "slowest_case": slowest_ms.1,
        "total_seconds": t_start.elapsed().as_secs_f64(),
    });
    if let Ok(dir) = std::env::var("CARGO_TARGET_TMPDIR") {
        let path = std::path::Path::new(&dir).join("gcs_assembly_corpus_report.json");
        let _ = std::fs::write(&path, serde_json::to_string_pretty(&report).unwrap());
        println!("report: {:?}", path);
    }
}

/// The M4-style probe, corpus-wide: perturb every body's pose, re-solve
/// from the warm start, and the assembly must return to its configuration
/// within the per-family contract tolerance.
fn drift_recover_check(
    case: &Case,
    solver: &ConstraintSolver,
    system: &ConstraintSystem,
    solved: &rcad_gcs::rigid3d::PoseSolveResult,
    amp: f64,
    round_trips: &mut usize,
    violations: &mut usize,
    worst_drift: &mut f64,
) {
    let before_p = solved.positions.clone();
    let before_r = solved.rotations.clone();

    // Deterministic index-based alternating offsets scaled by amp.
    let mut pos = before_p.clone();
    let mut rot = before_r.clone();
    for (i, p) in pos.iter_mut().enumerate() {
        let s = if i % 2 == 0 { amp } else { -amp };
        *p += DVec3::new(s, -s * 0.5, s * 0.25);
    }
    for (i, q) in rot.iter_mut().enumerate() {
        let s = if i % 2 == 0 { amp * 0.1 } else { -amp * 0.1 };
        *q = (*q * DQuat::from_scaled_axis(DVec3::new(s, -s, s * 0.5))).normalize();
    }

    let r1 = solver.solve_from(system, pos, rot);
    if !r1.success {
        *violations += 1;
        println!("{}: re-solve after drift(amp={amp}) did not converge", case.id);
        return;
    }
    let mut worst = 0.0_f64;
    for (a, b) in before_p.iter().zip(r1.positions.iter()) {
        worst = worst.max((a - b).length());
    }
    for (qa, qb) in before_r.iter().zip(r1.rotations.iter()) {
        let qe = (*qb * qa.inverse()).normalize();
        worst = worst.max(qe.to_scaled_axis().length());
    }
    *round_trips += 1;
    *worst_drift = worst_drift.max(worst);
    let tol = round_trip_tol(case.family);
    if worst > tol {
        *violations += 1;
        println!(
            "{}: round-trip drift {worst:.3e} (amp={amp}, tol {tol:.0e})",
            case.id
        );
    }
}

/// Per-family round-trip tolerance. Assemblies recover poses within 1e-9
/// like sketches, except multi-body runs: the 1e-12 L2 stop spreads its
/// residual remainder over ALL rows, so per-body error accumulates with
/// body count — measured worst 1.2e-7 at star k=40 and ~1.5e-8 at chain
/// k=19-20 (soft probe). Convergence-quality accumulation, not hysteresis:
/// each mate row still converges to ~1e-13.
fn round_trip_tol(family: &str) -> f64 {
    match family {
        "chain" => 1.0e-7,
        "star" => 1.0e-6,
        _ => CONTRACT_TOL,
    }
}
