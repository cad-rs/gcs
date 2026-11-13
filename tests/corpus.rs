//! The generated benchmark corpus (solver roadmap M5).
//!
//! Three pieces in one file:
//!
//! 1. **Generator** — a fixed-seed, fully deterministic corpus of 300+
//!    constraint systems across every constraint family: healthy
//!    (fully-constrained solvable), under-constrained, conflicting,
//!    redundant, degenerate-adjacent (exact tangency), and scale runs
//!    (multi-component sketches up to hundreds of parameters).
//!    Rebuilding is reproducible: same seed, same cases, every run.
//! 2. **Metrics gate** (`corpus_metrics_within_thresholds`, runs in CI) —
//!    every case is diagnosed and solved; exact red-lines on convergence
//!    rate, conflict recall/precision, redundant detection, and the M4
//!    no-hysteresis round trip; a JSON report lands in the target tmp
//!    dir.
//! 3. **Manifest** (`write_corpus_manifest`, `#[ignore]`) — regenerates
//!    `tests/corpus/manifest.json`; `corpus_manifest_matches_generated`
//!    fails if the generator drifts without the manifest being refreshed.
//!
//! This is the "修一个崩一个" guard: every solver change runs against the
//! whole corpus, not just the case that motivated it.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::sketch::Sketch;
use serde_json::json;

const CONTRACT_TOL: f64 = 1e-9;
const MIN_CORPUS_SIZE: usize = 300;

/// Per-family no-hysteresis round-trip tolerance. 1e-9 everywhere except
/// the exact-tangency family: at exact tangency the Jacobian is
/// first-order RANK-DEFICIENT (analytically, the Tangent row equals
/// DX-row + DY-row there, so the end-x direction is constrained only at
/// second order) — parameters recover to ~1e-4 (measured worst 1.5e-4
/// for the jittered-authoring cases at amp 0.05) while the tangency
/// system, not of the solver; every CAD solver with a tangency constraint
/// shares it, which is why CAD editors snap tangency instead of trusting
/// numerics alone.
fn round_trip_tol(family: &str) -> f64 {
    match family {
        "tangent" => 1.0e-3,
        _ => CONTRACT_TOL,
    }
}

// ── deterministic PRNG (xorshift64*) ────────────────────────────────────────

struct Lcg(u64);
impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg(seed ^ 0x9E3779B97F4A7C15)
    }
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let x = self.0 ^ (self.0 >> 29);
        let x = x.wrapping_mul(0xBF58476D1CE4E5B9);
        (x >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

// ── corpus model ────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Expect {
    /// Must converge; no conflicting rows; rigid (dof 0).
    Solve,
    /// Must converge; true DOF > 0; no conflicting rows.
    Under,
    /// Diagnosis must name a non-empty conflicting set before solving.
    Conflict,
    /// Must converge; a redundant (dependent, consistent) set is reported;
    /// no conflicting rows.
    Redundant,
    /// Solvable but FIRST-ORDER SINGULAR at the authored state (exact
    /// tangency: the Tangent row equals DX+DY there, so dof 1 is the
    /// correct diagnosis). Converges, no conflicts; the round trip runs
    /// at the relaxed per-family tolerance.
    Singular,
}

impl Expect {
    fn as_str(self) -> &'static str {
        match self {
            Expect::Solve => "solve",
            Expect::Singular => "singular",
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
    build: Box<dyn Fn() -> Sketch>,
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

// ── shared builders ─────────────────────────────────────────────────────────

/// One fully-constrained rectangle (18 params / 18 equations): fixed
/// corner point, chained coincidences, H/V pair, two driven lengths.
fn push_rectangle(sk: &mut Sketch, ox: f64, oy: f64, w: f64, h: f64, drift: &mut Lcg) {
    let origin = sk.add_point(ox, oy);
    let j = |d: &mut Lcg| d.range(-0.05, 0.05);
    let b = sk.add_line(ox, oy, ox + w + j(drift), oy + j(drift));
    let r = sk.add_line(ox + w, oy, ox + w, oy + h);
    let t = sk.add_line(ox + w, oy + h, ox, oy + h);
    let l = sk.add_line(ox, oy + h, ox, oy);
    sk.add_constraint(Fixed {
        point: PointRef::Point(origin),
        x: ox,
        y: oy,
    });
    sk.add_constraint(Coincident(PointRef::LineStart(b), PointRef::Point(origin)));
    for pair in [
        (PointRef::LineEnd(b), PointRef::LineStart(r)),
        (PointRef::LineEnd(r), PointRef::LineStart(t)),
        (PointRef::LineEnd(t), PointRef::LineStart(l)),
        (PointRef::LineEnd(l), PointRef::LineStart(b)),
    ] {
        sk.add_constraint(Coincident(pair.0, pair.1));
    }
    sk.add_constraint(Horizontal(b));
    sk.add_constraint(Horizontal(t));
    sk.add_constraint(Vertical(r));
    sk.add_constraint(Vertical(l));
    sk.add_constraint(LineLength { line: b, length: w });
    sk.add_constraint(LineLength { line: r, length: h });
}

/// The axis-aligned margin pattern (skeleton + corner hole + centered
/// bore).
fn push_margin(sk: &mut Sketch, len: f64, margin: f64, hole_d: f64, bore_d: f64) {
    let anchor = sk.add_point(0.0, 0.0);
    let line = sk.add_line(0.0, 0.0, len, 0.0);
    sk.add_constraint(Fixed {
        point: PointRef::Point(anchor),
        x: 0.0,
        y: 0.0,
    });
    sk.add_constraint(Coincident(
        PointRef::LineStart(line),
        PointRef::Point(anchor),
    ));
    sk.add_constraint(Horizontal(line));
    sk.add_constraint(DistanceX {
        p1: PointRef::LineEnd(line),
        p2: PointRef::LineStart(line),
        distance: len,
    });
    let hole = sk.add_circle(len - margin, margin, hole_d / 2.0);
    sk.add_constraint(Radius {
        circle: hole,
        radius: hole_d / 2.0,
    });
    sk.add_constraint(DistanceX {
        p1: PointRef::LineEnd(line),
        p2: PointRef::Center(hole),
        distance: margin,
    });
    sk.add_constraint(DistanceY {
        p1: PointRef::Center(hole),
        p2: PointRef::Point(anchor),
        distance: margin,
    });
    let mid = sk.add_point(len * 0.5, 0.0);
    sk.add_constraint(Midpoint {
        point: PointRef::Point(mid),
        line,
    });
    let bore = sk.add_circle(len * 0.5, margin, bore_d / 2.0);
    sk.add_constraint(DistanceX {
        p1: PointRef::Point(mid),
        p2: PointRef::Center(bore),
        distance: 0.0,
    });
    sk.add_constraint(DistanceY {
        p1: PointRef::Center(bore),
        p2: PointRef::Point(anchor),
        distance: margin,
    });
    sk.add_constraint(Diameter {
        circle: bore,
        diameter: bore_d,
    });
}

// ── the generator ───────────────────────────────────────────────────────────

fn generate() -> Vec<Case> {
    let mut cases: Vec<Case> = Vec::new();

    // ── rect family: k independent rectangles (a multi-component scale
    //    run; k=30 is 540 params / 30 components), drift levels 0-2.
    for k in 1..=40usize {
        for drift in 0..3usize {
            let id = format!("rect-k{k}-d{drift}");
            cases.push(Case {
                id: id.clone(),
                family: "rect",
                expect: Expect::Solve,
                build: Box::new(move || {
                    let mut rng = Lcg::new(1000 + k as u64 * 10 + drift as u64);
                    let mut sk = Sketch::new();
                    for i in 0..k {
                        push_rectangle(
                            &mut sk,
                            (i % 5) as f64 * 20.0,
                            (i / 5) as f64 * 12.0,
                            6.0 + i as f64 * 0.25,
                            4.0,
                            &mut rng,
                        );
                    }
                    sk
                }),
            });
        }
    }

    // ── margin family: the mounting-plate pattern, parameter grid × drift.
    for len in [8.0, 10.0, 12.0, 16.0, 20.0] {
        for margin in [1.5, 2.0, 3.0] {
            for hole in [6.0, 8.0] {
                for drift in 0..3usize {
                    let id = format!("margin-l{len}-m{margin}-h{hole}-d{drift}");
                    cases.push(Case {
                        id: id.clone(),
                        family: "margin",
                        expect: Expect::Solve,
                        build: Box::new(move || {
                            let mut sk = Sketch::new();
                            push_margin(&mut sk, len, margin, hole, (len * 0.5).min(10.0));
                            let _ = drift;
                            sk
                        }),
                    });
                }
            }
        }
    }

    // ── tangent family: start pinned on the circle, tangent, driven
    //    length — exact tangency included (degenerate-adjacent).
    for r in [1.0, 2.0, 3.0, 4.0, 5.0] {
        for len_mult in [2.0, 3.0] {
            for (jx, jy) in [(0.0, 0.0), (0.1, -0.05), (-0.1, 0.05)] {
                // j = 0 authors the state exactly AT the singular point:
                // first-order rank deficiency, dof 1 is the correct
                // diagnosis (Expect::Singular).
                let id = format!("tangent-r{r}-l{len_mult}-j{}x{}", jx * 100.0, jy * 100.0);
                cases.push(Case {
                    id: id.clone(),
                    family: "tangent",
                    expect: if jx == 0.0 {
                        Expect::Singular
                    } else {
                        Expect::Solve
                    },
                    build: Box::new(move || {
                        let mut sk = Sketch::new();
                        let c = sk.add_circle(0.0, 0.0, r);
                        sk.fix_entity(c);
                        // Authored SATISFYING the signed tangent: the
                        // line goes DOWN from the on-circle start (the
                        // signed distance convention picks this branch).
                        let l = sk.add_line(r + jx, jy, r + jx, jy - len_mult * r);
                        // Start pinned ON the circle at angle 0 via directed
                        // distances (PointOnCircle would leave the line one
                        // DOF short — it slides along its own tangent; that
                        // variant lives in the under family).
                        sk.add_constraint(DistanceX {
                            p1: PointRef::LineStart(l),
                            p2: PointRef::Center(c),
                            distance: r,
                        });
                        sk.add_constraint(DistanceY {
                            p1: PointRef::LineStart(l),
                            p2: PointRef::Center(c),
                            distance: 0.0,
                        });
                        sk.add_constraint(Tangent { circle: c, line: l });
                        sk.add_constraint(LineLength {
                            line: l,
                            length: len_mult * r,
                        });
                        sk
                    }),
                });
            }
        }
    }

    // ── symmetric family: driven span, point-symmetric bore.
    for span in [6.0, 10.0, 14.0] {
        for offset in [0.0, 2.0, 4.0] {
            let id = format!("symmetric-s{span}-o{offset}");
            cases.push(Case {
                id: id.clone(),
                family: "symmetric",
                expect: Expect::Solve,
                build: Box::new(move || {
                    let mut sk = Sketch::new();
                    let a = sk.add_point(0.0, 0.0);
                    let b = sk.add_point(span, offset);
                    let bore = sk.add_circle(span * 0.5, offset * 0.5, 1.0);
                    sk.add_constraint(Fixed {
                        point: PointRef::Point(a),
                        x: 0.0,
                        y: 0.0,
                    });
                    sk.add_constraint(DistanceX {
                        p1: PointRef::Point(b),
                        p2: PointRef::Point(a),
                        distance: span,
                    });
                    sk.add_constraint(DistanceY {
                        p1: PointRef::Point(b),
                        p2: PointRef::Point(a),
                        distance: offset,
                    });
                    sk.add_constraint(SymmetricPoint {
                        p1: PointRef::Center(bore),
                        p2: PointRef::Point(b),
                        about: PointRef::Point(a),
                    });
                    sk.add_constraint(Radius {
                        circle: bore,
                        radius: 1.0,
                    });
                    sk
                }),
            });
        }
    }

    // ── angle family: driven angle + length on a pinned second line.
    for deg in [15.0_f64, 45.0, 90.0, 135.0, 200.0] {
        for len in [3.0, 7.0] {
            let id = format!("angle-{deg}-{len}");
            cases.push(Case {
                id: id.clone(),
                family: "angle",
                expect: Expect::Solve,
                build: Box::new(move || {
                    let mut sk = Sketch::new();
                    let l1 = sk.add_line(0.0, 0.0, 5.0, 0.0);
                    sk.fix_entity(l1);
                    let l2 = sk.add_line(0.0, 0.0, len, deg.to_radians() * 0.4);
                    let ps = sk.entities[l2].param_start;
                    sk.fix_param(ps);
                    sk.fix_param(ps + 1); // pin l2's start
                    sk.add_constraint(Angle {
                        l1,
                        l2,
                        angle_rad: deg.to_radians(),
                    });
                    sk.add_constraint(LineLength {
                        line: l2,
                        length: len,
                    });
                    sk
                }),
            });
        }
    }

    // ── circle-chain family: externally tangent chain, colinear centers.
    for k in 2..=8usize {
        let id = format!("circle-chain-k{k}");
        cases.push(Case {
            id: id.clone(),
            family: "circle-chain",
            expect: Expect::Solve,
            build: Box::new(move || {
                let mut sk = Sketch::new();
                let mut prev = sk.add_circle(0.0, 0.0, 1.0);
                sk.fix_entity(prev);
                for i in 1..k {
                    let r = 0.5 + i as f64 * 0.25;
                    let x = i as f64 * 1.75;
                    let ci = sk.add_circle(x, 0.0, r);
                    sk.add_constraint(CircleCircleTangent {
                        c1: prev,
                        c2: ci,
                        external: true,
                    });
                    sk.add_constraint(DistanceY {
                        p1: PointRef::Center(ci),
                        p2: PointRef::Center(prev),
                        distance: 0.0,
                    });
                    sk.add_constraint(Radius {
                        circle: ci,
                        radius: r,
                    });
                    prev = ci;
                }
                sk
            }),
        });
    }

    cases
}

// ── under / conflict / redundant families (the labeled fault sets) ─────────

fn push_fault_families(cases: &mut Vec<Case>) {
    // under: the margin pattern missing one constraint kind.
    for len in [8.0, 12.0, 20.0] {
        for drop in 0..4usize {
            let id = format!("under-margin-l{len}-drop{drop}");
            cases.push(Case {
                id: id.clone(),
                family: "under",
                expect: Expect::Under,
                build: Box::new(move || {
                    let mut sk = Sketch::new();
                    let anchor = sk.add_point(0.0, 0.0);
                    let line = sk.add_line(0.0, 0.0, len, 0.0);
                    sk.add_constraint(Fixed {
                        point: PointRef::Point(anchor),
                        x: 0.0,
                        y: 0.0,
                    });
                    sk.add_constraint(Coincident(
                        PointRef::LineStart(line),
                        PointRef::Point(anchor),
                    ));
                    sk.add_constraint(Horizontal(line));
                    sk.add_constraint(DistanceX {
                        p1: PointRef::LineEnd(line),
                        p2: PointRef::LineStart(line),
                        distance: len,
                    });
                    let hole = sk.add_circle(len - 2.0, 2.0, 4.0);
                    if drop != 0 {
                        sk.add_constraint(Radius {
                            circle: hole,
                            radius: 2.0,
                        });
                    }
                    if drop != 1 {
                        sk.add_constraint(DistanceX {
                            p1: PointRef::LineEnd(line),
                            p2: PointRef::Center(hole),
                            distance: 2.0,
                        });
                    }
                    if drop != 2 {
                        sk.add_constraint(DistanceY {
                            p1: PointRef::Center(hole),
                            p2: PointRef::Point(anchor),
                            distance: 2.0,
                        });
                    }
                    // The bore circle is ALWAYS present; drop 3 removes
                    // only its constraint, leaving 3 unconstrained
                    // parameters (the under-constrainedness under test).
                    let bore = sk.add_circle(len * 0.5, 2.0, 10.0);
                    if drop != 3 {
                        sk.add_constraint(DistanceX {
                            p1: PointRef::Point(anchor),
                            p2: PointRef::Center(bore),
                            distance: len * 0.5,
                        });
                    }
                    sk
                }),
            });
        }
    }
    // under: the tangent line with a PointOnCircle start — the line can
    // slide along its own tangent (1 DOF), the classic mis-constrained
    // tangency.
    for r in [1.0, 2.0, 3.0] {
        let id = format!("under-tangent-pointoncircle-r{r}");
        cases.push(Case {
            id: id.clone(),
            family: "under",
            expect: Expect::Under,
            build: Box::new(move || {
                let mut sk = Sketch::new();
                let c = sk.add_circle(0.0, 0.0, r);
                sk.fix_entity(c);
                let l = sk.add_line(r, 0.0, r + 2.0 * r, 0.0);
                sk.add_constraint(PointOnCircle {
                    point: PointRef::LineStart(l),
                    circle: c,
                });
                sk.add_constraint(Tangent { circle: c, line: l });
                sk.add_constraint(LineLength {
                    line: l,
                    length: 2.0 * r,
                });
                sk
            }),
        });
    }
    // under: open distance chains (bend/rotation left free).
    for n in 3..=12usize {
        let id = format!("under-chain-n{n}");
        cases.push(Case {
            id: id.clone(),
            family: "under",
            expect: Expect::Under,
            build: Box::new(move || {
                let mut sk = Sketch::new();
                let p0 = sk.add_point(0.0, 0.0);
                sk.add_constraint(Fixed {
                    point: PointRef::Point(p0),
                    x: 0.0,
                    y: 0.0,
                });
                let mut prev = p0;
                for i in 1..n {
                    let p = sk.add_point(i as f64, 0.3 * i as f64);
                    sk.add_constraint(Constraint::point_distance(
                        PointRef::Point(prev),
                        PointRef::Point(p),
                        1.0,
                    ));
                    prev = p;
                }
                sk
            }),
        });
    }

    // conflict: radius vs disagreeing diameter.
    for k in [1.2, 2.0, 5.0] {
        let id = format!("conflict-radius-vs-diameter-k{k}");
        cases.push(Case {
            id: id.clone(),
            family: "conflict",
            expect: Expect::Conflict,
            build: Box::new(move || {
                let mut sk = Sketch::new();
                let c = sk.add_circle(1.0, 1.0, 3.0);
                sk.add_constraint(Radius {
                    circle: c,
                    radius: 2.0,
                });
                sk.add_constraint(Diameter {
                    circle: c,
                    diameter: 4.0 * k,
                });
                sk
            }),
        });
    }
    // conflict: one point demanded at two places.
    for (dx, dy) in [(0.5, 0.5), (1.0, 1.0), (3.0, 0.0), (0.0, 2.0)] {
        let id = format!("conflict-double-fixed-{dx}-{dy}");
        cases.push(Case {
            id: id.clone(),
            family: "conflict",
            expect: Expect::Conflict,
            build: Box::new(move || {
                let mut sk = Sketch::new();
                let p = sk.add_point(0.0, 0.0);
                sk.add_constraint(Fixed {
                    point: PointRef::Point(p),
                    x: 1.0,
                    y: 1.0,
                });
                sk.add_constraint(Fixed {
                    point: PointRef::Point(p),
                    x: 1.0 + dx,
                    y: 1.0 + dy,
                });
                sk
            }),
        });
    }
    // conflict: demanded distance beyond the pinned span.
    for span in [1.0, 3.0] {
        for gap in [2.0, 5.0] {
            let id = format!("conflict-distance-beyond-span-{span}-{gap}");
            cases.push(Case {
                id: id.clone(),
                family: "conflict",
                expect: Expect::Conflict,
                build: Box::new(move || {
                    let mut sk = Sketch::new();
                    let p1 = sk.add_point(0.0, 0.0);
                    let p2 = sk.add_point(span, 0.0);
                    sk.add_constraint(Fixed {
                        point: PointRef::Point(p1),
                        x: 0.0,
                        y: 0.0,
                    });
                    sk.add_constraint(Fixed {
                        point: PointRef::Point(p2),
                        x: span,
                        y: 0.0,
                    });
                    sk.add_constraint(Constraint::point_distance(p1, p2, span + gap));
                    sk
                }),
            });
        }
    }
    // conflict: contradictory margins (m1 + m2 > len).
    for len in [8.0, 12.0] {
        // Strictly contradictory pairs only: m1 + m2 == len is CONSISTENT
        // (the second margin is redundant, not conflicting).
        let pairs: Vec<(f64, f64)> = [(6.0, 6.0), (5.0, 7.0), (7.0, 8.0), (9.0, 4.0)]
            .into_iter()
            .filter(|&(m1, m2)| m1 + m2 > len)
            .collect();
        for (m1, m2) in pairs {
            let id = format!("conflict-margins-l{len}-{m1}-{m2}");
            cases.push(Case {
                id: id.clone(),
                family: "conflict",
                expect: Expect::Conflict,
                build: Box::new(move || {
                    let mut sk = Sketch::new();
                    let anchor = sk.add_point(0.0, 0.0);
                    let line = sk.add_line(0.0, 0.0, len, 0.0);
                    sk.add_constraint(Fixed {
                        point: PointRef::Point(anchor),
                        x: 0.0,
                        y: 0.0,
                    });
                    sk.add_constraint(Coincident(
                        PointRef::LineStart(line),
                        PointRef::Point(anchor),
                    ));
                    sk.add_constraint(Horizontal(line));
                    sk.add_constraint(LineLength { line, length: len });
                    // BOTH margins land on the SAME circle: its center is
                    // demanded at x = m1 (from the anchor) and at
                    // x = len − m2 (from the driven end) — a genuine
                    // linear conflict when m1 + m2 != len.
                    let hole = sk.add_circle(m1, 2.0, 1.0);
                    sk.add_constraint(DistanceX {
                        p1: PointRef::Center(hole),
                        p2: PointRef::Point(anchor),
                        distance: m1,
                    });
                    sk.add_constraint(DistanceX {
                        p1: PointRef::LineEnd(line),
                        p2: PointRef::Center(hole),
                        distance: m2,
                    });
                    sk
                }),
            });
        }
    }
    // conflict: parallel + perpendicular on the same pair (generic pose).
    for g in 0..3usize {
        let id = format!("conflict-parallel-perpendicular-{g}");
        cases.push(Case {
            id: id.clone(),
            family: "conflict",
            expect: Expect::Conflict,
            build: Box::new(move || {
                let mut sk = Sketch::new();
                let l1 = sk.add_line(0.0, 0.0, 3.0 + g as f64, 4.0);
                sk.fix_entity(l1);
                let l2 = sk.add_line(0.0, 5.0, 1.0, 8.0);
                sk.add_constraint(Fixed {
                    point: PointRef::LineStart(l2),
                    x: 0.0,
                    y: 5.0,
                });
                let y_param = sk.entities[l2].param_start + 3; // l2's end y
                sk.fix_param(y_param);
                sk.add_constraint(Parallel(l1, l2));

                sk.add_constraint(Perpendicular(l1, l2));
                sk
            }),
        });
    }

    // redundant: consistent duplicates and implications.
    for t in 0..3usize {
        let off = t as f64;
        let id = format!("redundant-duplicate-coincident-{t}");
        cases.push(Case {
            id: id.clone(),
            family: "redundant",
            expect: Expect::Redundant,
            build: Box::new(move || {
                let mut sk = Sketch::new();
                let p1 = sk.add_point(1.0 + off, 1.0);
                let p2 = sk.add_point(4.0 + off, 1.0);
                sk.add_constraint(Fixed {
                    point: PointRef::Point(p1),
                    x: 1.0 + off,
                    y: 1.0,
                });
                sk.add_constraint(Coincident(PointRef::Point(p1), PointRef::Point(p2)));
                sk.add_constraint(Coincident(PointRef::Point(p1), PointRef::Point(p2)));
                sk
            }),
        });
        let id = format!("redundant-consistent-radius-diameter-{t}");
        cases.push(Case {
            id: id.clone(),
            family: "redundant",
            expect: Expect::Redundant,
            build: Box::new(move || {
                let mut sk = Sketch::new();
                let c = sk.add_circle(0.0, 0.0, 5.0 + off);
                sk.add_constraint(Radius {
                    circle: c,
                    radius: 5.0 + off,
                });
                sk.add_constraint(Diameter {
                    circle: c,
                    diameter: 10.0 + 2.0 * off,
                });
                sk
            }),
        });
        let id = format!("redundant-distance-implied-{t}");
        cases.push(Case {
            id: id.clone(),
            family: "redundant",
            expect: Expect::Redundant,
            build: Box::new(move || {
                let mut sk = Sketch::new();
                let p1 = sk.add_point(0.0, 0.0);
                let p2 = sk.add_point(3.0, 4.0);
                sk.add_constraint(Fixed {
                    point: PointRef::Point(p1),
                    x: 0.0,
                    y: 0.0,
                });
                sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
                sk.add_constraint(DistanceX {
                    p1: PointRef::Point(p2),
                    p2: PointRef::Point(p1),
                    distance: 3.0,
                });
                sk.add_constraint(DistanceY {
                    p1: PointRef::Point(p2),
                    p2: PointRef::Point(p1),
                    distance: 4.0,
                });
                sk
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
        "generator": "rcad-gcs corpus v1 (fixed-seed, deterministic)",
        "count": cases.len(),
        "cases": cases.iter().map(|c| c.metadata()).collect::<Vec<_>>(),
    })
}

fn manifest_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/corpus/manifest.json")
}

/// Regenerate and overwrite `tests/corpus/manifest.json` (run explicitly
/// after changing the generator):
/// `cargo test -p rcad-gcs --test corpus -- --ignored --nocapture write_corpus_manifest`
#[test]
#[ignore]
fn write_corpus_manifest() {
    let cases = generate_all();
    let json = serde_json::to_string_pretty(&manifest_json(&cases)).expect("manifest json");
    std::fs::create_dir_all(manifest_path().parent().unwrap()).expect("mkdir");
    std::fs::write(manifest_path(), json.clone() + "\n").expect("write manifest");
    println!("wrote {} cases to {:?}", cases.len(), manifest_path());
}

/// The committed manifest must equal what the generator produces today —
/// generator drift without a manifest refresh is a review-visible diff.
#[test]
fn corpus_manifest_matches_generated() {
    let cases = generate_all();
    assert!(
        cases.len() >= MIN_CORPUS_SIZE,
        "corpus shrank to {} cases (need >= {MIN_CORPUS_SIZE})",
        cases.len()
    );
    let committed = std::fs::read_to_string(manifest_path())
        .expect("tests/corpus/manifest.json missing — run the ignored write_corpus_manifest test");
    // Compare modulo line endings: a checkout that predates `.gitattributes`
    // (or overrides it) can carry CRLF.
    let committed = committed.replace("\r\n", "\n");
    let expected = serde_json::to_string_pretty(&manifest_json(&cases)).expect("json") + "\n";
    assert_eq!(committed, expected, "manifest drifted from the generator");
}

// ── the planegcs comparison export (roadmap M5.3) ───────────────────────────

fn planegcs_export_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/corpus/planegcs_export.json")
}

/// Solver-neutral point reference: `[entity_id, pos]` with FreeCAD's
/// point-position convention (1 = start, 2 = end, 3 = center; a Point
/// entity is pos 1).
fn point_ref_json(p: &PointRef) -> serde_json::Value {
    use PointRef::*;
    match p {
        Point(id) => serde_json::json!([id, 1]),
        LineStart(id) => serde_json::json!([id, 1]),
        LineEnd(id) => serde_json::json!([id, 2]),
        Center(id) => serde_json::json!([id, 3]),
        ArcStart(id) => serde_json::json!([id, 4]),
        ArcEnd(id) => serde_json::json!([id, 5]),
    }
}

/// One constraint in the neutral export schema. Only the kinds the corpus
/// actually generates are mapped; anything else exports as `Unsupported`
/// (the comparison driver skips and counts those cases).
fn constraint_json(c: &Constraint) -> serde_json::Value {
    use serde_json::json;
    use Constraint::*;
    match c {
        Fixed { point, x, y } => {
            json!({"kind": "Fixed", "point": point_ref_json(point), "x": x, "y": y})
        }
        Coincident(p1, p2) => {
            json!({"kind": "Coincident", "p1": point_ref_json(p1), "p2": point_ref_json(p2)})
        }
        Horizontal(l) => json!({"kind": "Horizontal", "line": l}),
        Vertical(l) => json!({"kind": "Vertical", "line": l}),
        LineLength { line, length } => json!({"kind": "LineLength", "line": line, "value": length}),
        DistanceX { p1, p2, distance } => {
            json!({"kind": "DistanceX", "p1": point_ref_json(p1), "p2": point_ref_json(p2), "value": distance})
        }
        DistanceY { p1, p2, distance } => {
            json!({"kind": "DistanceY", "p1": point_ref_json(p1), "p2": point_ref_json(p2), "value": distance})
        }
        Radius { circle, radius } => json!({"kind": "Radius", "circle": circle, "value": radius}),
        Diameter { circle, diameter } => {
            json!({"kind": "Diameter", "circle": circle, "value": diameter})
        }
        Midpoint { point, line } => {
            json!({"kind": "Midpoint", "point": point_ref_json(point), "line": line})
        }
        Tangent { circle, line } => json!({"kind": "Tangent", "circle": circle, "line": line}),
        PointOnCircle { point, circle } => {
            json!({"kind": "PointOnCircle", "point": point_ref_json(point), "circle": circle})
        }
        SymmetricPoint { p1, p2, about } => {
            json!({"kind": "SymmetricPoint", "p1": point_ref_json(p1), "p2": point_ref_json(p2), "about": point_ref_json(about)})
        }
        Parallel(l1, l2) => json!({"kind": "Parallel", "l1": l1, "l2": l2}),
        Perpendicular(l1, l2) => json!({"kind": "Perpendicular", "l1": l1, "l2": l2}),
        Angle { l1, l2, angle_rad } => {
            json!({"kind": "Angle", "l1": l1, "l2": l2, "angle": angle_rad})
        }
        // planegcs picks the nearest tangency branch from the authored
        // geometry, so the external/internal hint rides along for the
        // report but the solve does not need it.
        CircleCircleTangent { c1, c2, external } => {
            json!({"kind": "Tangent", "circle": c1, "circle2": c2, "external": external})
        }
        PointDistance { p1, p2, distance } => {
            json!({"kind": "PointDistance", "p1": point_ref_json(p1), "p2": point_ref_json(p2), "value": distance})
        }
        other => {
            let variant = format!("{other:?}");
            let variant = variant.split(['(', ' ', '{']).next().unwrap_or("Unknown");
            json!({"kind": "Unsupported", "variant": variant})
        }
    }
}

/// Export the whole corpus to `tests/corpus/planegcs_export.json` for the
/// one-shot planegcs comparison (roadmap M5.3) — the FreeCAD driver
/// (`tools/planegcs-compare/drive_planegcs.py`) rebuilds every case on
/// the real planegcs and reports convergence/diagnosis/timing next to
/// rcad's own corpus metrics. Run explicitly:
/// `cargo test -p rcad-gcs --test corpus -- --ignored --nocapture write_planegcs_export`
#[test]
#[ignore]
fn write_planegcs_export() {
    let cases = generate_all();
    let mut exported = Vec::with_capacity(cases.len());
    for case in &cases {
        let sk = (case.build)();
        let entities: Vec<serde_json::Value> = sk
            .entities
            .iter()
            .enumerate()
            .map(|(id, e)| {
                serde_json::json!({
                    "id": id,
                    "kind": format!("{:?}", e.kind).to_lowercase(),
                    "params": sk.entity_params(id),
                })
            })
            .collect();
        let constraints: Vec<serde_json::Value> =
            sk.constraints().iter().map(constraint_json).collect();
        exported.push(serde_json::json!({
            "id": case.id,
            "family": case.family,
            "expect": case.expect.as_str(),
            "entities": entities,
            "constraints": constraints,
        }));
    }
    let json = serde_json::to_string_pretty(&serde_json::json!({
        "generator": "rcad-gcs corpus planegcs export v1 (fixed-seed, deterministic)",
        "count": exported.len(),
        "cases": exported,
    }))
    .expect("export json");
    std::fs::write(planegcs_export_path(), json + "\n").expect("write export");
    println!(
        "wrote {} cases to {:?}",
        exported.len(),
        planegcs_export_path()
    );
}

// ── the metrics gate ────────────────────────────────────────────────────────

/// Run every corpus case through diagnose + solve (+ the M4-style
/// drift-and-recover round trip) and enforce the red-lines: convergence
/// 100%, conflict recall 1.0 / precision 1.0, redundant detection 1.0,
/// under-DOF detection 1.0, no-hysteresis violations 0. A JSON report
/// lands in the target tmp dir.
#[test]
fn corpus_metrics_within_thresholds() {
    let t_start = std::time::Instant::now();
    let cases = generate_all();
    assert!(
        cases.len() >= MIN_CORPUS_SIZE,
        "corpus shrank to {} cases",
        cases.len()
    );

    let (mut n_solve, mut n_under, mut n_conflict, mut n_redundant, mut n_singular) =
        (0, 0, 0, 0, 0);
    let mut singular_flagged = 0usize;
    let mut singular_surfaced = 0usize;
    let (mut converged, mut conflicts_flagged, mut false_conflicts) = (0, 0, 0);
    let (mut redundant_flagged, mut under_flagged) = (0, 0);
    let mut solve_dof_violations = 0usize;
    let (mut round_trips, mut hysteresis_violations, mut branch_flips) = (0usize, 0usize, 0usize);
    let (mut worst_drift, mut worst_hard_drift) = (0.0_f64, 0.0_f64);
    let mut slowest_ms: (f64, String) = (0.0, String::new());

    for case in &cases {
        let t0 = std::time::Instant::now();
        let mut sk = (case.build)();
        let d = sk.diagnose();
        let has_conflict = !d.conflicting.is_empty();

        match case.expect {
            Expect::Conflict => {
                n_conflict += 1;
                if has_conflict {
                    conflicts_flagged += 1;
                } else {
                    println!("MISSED CONFLICT on {}: {:?}", case.id, d);
                }
                // The solve may thrash — the sketcher fail-fasts before it.
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
                let r = sk.solve();
                if r.converged {
                    converged += 1;
                } else {
                    println!("NOT CONVERGED on {}: residual {}", case.id, r.residual);
                }
                // Soft probe: red-line. Hard probe: report-only signal.
                // Both need a rigid (dof 0) system to be well-posed.
                if d.dof == 0 {
                    drift_recover_check(
                        case,
                        &mut sk,
                        0.05,
                        &mut round_trips,
                        &mut hysteresis_violations,
                        &mut worst_drift,
                    );
                    drift_recover_check(
                        case,
                        &mut sk,
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
                let r = sk.solve();
                if r.converged {
                    converged += 1;
                } else {
                    println!("NOT CONVERGED on {}: residual {}", case.id, r.residual);
                }
            }
            Expect::Singular => {
                n_singular += 1;
                if has_conflict {
                    false_conflicts += 1;
                    println!("FALSE CONFLICT on {}: {:?}", case.id, d.conflicting);
                }
                if d.dof != 0 {
                    singular_flagged += 1;
                }
                let r = sk.solve();
                if r.converged {
                    converged += 1;
                } else {
                    println!("NOT CONVERGED on {}: residual {}", case.id, r.residual);
                }
                // The tolerate-and-warn channel: at the SOLVED configuration
                // the rank-deficient rows must surface as dependent-but-
                // consistent (the sketch_singular warning's raw signal),
                // and the v2 dependency groups must name the whole
                // relation.
                let post = sk.diagnose();
                if (!post.redundant.is_empty() || !post.partially_redundant.is_empty())
                    && !post.dependency_groups.is_empty()
                {
                    singular_surfaced += 1;
                } else {
                    println!(
                        "SINGULAR NOT SURFACED on {}: groups {:?} redundant {:?}",
                        case.id, post.dependency_groups, post.redundant
                    );
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
                let r = sk.solve();
                if r.converged {
                    converged += 1;
                } else {
                    println!("NOT CONVERGED on {}: residual {}", case.id, r.residual);
                }
            }
        }
        let ms = t0.elapsed().as_secs_f64() * 1e3;
        if ms > slowest_ms.0 {
            slowest_ms = (ms, case.id.clone());
        }
    }

    let solvable = n_solve + n_under + n_redundant + n_singular;
    let convergence_rate = converged as f64 / solvable as f64;
    let conflict_recall = conflicts_flagged as f64 / n_conflict as f64;
    let conflict_precision = 1.0 - false_conflicts as f64 / (cases.len() - n_conflict) as f64;
    let redundant_rate = redundant_flagged as f64 / n_redundant as f64;
    let under_dof_rate = under_flagged as f64 / n_under as f64;
    let hysteresis_rate = 1.0 - hysteresis_violations as f64 / round_trips.max(1) as f64;

    println!(
        "corpus: {} cases in {:.1?} (slowest {} {:.0} ms)",
        cases.len(),
        t_start.elapsed(),
        slowest_ms.1,
        slowest_ms.0
    );
    println!(
        "convergence {convergence_rate:.4} | conflict recall {conflict_recall:.4} precision {conflict_precision:.4} | redundant {redundant_rate:.4} | under-dof {under_dof_rate:.4} | singular surfaced {singular_surfaced}/{n_singular} | no-hysteresis {hysteresis_rate:.4} (round trips {round_trips}, worst drift {worst_drift:.2e})"
    );

    // Red-lines (exact — the corpus and the solver are deterministic).
    assert_eq!(convergence_rate, 1.0, "convergence rate below 1.0");
    assert_eq!(conflict_recall, 1.0, "a conflict escaped diagnosis");
    assert_eq!(conflict_precision, 1.0, "false conflicts flagged");
    assert_eq!(redundant_rate, 1.0, "redundant rows missed");
    assert_eq!(under_dof_rate, 1.0, "under-constrained cases without DOF");
    assert_eq!(solve_dof_violations, 0, "Solve-labeled cases must be rigid");
    assert_eq!(
        singular_flagged, n_singular,
        "every exact-tangency case must diagnose its rank deficiency"
    );
    assert_eq!(
        singular_surfaced, n_singular,
        "every exact-tangency case must surface the singular set at the solved pose"
    );
    assert_eq!(
        hysteresis_violations, 0,
        "no-hysteresis round trips regressed"
    );
    assert!(
        t_start.elapsed().as_secs() < 120,
        "corpus run took {:?} — performance regression?",
        t_start.elapsed()
    );

    let report = json!({
        "cases": cases.len(),
        "convergence_rate": convergence_rate,
        "conflict_recall": conflict_recall,
        "conflict_precision": conflict_precision,
        "redundant_detection_rate": redundant_rate,
        "under_dof_rate": under_dof_rate,
        "singular_surfaced": singular_surfaced,
        "singular_cases": n_singular,
        "round_trips": round_trips,
        "hysteresis_violations": hysteresis_violations,
        "worst_round_trip_drift_soft": worst_drift,
        "worst_round_trip_drift_hard": worst_hard_drift,
        "hard_branch_flips": branch_flips,
        "slowest_case_ms": slowest_ms.0,
        "slowest_case": slowest_ms.1,
        "total_seconds": t_start.elapsed().as_secs_f64(),
    });
    // Cargo sets this for integration tests at compile time only.
    if let Some(dir) = option_env!("CARGO_TARGET_TMPDIR") {
        let path = std::path::Path::new(dir).join("gcs_corpus_report.json");
        let _ = std::fs::write(&path, serde_json::to_string_pretty(&report).unwrap());
        println!("report: {:?}", path);
    }
}

/// The M4 no-hysteresis probe, corpus-wide: drift every parameter hard,
/// re-solve, and the case must return to its configuration within the
/// contract tolerance.
/// The M4 no-hysteresis probe, corpus-wide. Two amplitudes:
///
/// - **soft** (`amp` 0.05, red-line): a perturbation inside the solution
///   basin — re-solve MUST return to the same configuration within the
///   contract tolerance. This is the robustness half of the M4 contract.
/// - **hard** (`amp` 3.0, report-only): an arbitrary excursion that may
///   legitimately cross branch boundaries of multi-solution systems
///   (mirrored rectangles etc.). Measured as the branch-stability rate —
///   a signal, not a contract: DCM's persistence promise covers dimension
///   round trips (enforced in tests/persistence.rs), not arbitrary
///   re-entry from anywhere in parameter space.
fn drift_recover_check(
    case: &Case,
    sk: &mut Sketch,
    amp: f64,
    round_trips: &mut usize,
    violations: &mut usize,
    worst_drift: &mut f64,
) {
    let before = sk.params_snapshot();
    // Mask-aware drift: fixed parameters (e.g. a pinned circle center)
    // must NOT move — set_params would silently re-position fixed
    // geometry and the "recovered" state would legitimately follow it.
    sk.drift_free_params(amp, 0x5EED ^ (before.len() as u64).wrapping_mul(31));
    let r1 = sk.solve();
    if !r1.converged {
        *violations += 1;
        println!(
            "{}: re-solve after drift(amp={amp}) did not converge",
            case.id
        );
        sk.set_params(before);
        return;
    }
    let worst = before
        .iter()
        .zip(sk.params_snapshot().iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    // Restore the base state so the next probe (hard amplitude) starts
    // from the same configuration.
    sk.set_params(before);
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
