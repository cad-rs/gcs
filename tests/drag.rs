//! Scripted drags (P5.1): the cursor as a soft target.
//!
//! Six drags over four fixtures, each a closed path so that the round trip
//! can be looked at: a free rectangle's corner (two degrees of freedom), a
//! four-bar linkage's coupler point (one), a chain of tangent circles
//! (four), and forty free rectangles in one component (200 entities). What
//! each drag must deliver:
//!
//! - the hard constraints stay satisfied on every frame;
//! - a handle with the freedom to reach its cursor reaches it, and the
//!   frame says so;
//! - a cursor outside what the constraints can reach leaves the geometry
//!   at the nearest configuration the frame found, says `reached` false,
//!   and still satisfies the constraints;
//! - a closed drag comes back to the configuration it started from (no
//!   hysteresis), and the chirality of the fixture does not flip along the
//!   way.
//!
//! The one exception is the rectangle dragged through zero width, where
//! the cursor itself asks for the mirrored configuration: the flip is
//! recorded, not forbidden. Guarding it is P5.3's job (`BranchFlip`).
//!
//! `report_drag_frames` (ignored by default) prints the same drags under
//! the two things a host could do before this API existed, which is where
//! the numbers in the plan's P5 baseline come from:
//!
//! ```text
//! cargo test --release --all-features --test drag -- --ignored --nocapture
//! ```

use rcad_gcs::constraint::Constraint;
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::{EntityId, PointRef};
use rcad_gcs::{Sketch, SolveStatus};

/// A scripted drag: the sketch, the point under the cursor, the triples
/// whose orientation must not flip, and the cursor path.
struct DragCase {
    id: &'static str,
    sketch: Sketch,
    handle: PointRef,
    triples: Vec<(PointRef, PointRef, PointRef)>,
    /// Cursor positions, frame by frame; the last one is the first.
    path: Vec<(f64, f64)>,
    /// Whether every cursor on the path is reachable.
    reachable: bool,
    /// Whether the handle's own position determines the rest of the
    /// configuration: with more free motions than the cursor pins, where
    /// the interior ends up after a closed path is a property of the
    /// path, not of the solver (a chain dragged round a loop and back
    /// does not re-straighten itself), so only these cases owe the round
    /// trip.
    determines_the_rest: bool,
}

/// The sign pattern of the orientation of each watched triple: a
/// continuous drag inside the feasible region never turns one over.
fn chirality(sk: &Sketch, triples: &[(PointRef, PointRef, PointRef)]) -> Vec<i8> {
    triples
        .iter()
        .map(|(a, b, c)| {
            let (a, b, c) = (
                sk.point_coords(*a),
                sk.point_coords(*b),
                sk.point_coords(*c),
            );
            let cross = (b - a).perp_dot(c - a);
            if cross > 1e-9 {
                1
            } else if cross < -1e-9 {
                -1
            } else {
                0
            }
        })
        .collect()
}

/// A closed loop of `frames` cursor positions, starting where the handle
/// is: around a centre `radius` to its left.
fn loop_path(x: f64, y: f64, radius: f64, frames: usize) -> Vec<(f64, f64)> {
    (0..=frames)
        .map(|i| {
            let t = i as f64 / frames as f64 * std::f64::consts::TAU;
            (x + radius * t.cos() - radius, y + radius * t.sin())
        })
        .collect()
}

/// Out and back along a straight line, `frames` positions each way.
fn sweep_path(from: (f64, f64), to: (f64, f64), frames: usize) -> Vec<(f64, f64)> {
    let step = |i: usize| {
        let t = i as f64 / frames as f64;
        (from.0 + (to.0 - from.0) * t, from.1 + (to.1 - from.1) * t)
    };
    (0..=frames)
        .map(step)
        .chain((0..frames).rev().map(step))
        .collect()
}

/// A rectangle with horizontal and vertical sides and no dimensions: width
/// and height are free, so a corner can be dragged. Returns the sketch,
/// its origin point and its four lines (bottom, right, top, left).
fn free_rectangle(w: f64, h: f64) -> (Sketch, EntityId, [EntityId; 4]) {
    let mut sk = Sketch::new();
    let origin = sk.add_point(0.0, 0.0);
    let bottom = sk.add_line(0.0, 0.0, w, 0.0);
    let right = sk.add_line(w, 0.0, w, h);
    let top = sk.add_line(w, h, 0.0, h);
    let left = sk.add_line(0.0, h, 0.0, 0.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Point(origin), 0.0, 0.0));
    sk.add_constraint(Coincident(
        PointRef::LineStart(bottom),
        PointRef::Point(origin),
    ));
    for (end, start) in [(bottom, right), (right, top), (top, left), (left, bottom)] {
        sk.add_constraint(Coincident(
            PointRef::LineEnd(end),
            PointRef::LineStart(start),
        ));
    }
    for line in [bottom, top] {
        sk.add_constraint(Horizontal(line));
    }
    for line in [right, left] {
        sk.add_constraint(Vertical(line));
    }
    (sk, origin, [bottom, right, top, left])
}

/// A four-bar linkage: ground pivots at `(0,0)` and `(span,0)`, a crank, a
/// coupler and a rocker of the given lengths. One degree of freedom, so a
/// cursor anywhere off its curve is unreachable — the case a hard target
/// cannot drag at all.
fn four_bar(span: f64, crank: f64, coupler: f64, rocker: f64) -> (Sketch, [EntityId; 4]) {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(0.0, crank);
    let c = sk.add_point(span - 1.0, rocker * 0.9);
    let d = sk.add_point(span, 0.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Point(a), 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(PointRef::Point(d), span, 0.0));
    for (p, q, length) in [(a, b, crank), (b, c, coupler), (c, d, rocker)] {
        sk.add_constraint(Constraint::point_distance(
            PointRef::Point(p),
            PointRef::Point(q),
            length,
        ));
    }
    (sk, [a, b, c, d])
}

/// `count` circles of one radius in a tangent chain, the first pinned.
///
/// The chain is authored turning, a sixth of a turn at each circle, not
/// straight: a straight chain is at the end of its reach, where the fold
/// it could make is a singular motion and a cursor just outside it is
/// unreachable. Four free degrees of freedom.
fn circle_chain(count: usize, radius: f64) -> (Sketch, Vec<EntityId>) {
    let mut sk = Sketch::new();
    let mut at = (0.0, 0.0);
    let mut ids: Vec<EntityId> = Vec::new();
    for i in 0..count {
        ids.push(sk.add_circle(at.0, at.1, radius));
        let heading = i as f64 * std::f64::consts::FRAC_PI_6;
        at = (
            at.0 + 2.0 * radius * heading.cos(),
            at.1 + 2.0 * radius * heading.sin(),
        );
    }
    sk.add_constraint(Constraint::fix_point(PointRef::Center(ids[0]), 0.0, 0.0));
    for &id in &ids {
        sk.add_constraint(Radius { circle: id, radius });
    }
    for pair in ids.windows(2) {
        sk.add_constraint(CircleCircleTangent {
            c1: pair[0],
            c2: pair[1],
            external: true,
        });
    }
    (sk, ids)
}

/// `count` free rectangles, each hanging off the previous one's corner:
/// one component, five entities each, nothing dimensioned.
fn free_rectangle_chain(count: usize, w: f64, h: f64) -> (Sketch, EntityId) {
    let (mut sk, _, lines) = free_rectangle(w, h);
    let mut previous = lines[1];
    for i in 1..count {
        let x = i as f64 * w;
        let bottom = sk.add_line(x, 0.0, x + w, 0.0);
        let right = sk.add_line(x + w, 0.0, x + w, h);
        let top = sk.add_line(x + w, h, x, h);
        let left = sk.add_line(x, h, x, 0.0);
        sk.add_constraint(Coincident(
            PointRef::LineStart(bottom),
            PointRef::LineStart(previous),
        ));
        for (end, start) in [(bottom, right), (right, top), (top, left), (left, bottom)] {
            sk.add_constraint(Coincident(
                PointRef::LineEnd(end),
                PointRef::LineStart(start),
            ));
        }
        for line in [bottom, top] {
            sk.add_constraint(Horizontal(line));
        }
        for line in [right, left] {
            sk.add_constraint(Vertical(line));
        }
        previous = right;
    }
    (sk, previous)
}

fn drag_cases() -> Vec<DragCase> {
    let mut cases = Vec::new();

    // 1. A rectangle's corner round a loop that stays well inside the
    //    feasible region.
    let (sk, origin, lines) = free_rectangle(4.0, 3.0);
    let corner_triples = vec![(
        PointRef::Point(origin),
        PointRef::LineEnd(lines[0]),
        PointRef::LineEnd(lines[1]),
    )];
    cases.push(DragCase {
        id: "rect-corner-loop",
        sketch: sk,
        handle: PointRef::LineEnd(lines[1]),
        triples: corner_triples.clone(),
        path: loop_path(4.0, 3.0, 1.5, 64),
        reachable: true,
        determines_the_rest: true,
    });

    // 2. The same corner dragged through zero width and back: the cursor
    //    asks for the mirrored configuration itself.
    let (sk, _, lines) = free_rectangle(4.0, 3.0);
    cases.push(DragCase {
        id: "rect-corner-through-zero",
        sketch: sk,
        handle: PointRef::LineEnd(lines[1]),
        triples: corner_triples,
        path: sweep_path((4.0, 3.0), (-3.0, 3.0), 48),
        reachable: true,
        determines_the_rest: true,
    });

    // 3. A four-bar's coupler point swung along its own curve: every
    //    cursor but the endpoints is off the curve, so the handle tracks
    //    only as far as one degree of freedom allows.
    let (mut sk, p) = four_bar(10.0, 4.0, 7.0, 5.0);
    assert!(sk.solve().converged, "four-bar fixture");
    let c = sk.point_coords(PointRef::Point(p[2]));
    let linkage_triples = vec![
        (
            PointRef::Point(p[0]),
            PointRef::Point(p[1]),
            PointRef::Point(p[2]),
        ),
        (
            PointRef::Point(p[1]),
            PointRef::Point(p[2]),
            PointRef::Point(p[3]),
        ),
    ];
    cases.push(DragCase {
        id: "four-bar-arc",
        sketch: sk,
        handle: PointRef::Point(p[2]),
        triples: linkage_triples.clone(),
        path: sweep_path((c.x, c.y), (c.x + 1.0, c.y - 3.5), 48),
        reachable: false,
        determines_the_rest: true,
    });

    // 4. The same coupler point pulled far outside the reachable region.
    let (mut sk, p) = four_bar(10.0, 4.0, 7.0, 5.0);
    assert!(sk.solve().converged, "four-bar fixture");
    let c = sk.point_coords(PointRef::Point(p[2]));
    cases.push(DragCase {
        id: "four-bar-unreachable",
        sketch: sk,
        handle: PointRef::Point(p[2]),
        triples: linkage_triples,
        path: sweep_path((c.x, c.y), (c.x + 14.0, c.y + 6.0), 48),
        reachable: false,
        determines_the_rest: true,
    });

    // 5. The last centre of a tangent circle chain round a loop. The
    //    chain has four degrees of freedom, so the loop is reachable.
    let (mut sk, ids) = circle_chain(5, 1.0);
    assert!(sk.solve().converged, "circle-chain fixture");
    let last = *ids.last().expect("a chain of five");
    let centre = sk.point_coords(PointRef::Center(last));
    // Nothing to watch for a flip here: three consecutive centres going
    // collinear is an ordinary configuration of a chain held only by
    // distances — the Jacobian is not singular there and the chain bends
    // the other way continuously — so the orientation of such a triple is
    // no branch signature. What a drag of the chain owes is the round
    // trip, which this case checks.
    cases.push(DragCase {
        id: "circle-chain-loop",
        sketch: sk,
        handle: PointRef::Center(last),
        triples: Vec::new(),
        path: loop_path(centre.x, centre.y, 1.2, 64),
        reachable: true,
        determines_the_rest: false,
    });

    // 6. The frame budget: 40 free rectangles (200 entities) in one
    //    component, the far corner round a small loop.
    let (mut sk, right) = free_rectangle_chain(40, 4.0, 3.0);
    assert!(sk.solve().converged, "200-entity fixture");
    let corner = sk.point_coords(PointRef::LineEnd(right));
    cases.push(DragCase {
        id: "chain-200-entities",
        sketch: sk,
        handle: PointRef::LineEnd(right),
        triples: Vec::new(),
        path: loop_path(corner.x, corner.y, 1.0, 32),
        reachable: true,
        determines_the_rest: false,
    });

    cases
}

/// What one scripted drag did.
struct Run {
    frames: usize,
    unconverged: usize,
    /// Frames that reported the cursor as reached.
    reached: usize,
    /// Genuine +1 ↔ −1 turn-overs of a watched triple (a pass through
    /// collinear is not one unless it comes out the other side).
    flips: usize,
    /// The largest residual any frame left behind.
    worst_residual: f64,
    /// The largest tracking error of a frame whose cursor was reachable.
    worst_tracking: f64,
    /// The largest parameter difference across the closed path.
    round_trip: f64,
    /// Whether the chirality came back to what it was.
    returned: bool,
    times_ms: Vec<f64>,
}

/// Run `case` through the drag API, frame by frame.
fn run_drag(case: &DragCase) -> Run {
    let mut sk = case.sketch.clone();
    let start = sk.params_snapshot();
    let before = chirality(&sk, &case.triples);
    let mut signature = before.clone();
    let mut run = Run {
        frames: 0,
        unconverged: 0,
        reached: 0,
        flips: 0,
        worst_residual: 0.0,
        worst_tracking: 0.0,
        round_trip: 0.0,
        returned: true,
        times_ms: Vec::new(),
    };

    let drag = sk.begin_drag(&[case.handle]);
    for &cursor in &case.path {
        let t0 = std::time::Instant::now();
        let frame = sk.drag_to(&drag, &[cursor]);
        run.times_ms.push(t0.elapsed().as_secs_f64() * 1e3);
        run.frames += 1;
        if frame.status > SolveStatus::ConvergedToPrecision {
            run.unconverged += 1;
            println!("{}: {:?} at {cursor:?}", case.id, frame.status);
        }
        if frame.reached {
            run.reached += 1;
        }
        run.worst_residual = run.worst_residual.max(frame.residual);
        if case.reachable {
            run.worst_tracking = run.worst_tracking.max(frame.tracking_error);
        }
        let now = chirality(&sk, &case.triples);
        if signature
            .iter()
            .zip(&now)
            .any(|(was, is)| *was != 0 && *is != 0 && was != is)
        {
            run.flips += 1;
        }
        for (slot, s) in signature.iter_mut().zip(&now) {
            if *s != 0 {
                *slot = *s;
            }
        }
    }
    let end = sk.end_drag(drag);
    assert!(
        end.converged,
        "{}: the final solve left {:.3e}",
        case.id, end.residual
    );
    run.round_trip = start
        .iter()
        .zip(sk.params_snapshot().iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    run.returned = chirality(&sk, &case.triples) == before;
    run
}

fn percentile(sorted: &[f64], q: f64) -> f64 {
    let i = ((sorted.len() - 1) as f64 * q).round() as usize;
    sorted[i]
}

/// The hard constraints hold on every frame of every scripted drag, and a
/// reachable cursor is reached.
#[test]
fn every_frame_keeps_the_constraints() {
    for case in &drag_cases() {
        let run = run_drag(case);
        assert_eq!(
            run.unconverged, 0,
            "{}: {} frames of {} did not converge",
            case.id, run.unconverged, run.frames
        );
        assert!(
            run.worst_residual <= 1e-9,
            "{}: worst frame residual {:.3e}",
            case.id,
            run.worst_residual
        );
        if case.reachable {
            assert_eq!(
                run.reached, run.frames,
                "{}: {} frames of {} reached the cursor",
                case.id, run.reached, run.frames
            );
            assert!(
                run.worst_tracking <= 1e-7,
                "{}: worst tracking error {:.3e}",
                case.id,
                run.worst_tracking
            );
        }
    }
}

/// A closed drag comes back to the configuration it started from, and the
/// fixture's chirality with it.
///
/// Two of the drags still flip, and are asserted to: the rectangle dragged
/// through zero width, where the cursor asks for the mirrored
/// configuration itself, and the four-bar pulled right out of its reach,
/// where the walk arrives at the configuration where the linkage is
/// straight and crosses it. Guarding both is P5.3's job (a chirality
/// signature, a shortened step, `BranchFlip`); this test is the record of
/// where it has to arrive, so it will fail — here — when the guard lands.
#[test]
fn a_closed_drag_comes_back() {
    for case in &drag_cases() {
        let run = run_drag(case);
        if matches!(case.id, "rect-corner-through-zero" | "four-bar-unreachable") {
            assert!(
                run.flips > 0,
                "{}: used to flip without a chirality guard",
                case.id
            );
            continue;
        }
        assert_eq!(run.flips, 0, "{}: {} chirality flips", case.id, run.flips);
        assert!(run.returned, "{}: the chirality did not come back", case.id);
        if !case.determines_the_rest {
            // The interior is free to end up elsewhere, but not anywhere:
            // a minimum-move step never takes it further than the cursor
            // itself travelled.
            let travel = case
                .path
                .iter()
                .map(|&(x, y)| (x - case.path[0].0).hypot(y - case.path[0].1))
                .fold(0.0_f64, f64::max);
            assert!(
                run.round_trip <= 2.0 * travel,
                "{}: the interior drifted {:.3e} over a cursor travel of {travel:.3}",
                case.id,
                run.round_trip
            );
            continue;
        }
        if case.reachable {
            assert!(
                run.round_trip <= 1e-7,
                "{}: round-trip drift {:.3e}",
                case.id,
                run.round_trip
            );
        }
    }
}

/// A cursor the constraints cannot reach: the frame says so, the geometry
/// stays feasible, and the handle ends up nearer the cursor than it began
/// (it walked as far as it could).
#[test]
fn an_unreachable_cursor_walks_as_far_as_it_can() {
    let cases = drag_cases();
    let case = cases
        .iter()
        .find(|c| c.id == "four-bar-unreachable")
        .expect("the unreachable four-bar");
    let mut sk = case.sketch.clone();
    let drag = sk.begin_drag(&[case.handle]);
    let far = *case.path.iter().max_by(|a, b| a.0.total_cmp(&b.0)).unwrap();
    let was =
        (sk.point_coords(case.handle).x - far.0).hypot(sk.point_coords(case.handle).y - far.1);
    let frame = sk.drag_to(&drag, &[far]);
    assert!(
        !frame.reached,
        "a cursor off the one-degree-of-freedom curve was reported as reached"
    );
    assert!(
        frame.status <= SolveStatus::ConvergedToPrecision,
        "the frame failed: {:?}",
        frame.status
    );
    assert!(
        frame.residual <= 1e-9,
        "the constraints were left at {:.3e}",
        frame.residual
    );
    assert!(
        frame.tracking_error < was,
        "the handle did not get nearer: {:.3} from {was:.3}",
        frame.tracking_error
    );
    let result = sk.end_drag(drag);
    assert!(
        result.converged,
        "the final solve left {:.3e}",
        result.residual
    );
}

/// A point no constraint touches follows its cursor exactly; a point whose
/// parameters are held does not move at all.
#[test]
fn a_loose_point_follows_and_a_held_one_does_not() {
    let mut sk = Sketch::new();
    let loose = sk.add_point(1.0, 1.0);
    let held = sk.add_point(5.0, 5.0);
    sk.fix_entity(held);
    let drag = sk.begin_drag(&[sk.point(loose), sk.point(held)]);
    let frame = sk.drag_to(&drag, &[(3.0, -2.0), (0.0, 0.0)]);
    assert_eq!(sk.point_coords(sk.point(loose)).to_array(), [3.0, -2.0]);
    assert_eq!(sk.point_coords(sk.point(held)).to_array(), [5.0, 5.0]);
    assert!(!frame.reached, "the held point cannot reach its cursor");
    assert!(
        (frame.tracking_error - (5.0_f64).hypot(5.0)).abs() < 1e-12,
        "tracking error {:.3}",
        frame.tracking_error
    );
    sk.end_drag(drag);
}

/// A frame that jumps a long way across a strongly curved constraint set
/// still arrives: the chain of tangent circles is folded in one frame from
/// its authored arc to a cursor four radii away.
#[test]
fn a_jump_across_a_chain_arrives() {
    let (mut sk, ids) = circle_chain(5, 1.0);
    assert!(sk.solve().converged, "the fixture");
    let last = PointRef::Center(*ids.last().expect("a chain of five"));
    let from = sk.point_coords(last);
    let to = (from.x - 4.0, from.y - 1.0);
    let drag = sk.begin_drag(&[last]);
    let frame = sk.drag_to(&drag, &[to]);
    assert!(
        frame.reached,
        "the jump stopped {:.3e} short, {:?}",
        frame.tracking_error, frame.status
    );
    assert!(frame.residual <= 1e-9, "residual {:.3e}", frame.residual);
    let result = sk.end_drag(drag);
    assert!(
        result.converged,
        "the final solve left {:.3e}",
        result.residual
    );
}

/// A point the constraints hold — triangulated by two distances from two
/// fixed points — does not move, however far the cursor goes, and the
/// frame says it did not arrive. The constraints are left exactly as they
/// were.
///
/// Distances rather than a coincidence with a fixed point on purpose: a
/// coincidence is eliminated as an alias before the solve, and the alias
/// snap would absorb a bad step on its own, hiding whether the drag
/// declined to take one.
#[test]
fn a_point_the_constraints_hold_does_not_move() {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(5.0, 0.0);
    let held = sk.add_point(4.0, 3.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Point(a), 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(PointRef::Point(b), 5.0, 0.0));
    sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(held),
        5.0,
    ));
    sk.add_constraint(Constraint::point_distance(
        PointRef::Point(b),
        PointRef::Point(held),
        10.0_f64.sqrt(),
    ));
    assert!(sk.solve().converged, "the fixture");
    let before = sk.params_snapshot();
    let drag = sk.begin_drag(&[sk.point(held)]);
    let frame = sk.drag_to(&drag, &[(9.0, -4.0)]);
    assert!(!frame.reached, "a held point cannot reach a cursor");
    assert!(
        frame.status <= SolveStatus::ConvergedToPrecision,
        "the frame failed: {:?}",
        frame.status
    );
    assert!(frame.residual <= 1e-12, "residual {:.3e}", frame.residual);
    assert_eq!(
        sk.params_snapshot(),
        before,
        "the sketch moved although nothing could"
    );
    // And it cost no solve: a handle the constraints hold is seen for what
    // it is before anything is stepped or corrected.
    assert_eq!(
        frame.iterations, 0,
        "a held handle was stepped and then corrected back"
    );
    sk.end_drag(drag);
}

/// One frame that jumps across the drawing arrives where a frame-by-frame
/// drag along the same route arrives: the walk re-linearises instead of
/// aiming once at a tangent that has already turned.
#[test]
fn a_jump_arrives_where_the_walk_arrives() {
    let (mut sk, p) = four_bar(10.0, 4.0, 7.0, 5.0);
    assert!(sk.solve().converged, "the fixture");
    let handle = PointRef::Point(p[2]);
    let start = sk.params_snapshot();

    // Walk the coupler point along its curve in small steps, by asking for
    // a cursor a little further each time; where it ends up is on the
    // curve by construction.
    let mut walked = sk.clone();
    let drag = walked.begin_drag(&[handle]);
    let from = walked.point_coords(handle);
    for i in 1..=40 {
        let t = i as f64 / 40.0;
        walked.drag_to(&drag, &[(from.x + 2.0 * t, from.y - 4.0 * t)]);
    }
    walked.end_drag(drag);
    let target = walked.point_coords(handle);
    let reached = walked.params_snapshot();

    // The same target in one frame.
    sk.set_params(start);
    let drag = sk.begin_drag(&[handle]);
    let frame = sk.drag_to(&drag, &[(target.x, target.y)]);
    sk.end_drag(drag);
    assert!(
        frame.reached,
        "the jump stopped {:.3e} short of a reachable target",
        frame.tracking_error
    );
    let worst = reached
        .iter()
        .zip(sk.params_snapshot().iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        worst <= 1e-6,
        "the jump landed {worst:.3e} from where the walk did"
    );
}

/// Two handles dragged at once, each to its own cursor.
#[test]
fn two_handles_drag_at_once() {
    let (mut sk, _, lines) = free_rectangle(4.0, 3.0);
    assert!(sk.solve().converged, "the fixture");
    // The bottom-right corner owns the width, the top-left the height.
    let (width, height) = (PointRef::LineEnd(lines[0]), PointRef::LineEnd(lines[2]));
    let drag = sk.begin_drag(&[width, height]);
    let frame = sk.drag_to(&drag, &[(6.0, 0.0), (0.0, 5.0)]);
    assert!(frame.reached, "both cursors are reachable");
    assert!(frame.residual <= 1e-9, "residual {:.3e}", frame.residual);
    let bottom_right = sk.point_coords(width);
    let top_left = sk.point_coords(height);
    assert!(
        (bottom_right.x - 6.0).abs() < 1e-7,
        "width {bottom_right:?}"
    );
    assert!((top_left.y - 5.0).abs() < 1e-7, "height {top_left:?}");
    sk.end_drag(drag);
}

/// A token from a drag that is over moves nothing.
#[test]
fn a_stale_token_drags_nothing() {
    let (mut sk, _, lines) = free_rectangle(4.0, 3.0);
    let handle = PointRef::LineEnd(lines[1]);
    let first = sk.begin_drag(&[handle]);
    let second = sk.begin_drag(&[handle]);
    let before = sk.params_snapshot();
    let frame = sk.drag_to(&first, &[(6.0, 4.0)]);
    assert_eq!(frame.status, SolveStatus::Stalled);
    assert!(!frame.reached);
    assert_eq!(
        sk.params_snapshot(),
        before,
        "a stale token moved the sketch"
    );
    let frame = sk.drag_to(&second, &[(6.0, 4.0)]);
    assert!(frame.reached, "the live token drags");
    sk.end_drag(second);
}

/// The two things a host could do before this API existed, for comparison:
/// write the cursor into the dragged point and re-solve, or add a
/// temporary `Fixed` at the cursor.
fn run_without_the_api(case: &DragCase, temp_fixed: bool) -> Run {
    let mut sk = case.sketch.clone();
    let start = sk.params_snapshot();
    let before = chirality(&sk, &case.triples);
    let mut signature = before.clone();
    let (xi, yi) = case.handle.param_indices(&sk.entities);
    let mut run = Run {
        frames: 0,
        unconverged: 0,
        reached: 0,
        flips: 0,
        worst_residual: 0.0,
        worst_tracking: 0.0,
        round_trip: 0.0,
        returned: true,
        times_ms: Vec::new(),
    };
    for &(cx, cy) in &case.path {
        let t0 = std::time::Instant::now();
        let result = if temp_fixed {
            let slot = sk.add_constraint(Fixed {
                point: case.handle,
                x: cx,
                y: cy,
            });
            let r = sk.solve();
            sk.remove_constraint(slot);
            r
        } else {
            let mut params = sk.params_snapshot();
            params[xi] = cx;
            params[yi] = cy;
            sk.set_params(params);
            sk.solve()
        };
        run.times_ms.push(t0.elapsed().as_secs_f64() * 1e3);
        run.frames += 1;
        if !result.converged {
            run.unconverged += 1;
        }
        let at = sk.point_coords(case.handle);
        let error = (at.x - cx).hypot(at.y - cy);
        if error <= 1e-9 {
            run.reached += 1;
        }
        if case.reachable {
            run.worst_tracking = run.worst_tracking.max(error);
        }
        run.worst_residual = run.worst_residual.max(result.residual);
        let now = chirality(&sk, &case.triples);
        if signature
            .iter()
            .zip(&now)
            .any(|(was, is)| *was != 0 && *is != 0 && was != is)
        {
            run.flips += 1;
        }
        for (slot, s) in signature.iter_mut().zip(&now) {
            if *s != 0 {
                *slot = *s;
            }
        }
    }
    run.round_trip = start
        .iter()
        .zip(sk.params_snapshot().iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    run.returned = chirality(&sk, &case.triples) == before;
    run
}

#[test]
#[ignore = "measurement, not an assertion"]
fn report_drag_frames() {
    for case in &drag_cases() {
        let runs = [
            ("drag API", run_drag(case)),
            ("warm-set", run_without_the_api(case, false)),
            ("temp-Fixed", run_without_the_api(case, true)),
        ];
        for (label, run) in runs {
            let mut times = run.times_ms.clone();
            times.sort_by(f64::total_cmp);
            println!(
                "{:<24} {:<11} frames {:3}  unconverged {:3}  reached {:3}  flips {:3}  \
                 returned {:<5} round trip {:8.2e}  track {:8.2e}  residual {:8.2e}  \
                 ms p50 {:7.3} p99 {:7.3}",
                case.id,
                label,
                run.frames,
                run.unconverged,
                run.reached,
                run.flips,
                run.returned,
                run.round_trip,
                run.worst_tracking,
                run.worst_residual,
                percentile(&times, 0.5),
                percentile(&times, 0.99),
            );
        }
    }
}

/// A rigid triangle dragged by its apex: the apex follows the cursor
/// whatever happens, and what stays free is the rotation about it, which
/// trades one base point's motion against the other's. That is the choice
/// the weights make. Uniform weights translate the triangle, so both base
/// points move by exactly what the apex did; a heavier point gives less
/// and the other gives more.
#[test]
fn a_heavier_entity_gives_less() {
    let triangle = || {
        let mut sk = Sketch::new();
        let a = sk.add_point(0.0, 0.0);
        let b = sk.add_point(3.8, 0.0);
        let apex = sk.add_point(1.9, 0.6245);
        for (p, q, d) in [(a, apex, 2.0), (apex, b, 2.0), (a, b, 3.8)] {
            sk.add_constraint(Constraint::point_distance(
                PointRef::Point(p),
                PointRef::Point(q),
                d,
            ));
        }
        assert!(sk.solve().converged, "the triangle fixture");
        (sk, a, b, apex)
    };
    // How far each base point travels when the apex is lifted one unit.
    let lift = |weight_a: f64, weight_b: f64| {
        let (mut sk, a, b, apex) = triangle();
        assert!(sk.set_entity_weight(a, weight_a));
        assert!(sk.set_entity_weight(b, weight_b));
        assert_eq!(sk.entity_weight(a), weight_a);
        let before = sk.params_snapshot();
        let at = sk.point_coords(sk.point(apex));
        let drag = sk.begin_drag(&[sk.point(apex)]);
        let frame = sk.drag_to(&drag, &[(at.x, at.y + 1.0)]);
        assert!(frame.reached, "the apex can be lifted");
        assert!(frame.residual <= 1e-8, "residual {:.3e}", frame.residual);
        sk.end_drag(drag);
        let after = sk.params_snapshot();
        let moved = |id| {
            let r = sk.entity_param_range(id);
            (after[r.start] - before[r.start]).hypot(after[r.start + 1] - before[r.start + 1])
        };
        (moved(a), moved(b))
    };

    let (a_even, b_even) = lift(1.0, 1.0);
    assert!(
        (a_even - 1.0).abs() < 1e-6 && (b_even - 1.0).abs() < 1e-6,
        "uniform weights should translate the triangle: {a_even:.4} and {b_even:.4}"
    );
    let (a_heavy, b_heavy) = lift(4.0, 1.0);
    assert!(
        a_heavy < 0.6 && b_heavy > 1.7,
        "a four times heavier A should give less: {a_heavy:.4} and {b_heavy:.4}"
    );
    let (a_light, b_light) = lift(1.0, 4.0);
    assert!(
        (a_light - b_heavy).abs() < 1e-6 && (b_light - a_heavy).abs() < 1e-6,
        "the mirrored weights should give the mirrored motion: {a_light:.4} and {b_light:.4}"
    );
}

/// A weight is neither a hold nor a constraint, and `solve()` does not
/// read it: a very heavy entity still moves where the constraints say it
/// must, and the solution of a sketch is what it was.
#[test]
fn a_weight_is_not_a_hold() {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(1.0, 0.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Point(a), 0.0, 0.0));
    let span = sk.add_constraint(Constraint::point_distance(
        PointRef::Point(a),
        PointRef::Point(b),
        1.0,
    ));
    assert!(sk.solve().converged, "the fixture");
    let mut heavy = sk.clone();
    assert!(heavy.set_entity_weight(b, 1.0e6));
    // The same dimension change, with and without the weight.
    for sketch in [&mut sk, &mut heavy] {
        assert!(sketch.set_constraint_value(span, 4.0));
        assert!(sketch.solve().converged);
    }
    assert_eq!(
        sk.params_snapshot(),
        heavy.params_snapshot(),
        "the weight changed what solve() does"
    );
    assert!((sk.point_coords(sk.point(b)).x - 4.0).abs() < 1e-9);
    // Dragging it is still possible, since nothing holds it.
    let drag = heavy.begin_drag(&[heavy.point(b)]);
    let frame = heavy.drag_to(&drag, &[(0.0, 4.0)]);
    heavy.end_drag(drag);
    assert!(frame.reached, "a heavy point is not a held one");
}

/// A weight has to be a positive, finite number; anything else is refused
/// and changes nothing.
#[test]
fn a_weight_has_to_be_positive_and_finite() {
    let mut sk = Sketch::new();
    let p = sk.add_point(0.0, 0.0);
    assert_eq!(sk.entity_weight(p), 1.0, "the default weight");
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(!sk.set_entity_weight(p, bad), "{bad} was accepted");
        assert_eq!(sk.entity_weight(p), 1.0, "{bad} changed the weight");
    }
    assert!(
        !sk.set_entity_weight(p + 1, 2.0),
        "an entity that is not there"
    );
    assert!(sk.set_entity_weight(p, 2.5));
    assert_eq!(sk.entity_weight(p), 2.5);
}
