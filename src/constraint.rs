//! Geometric constraints and their residual equations.
//!
//! Each [`Constraint`] variant encodes one or more scalar equations f(x)=0.
//! The [`Constraint::residuals`] method evaluates those equations given the
//! current parameter vector, and [`Constraint::equation_count`] returns how
//! many equations the constraint contributes.

use crate::entity::{Entity, EntityId, EntityKind, PointRef};
use crate::solver::SOLVER_NORM_TOL;

/// A geometric constraint between sketch entities.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub enum Constraint {
    // ── Point constraints ────────────────────────────────────────────────────
    /// Two points occupy the same location. (2 equations)
    Coincident(PointRef, PointRef),

    /// A point is fixed at a specific (x, y) location. (2 equations)
    Fixed { point: PointRef, x: f64, y: f64 },

    // ── Distance / length ────────────────────────────────────────────────────
    /// Euclidean distance between two points equals `distance`. (1 equation)
    PointDistance {
        p1: PointRef,
        p2: PointRef,
        distance: f64,
    },

    /// Directed distance between two points along the X axis: `x1 − x2`
    /// equals `distance`. (1 equation)
    ///
    /// The residual is the signed difference, not `|Δx| − distance`: the
    /// linear form keeps the residual smooth at `Δx = 0` (an absolute-value
    /// residual has a vanishing/discontinuous Jacobian entry exactly at the
    /// symmetric initial guesses sketches commonly start from), and the
    /// finite-difference Jacobian integrates it exactly. An edge margin is
    /// expressed by the argument order — `DistanceX(edge_point, hole, m)`
    /// places the hole `m` units before `edge_point` along −X.
    DistanceX {
        p1: PointRef,
        p2: PointRef,
        distance: f64,
    },

    /// Directed distance between two points along the Y axis: `y1 − y2`
    /// equals `distance`. (1 equation)
    ///
    /// Signed for the same reason as [`Constraint::DistanceX`].
    DistanceY {
        p1: PointRef,
        p2: PointRef,
        distance: f64,
    },

    /// Length of a Line entity equals `length`. (1 equation)
    LineLength { line: EntityId, length: f64 },

    // ── Line orientation ─────────────────────────────────────────────────────
    /// A Line entity is horizontal (y1 == y2). (1 equation)
    Horizontal(EntityId),

    /// A Line entity is vertical (x1 == x2). (1 equation)
    Vertical(EntityId),

    // ── Line–line relations ──────────────────────────────────────────────────
    /// Two Line entities are parallel. (1 equation: the sine of the angle
    /// between them, so either sense of the second line satisfies it)
    Parallel(EntityId, EntityId),

    /// Two Line entities are perpendicular. (1 equation: the cosine of the
    /// angle between them)
    Perpendicular(EntityId, EntityId),

    /// Two Line entities have equal length. (1 equation: the difference of
    /// the lengths)
    EqualLength(EntityId, EntityId),

    /// Angle from line `l1` to line `l2` equals `angle_rad`. (1 equation:
    /// `sin(angle_rad − θ)` for the angle `θ` the lines make, which also
    /// holds at `angle_rad + π`)
    ///
    /// The three angular relations are stated on the unit directions of the
    /// lines and multiplied by the geometric mean of the two lengths, so
    /// that their residual is a length like that of every other 2D
    /// equation: roughly how far the end of a line is from where the
    /// relation would put it.
    Angle {
        l1: EntityId,
        l2: EntityId,
        angle_rad: f64,
    },

    // ── Circle / arc constraints ─────────────────────────────────────────────
    /// Two Circle or Arc entities have equal radius. (1 equation)
    EqualRadius(EntityId, EntityId),

    /// A point lies on a circle (or arc circle). (1 equation)
    PointOnCircle { point: PointRef, circle: EntityId },

    /// A point lies on a line (infinite extension). (1 equation)
    PointOnLine { point: PointRef, line: EntityId },

    /// A circle is tangent to a line. (1 equation)
    Tangent { circle: EntityId, line: EntityId },

    /// Two circles are externally or internally tangent. (1 equation)
    ///
    /// External tangency: dist(c1, c2) = r1 + r2.
    /// Internal tangency: dist(c1, c2) = |r1 - r2|.
    /// The solver finds whichever solution is closest to the initial guess.
    /// Use `external: true` for external tangency, `false` for internal.
    CircleCircleTangent {
        c1: EntityId,
        c2: EntityId,
        external: bool,
    },

    /// Circle radius equals a fixed value. (1 equation)
    Radius { circle: EntityId, radius: f64 },

    /// Two arcs (or arc and circle) are tangent. (1 equation)
    ///
    /// Uses the same residual as `CircleCircleTangent` — arcs share the same
    /// center/radius parameter layout as circles.
    ArcArcTangent {
        a1: EntityId,
        a2: EntityId,
        external: bool,
    },

    /// Two points are symmetric about a line. (2 equations)
    ///
    /// The midpoint of p1–p2 lies on the line, and the segment p1–p2 is
    /// perpendicular to the line direction.
    Symmetric {
        p1: PointRef,
        p2: PointRef,
        line: EntityId,
    },

    /// Two points are symmetric about a third point: `p1 + p2 = 2·about`.
    /// (2 equations)
    ///
    /// The point-centered counterpart of [`Constraint::Symmetric`] — e.g.
    /// keeping a bore centered between two skeleton endpoints.
    SymmetricPoint {
        p1: PointRef,
        p2: PointRef,
        about: PointRef,
    },

    // ── Circle / arc relations ────────────────────────────────────────────────
    /// Two Circle/Arc entities share the same center. (2 equations)
    Concentric(EntityId, EntityId),

    /// A point is the midpoint of a Line segment. (2 equations)
    Midpoint { point: PointRef, line: EntityId },

    /// Circle/Arc diameter equals a fixed value. (1 equation)
    Diameter { circle: EntityId, diameter: f64 },
}

impl Constraint {
    // ── Convenience constructors ─────────────────────────────────────────────

    pub fn fix_point(point: impl Into<PointRef>, x: f64, y: f64) -> Self {
        Constraint::Fixed {
            point: point.into(),
            x,
            y,
        }
    }

    pub fn point_distance(p1: impl Into<PointRef>, p2: impl Into<PointRef>, distance: f64) -> Self {
        Constraint::PointDistance {
            p1: p1.into(),
            p2: p2.into(),
            distance,
        }
    }

    pub fn distance_x(p1: impl Into<PointRef>, p2: impl Into<PointRef>, distance: f64) -> Self {
        Constraint::DistanceX {
            p1: p1.into(),
            p2: p2.into(),
            distance,
        }
    }

    pub fn distance_y(p1: impl Into<PointRef>, p2: impl Into<PointRef>, distance: f64) -> Self {
        Constraint::DistanceY {
            p1: p1.into(),
            p2: p2.into(),
            distance,
        }
    }

    pub fn coincident(p1: impl Into<PointRef>, p2: impl Into<PointRef>) -> Self {
        Constraint::Coincident(p1.into(), p2.into())
    }

    // ── Equation count ───────────────────────────────────────────────────────

    /// Number of scalar equations this constraint contributes.
    pub fn equation_count(&self) -> usize {
        match self {
            Constraint::Coincident(..)
            | Constraint::Fixed { .. }
            | Constraint::Symmetric { .. }
            | Constraint::SymmetricPoint { .. }
            | Constraint::Concentric(..)
            | Constraint::Midpoint { .. } => 2,
            _ => 1,
        }
    }

    // ── Residual evaluation ──────────────────────────────────────────────────

    /// Evaluate the constraint residuals f(x) into `out`.
    ///
    /// `out` must have length == `equation_count()`.
    /// Returns `false` if the constraint references an entity of the wrong kind.
    pub fn residuals(&self, params: &[f64], entities: &[Entity], out: &mut [f64]) -> bool {
        match self {
            // ── Coincident ───────────────────────────────────────────────────
            Constraint::Coincident(p1, p2) => {
                let (x1, y1) = p1.point_coords(params, entities);
                let (x2, y2) = p2.point_coords(params, entities);
                out[0] = x1 - x2;
                out[1] = y1 - y2;
            }

            // ── Fixed ────────────────────────────────────────────────────────
            Constraint::Fixed { point, x, y } => {
                let (px, py) = point.point_coords(params, entities);
                out[0] = px - x;
                out[1] = py - y;
            }

            // ── PointDistance ────────────────────────────────────────────────
            Constraint::PointDistance { p1, p2, distance } => {
                let (x1, y1) = p1.point_coords(params, entities);
                let (x2, y2) = p2.point_coords(params, entities);
                let dx = x1 - x2;
                let dy = y1 - y2;
                // Use sqrt form for better conditioning near solution.
                out[0] = (dx * dx + dy * dy).sqrt() - distance;
            }

            // ── DistanceX ────────────────────────────────────────────────────
            Constraint::DistanceX { p1, p2, distance } => {
                let (x1, _) = p1.point_coords(params, entities);
                let (x2, _) = p2.point_coords(params, entities);
                out[0] = x1 - x2 - distance;
            }

            // ── DistanceY ────────────────────────────────────────────────────
            Constraint::DistanceY { p1, p2, distance } => {
                let (_, y1) = p1.point_coords(params, entities);
                let (_, y2) = p2.point_coords(params, entities);
                out[0] = y1 - y2 - distance;
            }

            // ── LineLength ───────────────────────────────────────────────────
            Constraint::LineLength { line, length } => {
                let e = &entities[*line];
                let dx = params[e.param(2)] - params[e.param(0)];
                let dy = params[e.param(3)] - params[e.param(1)];
                out[0] = (dx * dx + dy * dy).sqrt() - length;
            }

            // ── Horizontal ───────────────────────────────────────────────────
            Constraint::Horizontal(id) => {
                let e = &entities[*id];
                out[0] = params[e.param(1)] - params[e.param(3)]; // y1 - y2
            }

            // ── Vertical ─────────────────────────────────────────────────────
            Constraint::Vertical(id) => {
                let e = &entities[*id];
                out[0] = params[e.param(0)] - params[e.param(2)]; // x1 - x2
            }

            // ── Line–line relations ───────────────────────────────────────────
            // Stated on the unit directions and measured in the lines' own
            // length (see `angular_length`), so that each is a length
            // whatever the size of the lines.
            Constraint::Parallel(l1, l2) => {
                let (a, b) = (
                    line_direction(params, &entities[*l1]),
                    line_direction(params, &entities[*l2]),
                );
                // sin of the angle from l1 to l2
                out[0] = angular_length(&a, &b) * (a.x * b.y - a.y * b.x);
            }
            Constraint::Perpendicular(l1, l2) => {
                let (a, b) = (
                    line_direction(params, &entities[*l1]),
                    line_direction(params, &entities[*l2]),
                );
                // cos of the angle from l1 to l2
                out[0] = angular_length(&a, &b) * (a.x * b.x + a.y * b.y);
            }
            Constraint::EqualLength(l1, l2) => {
                let (a, b) = (
                    line_direction(params, &entities[*l1]),
                    line_direction(params, &entities[*l2]),
                );
                out[0] = a.length - b.length;
            }
            Constraint::Angle { l1, l2, angle_rad } => {
                let (a, b) = (
                    line_direction(params, &entities[*l1]),
                    line_direction(params, &entities[*l2]),
                );
                // sin(angle − θ), θ the angle from l1 to l2:
                // sin(angle)·cos θ − cos(angle)·sin θ
                let dot = a.x * b.x + a.y * b.y;
                let cross = a.x * b.y - a.y * b.x;
                out[0] = angular_length(&a, &b) * (angle_rad.sin() * dot - angle_rad.cos() * cross);
            }

            // ── EqualRadius ───────────────────────────────────────────────────
            Constraint::EqualRadius(c1, c2) => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                // radius is param[2] for both Circle and Arc
                out[0] = params[e1.param(2)] - params[e2.param(2)];
            }

            // ── Radius ────────────────────────────────────────────────────────
            Constraint::Radius { circle, radius } => {
                let e = &entities[*circle];
                out[0] = params[e.param(2)] - radius;
            }

            // ── PointOnCircle ─────────────────────────────────────────────────
            Constraint::PointOnCircle { point, circle } => {
                let (px, py) = point.point_coords(params, entities);
                let ec = &entities[*circle];
                let cx = params[ec.param(0)];
                let cy = params[ec.param(1)];
                let r = params[ec.param(2)];
                let dx = px - cx;
                let dy = py - cy;
                out[0] = (dx * dx + dy * dy).sqrt() - r;
            }

            // ── PointOnLine ───────────────────────────────────────────────────
            Constraint::PointOnLine { point, line } => {
                let (px, py) = point.point_coords(params, entities);
                let el = &entities[*line];
                let x1 = params[el.param(0)];
                let y1 = params[el.param(1)];
                let x2 = params[el.param(2)];
                let y2 = params[el.param(3)];
                // Normalize by line length to improve conditioning.
                let len = ((x2 - x1).powi(2) + (y2 - y1).powi(2))
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                out[0] = ((x2 - x1) * (py - y1) - (y2 - y1) * (px - x1)) / len;
            }

            // ── Tangent (circle–line) ─────────────────────────────────────────
            Constraint::Tangent { circle, line } => {
                let ec = &entities[*circle];
                let el = &entities[*line];
                let cx = params[ec.param(0)];
                let cy = params[ec.param(1)];
                let r = params[ec.param(2)];
                let x1 = params[el.param(0)];
                let y1 = params[el.param(1)];
                let x2 = params[el.param(2)];
                let y2 = params[el.param(3)];
                let dx = x2 - x1;
                let dy = y2 - y1;
                let len = (dx * dx + dy * dy).sqrt().max(SOLVER_NORM_TOL);
                // dist(center, line) - r = 0
                let num = (cx - x1) * dy - (cy - y1) * dx;
                out[0] = num / len - r;
            }

            // ── CircleCircleTangent ───────────────────────────────────────────
            Constraint::CircleCircleTangent { c1, c2, external } => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                let cx1 = params[e1.param(0)];
                let cy1 = params[e1.param(1)];
                let r1 = params[e1.param(2)];
                let cx2 = params[e2.param(0)];
                let cy2 = params[e2.param(1)];
                let r2 = params[e2.param(2)];
                let dx = cx2 - cx1;
                let dy = cy2 - cy1;
                let dist = (dx * dx + dy * dy).sqrt();
                if *external {
                    // External tangency: dist = r1 + r2
                    out[0] = dist - (r1 + r2);
                } else {
                    // Internal tangency: dist = |r1 - r2|
                    out[0] = dist - (r1 - r2).abs();
                }
            }

            // ── ArcArcTangent ─────────────────────────────────────────────────
            Constraint::ArcArcTangent { a1, a2, external } => {
                let e1 = &entities[*a1];
                let e2 = &entities[*a2];
                let cx1 = params[e1.param(0)];
                let cy1 = params[e1.param(1)];
                let r1 = params[e1.param(2)];
                let cx2 = params[e2.param(0)];
                let cy2 = params[e2.param(1)];
                let r2 = params[e2.param(2)];
                let dx = cx2 - cx1;
                let dy = cy2 - cy1;
                let dist = (dx * dx + dy * dy).sqrt();
                if *external {
                    out[0] = dist - (r1 + r2);
                } else {
                    out[0] = dist - (r1 - r2).abs();
                }
            }

            // ── Symmetric ─────────────────────────────────────────────────────
            Constraint::Symmetric { p1, p2, line } => {
                let (x1, y1) = p1.point_coords(params, entities);
                let (x2, y2) = p2.point_coords(params, entities);
                let el = &entities[*line];
                let lx1 = params[el.param(0)];
                let ly1 = params[el.param(1)];
                let lx2 = params[el.param(2)];
                let ly2 = params[el.param(3)];
                let ldx = lx2 - lx1;
                let ldy = ly2 - ly1;
                let len = (ldx * ldx + ldy * ldy).sqrt().max(SOLVER_NORM_TOL);
                // Midpoint of p1–p2
                let mx = (x1 + x2) * 0.5;
                let my = (y1 + y2) * 0.5;
                // Eq 1: midpoint lies on the line (cross product = 0)
                out[0] = ((lx2 - lx1) * (my - ly1) - (ly2 - ly1) * (mx - lx1)) / len;
                // Eq 2: segment p1–p2 is perpendicular to line direction (dot = 0)
                let sdx = x2 - x1;
                let sdy = y2 - y1;
                out[1] = (sdx * ldx + sdy * ldy) / len;
            }

            // ── SymmetricPoint ───────────────────────────────────────────────
            Constraint::SymmetricPoint { p1, p2, about } => {
                let (x1, y1) = p1.point_coords(params, entities);
                let (x2, y2) = p2.point_coords(params, entities);
                let (xa, ya) = about.point_coords(params, entities);
                out[0] = x1 + x2 - 2.0 * xa;
                out[1] = y1 + y2 - 2.0 * ya;
            }

            // ── Concentric ────────────────────────────────────────────────────
            Constraint::Concentric(c1, c2) => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                out[0] = params[e1.param(0)] - params[e2.param(0)];
                out[1] = params[e1.param(1)] - params[e2.param(1)];
            }

            // ── Midpoint ──────────────────────────────────────────────────────
            Constraint::Midpoint { point, line } => {
                let (px, py) = point.point_coords(params, entities);
                let el = &entities[*line];
                let mx = (params[el.param(0)] + params[el.param(2)]) * 0.5;
                let my = (params[el.param(1)] + params[el.param(3)]) * 0.5;
                out[0] = px - mx;
                out[1] = py - my;
            }

            // ── Diameter ──────────────────────────────────────────────────────
            Constraint::Diameter { circle, diameter } => {
                let e = &entities[*circle];
                out[0] = params[e.param(2)] - diameter / 2.0;
            }
        }
        true
    }

    /// Analytic Jacobian of this constraint's equations — the companion of
    /// [`Constraint::residuals`], kept in lockstep with it by the generic
    /// `analytic_jacobian_matches_finite_difference` test (a variant added
    /// without derivatives fails CI).
    ///
    /// `out` must have length == `equation_count()`; `out[r]` receives the
    /// `(absolute parameter index, ∂rᵣ/∂param)` pairs for equation row `r`.
    /// Rows are sparse — a constraint touches 1–8 parameters. Parameter
    /// indices held by the solver's fixed mask are filtered by the caller.
    pub fn jacobian(&self, params: &[f64], entities: &[Entity], out: &mut [Vec<(usize, f64)>]) {
        // Push (param, derivative) entries; local shorthand for readability.
        macro_rules! row {
            ($r:expr, $($p:expr => $d:expr),+ $(,)?) => {{
                $( $r.push(($p, $d)); )+
            }};
        }
        match self {
            // ── Coincident ───────────────────────────────────────────────────
            Constraint::Coincident(p1, p2) => {
                let deps1 = point_param_deps(p1, params, entities);
                let deps2 = point_param_deps(p2, params, entities);
                // r0 = px1 - px2, r1 = py1 - py2
                accumulate_point_jacobian(&mut out[0], 1.0, 0.0, &deps1);
                accumulate_point_jacobian(&mut out[0], -1.0, 0.0, &deps2);
                accumulate_point_jacobian(&mut out[1], 0.0, 1.0, &deps1);
                accumulate_point_jacobian(&mut out[1], 0.0, -1.0, &deps2);
            }

            // ── Fixed ────────────────────────────────────────────────────────
            Constraint::Fixed { point, .. } => {
                let deps = point_param_deps(point, params, entities);
                // r0 = px - x, r1 = py - y
                accumulate_point_jacobian(&mut out[0], 1.0, 0.0, &deps);
                accumulate_point_jacobian(&mut out[1], 0.0, 1.0, &deps);
            }

            // ── PointDistance ────────────────────────────────────────────────
            Constraint::PointDistance { p1, p2, .. } => {
                let (x1, y1) = p1.point_coords(params, entities);
                let (x2, y2) = p2.point_coords(params, entities);
                let dx = x1 - x2;
                let dy = y1 - y2;
                // Guarded: at a coincident pair the residual's |·| kink makes
                // the derivative undefined; a zero gradient (absorbed by the
                // Tikhonov minimum-norm step) beats NaN poisoning the whole
                // normal-equation system.
                let l = (dx * dx + dy * dy).sqrt().max(SOLVER_NORM_TOL);
                let deps1 = point_param_deps(p1, params, entities);
                let deps2 = point_param_deps(p2, params, entities);
                accumulate_point_jacobian(&mut out[0], dx / l, dy / l, &deps1);
                accumulate_point_jacobian(&mut out[0], -dx / l, -dy / l, &deps2);
            }

            // ── DistanceX / DistanceY ────────────────────────────────────────
            Constraint::DistanceX { p1, p2, .. } => {
                let deps1 = point_param_deps(p1, params, entities);
                let deps2 = point_param_deps(p2, params, entities);
                // r0 = px1 - px2 - distance
                accumulate_point_jacobian(&mut out[0], 1.0, 0.0, &deps1);
                accumulate_point_jacobian(&mut out[0], -1.0, 0.0, &deps2);
            }
            Constraint::DistanceY { p1, p2, .. } => {
                let deps1 = point_param_deps(p1, params, entities);
                let deps2 = point_param_deps(p2, params, entities);
                // r0 = py1 - py2 - distance
                accumulate_point_jacobian(&mut out[0], 0.0, 1.0, &deps1);
                accumulate_point_jacobian(&mut out[0], 0.0, -1.0, &deps2);
            }

            // ── LineLength ───────────────────────────────────────────────────
            Constraint::LineLength { line, .. } => {
                let e = &entities[*line];
                let (p0, p1, p2, p3) = (e.param(0), e.param(1), e.param(2), e.param(3));
                let dx = params[p2] - params[p0];
                let dy = params[p3] - params[p1];
                let l = (dx * dx + dy * dy).sqrt();
                row!(out[0], p0 => -dx / l, p1 => -dy / l, p2 => dx / l, p3 => dy / l);
            }

            // ── Horizontal / Vertical ────────────────────────────────────────
            Constraint::Horizontal(id) => {
                let e = &entities[*id];
                row!(out[0], e.param(1) => 1.0, e.param(3) => -1.0);
            }
            Constraint::Vertical(id) => {
                let e = &entities[*id];
                row!(out[0], e.param(0) => 1.0, e.param(2) => -1.0);
            }

            // ── Line–line relations ───────────────────────────────────────────
            // Each residual is `m·c`, `c` bilinear in the two unit directions
            // `a`, `b` and `m` their `angular_length`; given `∂c/∂a` and
            // `∂c/∂b` the entries are `direction_row`'s.
            Constraint::Parallel(l1, l2) => {
                let (e1, e2) = (&entities[*l1], &entities[*l2]);
                let (a, b) = (line_direction(params, e1), line_direction(params, e2));
                // c = a.x·b.y − a.y·b.x
                let (c, m) = (a.x * b.y - a.y * b.x, angular_length(&a, &b));
                direction_row(&mut out[0], e1, &a, (b.y, -b.x), c, m);
                direction_row(&mut out[0], e2, &b, (-a.y, a.x), c, m);
            }
            Constraint::Perpendicular(l1, l2) => {
                let (e1, e2) = (&entities[*l1], &entities[*l2]);
                let (a, b) = (line_direction(params, e1), line_direction(params, e2));
                // c = a·b
                let (c, m) = (a.x * b.x + a.y * b.y, angular_length(&a, &b));
                direction_row(&mut out[0], e1, &a, (b.x, b.y), c, m);
                direction_row(&mut out[0], e2, &b, (a.x, a.y), c, m);
            }
            Constraint::EqualLength(l1, l2) => {
                let (e1, e2) = (&entities[*l1], &entities[*l2]);
                let (a, b) = (line_direction(params, e1), line_direction(params, e2));
                // r = |d1| − |d2|: each end moves its length along the line.
                row!(out[0],
                    e1.param(0) => -a.x, e1.param(2) => a.x,
                    e1.param(1) => -a.y, e1.param(3) => a.y,
                    e2.param(0) => b.x,  e2.param(2) => -b.x,
                    e2.param(1) => b.y,  e2.param(3) => -b.y,
                );
            }
            Constraint::Angle { l1, l2, angle_rad } => {
                let (e1, e2) = (&entities[*l1], &entities[*l2]);
                let (a, b) = (line_direction(params, e1), line_direction(params, e2));
                // c = sin·(a·b) − cos·(a×b)
                let (sin, cos) = (angle_rad.sin(), angle_rad.cos());
                let c = sin * (a.x * b.x + a.y * b.y) - cos * (a.x * b.y - a.y * b.x);
                let m = angular_length(&a, &b);
                let ga = (sin * b.x - cos * b.y, sin * b.y + cos * b.x);
                let gb = (sin * a.x + cos * a.y, sin * a.y - cos * a.x);
                direction_row(&mut out[0], e1, &a, ga, c, m);
                direction_row(&mut out[0], e2, &b, gb, c, m);
            }

            // ── EqualRadius / Radius / Diameter ───────────────────────────────
            Constraint::EqualRadius(c1, c2) => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                row!(out[0], e1.param(2) => 1.0, e2.param(2) => -1.0);
            }
            Constraint::Radius { circle, .. } => {
                let e = &entities[*circle];
                row!(out[0], e.param(2) => 1.0);
            }
            Constraint::Diameter { circle, .. } => {
                let e = &entities[*circle];
                row!(out[0], e.param(2) => 1.0);
            }

            // ── PointOnCircle ─────────────────────────────────────────────────
            Constraint::PointOnCircle { point, circle } => {
                let (px, py) = point.point_coords(params, entities);
                let ec = &entities[*circle];
                let (cx, cy, cr) = (ec.param(0), ec.param(1), ec.param(2));
                let dx = px - params[cx];
                let dy = py - params[cy];
                let l = (dx * dx + dy * dy).sqrt().max(SOLVER_NORM_TOL);
                let deps = point_param_deps(point, params, entities);
                accumulate_point_jacobian(&mut out[0], dx / l, dy / l, &deps);
                row!(out[0], cx => -dx / l, cy => -dy / l, cr => -1.0);
            }

            // ── PointOnLine ───────────────────────────────────────────────────
            Constraint::PointOnLine { point, line } => {
                let (px, py) = point.point_coords(params, entities);
                let el = &entities[*line];
                let (lx1, ly1, lx2, ly2) = (el.param(0), el.param(1), el.param(2), el.param(3));
                let a = params[lx2] - params[lx1];
                let b = params[ly2] - params[ly1];
                let ux = px - params[lx1];
                let uy = py - params[ly1];
                // r = (A·uy − B·ux)/L; ∂r/∂q = (∂N/∂q)/L − r·(∂L/∂q)/L
                let l = (a * a + b * b).sqrt().max(SOLVER_NORM_TOL);
                let n = a * uy - b * ux;
                let big_r = n / l;
                let deps = point_param_deps(point, params, entities);
                accumulate_point_jacobian(&mut out[0], -b / l, a / l, &deps);
                row!(out[0],
                    lx1 => (b - uy) / l + big_r * a / (l * l),
                    ly1 => (ux - a) / l + big_r * b / (l * l),
                    lx2 => uy / l - big_r * a / (l * l),
                    ly2 => -ux / l - big_r * b / (l * l),
                );
            }

            // ── Tangent (circle–line) ─────────────────────────────────────────
            Constraint::Tangent { circle, line } => {
                let ec = &entities[*circle];
                let el = &entities[*line];
                let (cx, cy, cr) = (ec.param(0), ec.param(1), ec.param(2));
                let (lx1, ly1, lx2, ly2) = (el.param(0), el.param(1), el.param(2), el.param(3));
                let a = params[lx2] - params[lx1];
                let b = params[ly2] - params[ly1];
                let ux = params[cx] - params[lx1];
                let uy = params[cy] - params[ly1];
                // r = (ux·B − uy·A)/L − rc; big_r = N/L
                // ∂N/∂ly1 = A − ux (B depends on ly1!), ∂L/∂ly1 = −B/L
                let l = (a * a + b * b).sqrt().max(SOLVER_NORM_TOL);
                let n = ux * b - uy * a;
                let big_r = n / l;
                row!(out[0],
                    cr => -1.0,
                    cx => b / l,
                    cy => -a / l,
                    lx1 => (uy - b) / l + big_r * a / (l * l),
                    ly1 => (a - ux) / l + big_r * b / (l * l),
                    lx2 => -uy / l - big_r * a / (l * l),
                    ly2 => ux / l - big_r * b / (l * l),
                );
            }

            // ── CircleCircleTangent / ArcArcTangent ───────────────────────────
            Constraint::CircleCircleTangent { .. } | Constraint::ArcArcTangent { .. } => {
                let (e1, e2, external) = match self {
                    Constraint::CircleCircleTangent { c1, c2, external } => {
                        (&entities[*c1], &entities[*c2], *external)
                    }
                    Constraint::ArcArcTangent { a1, a2, external } => {
                        (&entities[*a1], &entities[*a2], *external)
                    }
                    _ => unreachable!(),
                };
                let (cx1, cy1, r1i) = (e1.param(0), e1.param(1), e1.param(2));
                let (cx2, cy2, r2i) = (e2.param(0), e2.param(1), e2.param(2));
                let dx = params[cx2] - params[cx1];
                let dy = params[cy2] - params[cy1];
                let d = (dx * dx + dy * dy).sqrt().max(SOLVER_NORM_TOL);
                let mut row = Vec::with_capacity(6);
                row.push((cx1, -dx / d));
                row.push((cy1, -dy / d));
                row.push((cx2, dx / d));
                row.push((cy2, dy / d));
                if external {
                    row.push((r1i, -1.0));
                    row.push((r2i, -1.0));
                } else {
                    let s = (params[r1i] - params[r2i]).signum();
                    row.push((r1i, -s));
                    row.push((r2i, s));
                }
                out[0] = row;
            }

            // ── Symmetric ─────────────────────────────────────────────────────
            Constraint::Symmetric { p1, p2, line } => {
                let (x1, y1) = p1.point_coords(params, entities);
                let (x2, y2) = p2.point_coords(params, entities);
                let el = &entities[*line];
                let (lx1, ly1, lx2, ly2) = (el.param(0), el.param(1), el.param(2), el.param(3));
                let a = params[lx2] - params[lx1];
                let b = params[ly2] - params[ly1];
                let mx = (x1 + x2) * 0.5;
                let my = (y1 + y2) * 0.5;
                let sdx = x2 - x1;
                let sdy = y2 - y1;
                let l = (a * a + b * b).sqrt().max(SOLVER_NORM_TOL);
                let deps1 = point_param_deps(p1, params, entities);
                let deps2 = point_param_deps(p2, params, entities);
                // Eq 0: r0 = (A·(my − ly1) − B·(mx − lx1))/L
                let big_r0 = (a * (my - params[ly1]) - b * (mx - params[lx1])) / l;
                // Point derivatives for eq 0: dr0/dpx1 = -b/(2L), dr0/dpy1 = a/(2L)
                accumulate_point_jacobian(&mut out[0], -b / (2.0 * l), a / (2.0 * l), &deps1);
                accumulate_point_jacobian(&mut out[0], -b / (2.0 * l), a / (2.0 * l), &deps2);
                row!(out[0],
                    lx1 => (b - my + params[ly1]) / l + big_r0 * a / (l * l),
                    ly1 => (mx - params[lx1] - a) / l + big_r0 * b / (l * l),
                    lx2 => (my - params[ly1]) / l - big_r0 * a / (l * l),
                    ly2 => (params[lx1] - mx) / l - big_r0 * b / (l * l),
                );
                // Eq 1: r1 = (sdx·A + sdy·B)/L
                let big_r1 = (sdx * a + sdy * b) / l;
                // Point derivatives for eq 1: dr1/dpx1 = -a/L, dr1/dpy1 = -b/L
                accumulate_point_jacobian(&mut out[1], -a / l, -b / l, &deps1);
                accumulate_point_jacobian(&mut out[1], a / l, b / l, &deps2);
                row!(out[1],
                    lx1 => -sdx / l + big_r1 * a / (l * l),
                    ly1 => -sdy / l + big_r1 * b / (l * l),
                    lx2 => sdx / l - big_r1 * a / (l * l),
                    ly2 => sdy / l - big_r1 * b / (l * l),
                );
            }

            // ── SymmetricPoint ────────────────────────────────────────────────
            Constraint::SymmetricPoint { p1, p2, about } => {
                let deps1 = point_param_deps(p1, params, entities);
                let deps2 = point_param_deps(p2, params, entities);
                let deps_about = point_param_deps(about, params, entities);
                // r0 = px1 + px2 - 2*pxa, r1 = py1 + py2 - 2*pya
                accumulate_point_jacobian(&mut out[0], 1.0, 0.0, &deps1);
                accumulate_point_jacobian(&mut out[0], 1.0, 0.0, &deps2);
                accumulate_point_jacobian(&mut out[0], -2.0, 0.0, &deps_about);
                accumulate_point_jacobian(&mut out[1], 0.0, 1.0, &deps1);
                accumulate_point_jacobian(&mut out[1], 0.0, 1.0, &deps2);
                accumulate_point_jacobian(&mut out[1], 0.0, -2.0, &deps_about);
            }

            // ── Concentric ────────────────────────────────────────────────────
            Constraint::Concentric(c1, c2) => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                row!(out[0], e1.param(0) => 1.0, e2.param(0) => -1.0);
                row!(out[1], e1.param(1) => 1.0, e2.param(1) => -1.0);
            }

            // ── Midpoint ──────────────────────────────────────────────────────
            Constraint::Midpoint { point, line } => {
                let deps = point_param_deps(point, params, entities);
                let el = &entities[*line];
                // r0 = px - mx, r1 = py - my
                accumulate_point_jacobian(&mut out[0], 1.0, 0.0, &deps);
                accumulate_point_jacobian(&mut out[1], 0.0, 1.0, &deps);
                row!(out[0], el.param(0) => -0.5, el.param(2) => -0.5);
                row!(out[1], el.param(1) => -0.5, el.param(3) => -0.5);
            }
        }
    }
}

// ── From<EntityId> for PointRef (Point entities only) ────────────────────────

impl From<EntityId> for PointRef {
    fn from(id: EntityId) -> Self {
        PointRef::Point(id)
    }
}

/// Compute parameter dependencies of a point reference for Jacobian computation.
/// Returns a list of (param_idx, dpx/dparam, dpy/dparam) tuples.
///
/// For stored points (Point, LineStart, LineEnd, Center), the dependencies are
/// trivial: [(x_idx, 1, 0), (y_idx, 0, 1)].
///
/// For derived points (ArcStart, ArcEnd), the dependencies include the chain
/// rule through center, radius, and angle:
///   x = cx + r*cos(angle), y = cy + r*sin(angle)
///   dx/dcx = 1, dx/dcy = 0, dx/dr = cos(angle), dx/dangle = -r*sin(angle)
///   dy/dcx = 0, dy/dcy = 1, dy/dr = sin(angle), dy/dangle = r*cos(angle)
fn point_param_deps(
    pref: &PointRef,
    params: &[f64],
    entities: &[Entity],
) -> Vec<(usize, f64, f64)> {
    match *pref {
        PointRef::Point(id) => {
            let e = &entities[id];
            vec![(e.param(0), 1.0, 0.0), (e.param(1), 0.0, 1.0)]
        }
        PointRef::LineStart(id) => {
            let e = &entities[id];
            vec![(e.param(0), 1.0, 0.0), (e.param(1), 0.0, 1.0)]
        }
        PointRef::LineEnd(id) => {
            let e = &entities[id];
            vec![(e.param(2), 1.0, 0.0), (e.param(3), 0.0, 1.0)]
        }
        PointRef::Center(id) => {
            let e = &entities[id];
            vec![(e.param(0), 1.0, 0.0), (e.param(1), 0.0, 1.0)]
        }
        PointRef::ArcStart(id) => {
            let e = &entities[id];
            let r = params[e.param(2)];
            let angle = params[e.param(3)];
            let cos_a = angle.cos();
            let sin_a = angle.sin();
            vec![
                (e.param(0), 1.0, 0.0),
                (e.param(1), 0.0, 1.0),
                (e.param(2), cos_a, sin_a),
                (e.param(3), -r * sin_a, r * cos_a),
            ]
        }
        PointRef::ArcEnd(id) => {
            let e = &entities[id];
            let r = params[e.param(2)];
            let angle = params[e.param(4)];
            let cos_a = angle.cos();
            let sin_a = angle.sin();
            vec![
                (e.param(0), 1.0, 0.0),
                (e.param(1), 0.0, 1.0),
                (e.param(2), cos_a, sin_a),
                (e.param(4), -r * sin_a, r * cos_a),
            ]
        }
    }
}

/// The direction of a Line entity: its unit vector, its length and the
/// reciprocal of that.
struct Direction {
    x: f64,
    y: f64,
    length: f64,
    /// `1 / length`, or zero for a line of zero length.
    inverse: f64,
}

/// The [`Direction`] of line `e`. A line of zero length has no direction:
/// it reads as the zero vector with a zero gradient, which leaves the
/// relations it takes part in unmet and still (rather than NaN) until
/// something else gives it a length.
fn line_direction(params: &[f64], e: &Entity) -> Direction {
    let dx = params[e.param(2)] - params[e.param(0)];
    let dy = params[e.param(3)] - params[e.param(1)];
    // `hypot` rather than a sum of squares: no underflow on tiny lines.
    let length = dx.hypot(dy);
    if length > 0.0 {
        let inverse = 1.0 / length;
        Direction {
            x: dx * inverse,
            y: dy * inverse,
            length,
            inverse,
        }
    } else {
        Direction {
            x: 0.0,
            y: 0.0,
            length: 0.0,
            inverse: 0.0,
        }
    }
}

/// The length an angular relation between lines `a` and `b` is measured
/// in: the geometric mean of their lengths. A relation `c(a, b)` on the
/// unit directions is a pure number; stated as `√(|a|·|b|)·c` it is a
/// length, like every other 2D equation, so that one tolerance and one
/// regularisation suit them all at any size. (A product of the two lengths
/// instead would be an area, and would fall under the tolerance on small
/// geometry long before the relation held.)
fn angular_length(a: &Direction, b: &Direction) -> f64 {
    // Two roots rather than the root of the product: no underflow.
    a.length.sqrt() * b.length.sqrt()
}

/// The entries for line `e` of a residual `m·c`, where `c` depends on the
/// line through its unit direction `u` with `g = ∂c/∂u` and is bilinear in
/// the directions (so `g·u = c`), and `m` is the [`angular_length`].
///
/// Moving the end of the line by `δ` turns `u` by `(I − u·uᵀ)·δ / length`
/// and stretches `m` by `½·m·(u·δ) / length`, so the derivative with
/// respect to the end is `m·(g − ½·c·u) / length`. The start moves it the
/// other way.
fn direction_row(
    row: &mut Vec<(usize, f64)>,
    e: &Entity,
    u: &Direction,
    g: (f64, f64),
    c: f64,
    m: f64,
) {
    let dx = m * (g.0 - 0.5 * c * u.x) * u.inverse;
    let dy = m * (g.1 - 0.5 * c * u.y) * u.inverse;
    row.extend([
        (e.param(0), -dx),
        (e.param(1), -dy),
        (e.param(2), dx),
        (e.param(3), dy),
    ]);
}

/// Accumulate Jacobian entries for a point using the chain rule.
/// Given the point derivatives (dr/dpx, dr/dpy) and the parameter dependencies
/// from `point_param_deps`, accumulate the total derivative dr/dparam into `row`.
fn accumulate_point_jacobian(
    row: &mut Vec<(usize, f64)>,
    dr_dpx: f64,
    dr_dpy: f64,
    deps: &[(usize, f64, f64)],
) {
    for &(pidx, dpx, dpy) in deps {
        let deriv = dr_dpx * dpx + dr_dpy * dpy;
        // Merge with existing entry for the same parameter
        if let Some(entry) = row.iter_mut().find(|(p, _)| *p == pidx) {
            entry.1 += deriv;
        } else {
            row.push((pidx, deriv));
        }
    }
}

const LINE: &[EntityKind] = &[EntityKind::Line];
const CIRCULAR: &[EntityKind] = &[EntityKind::Circle, EntityKind::Arc];

impl Constraint {
    /// Every entity this constraint refers to, each with the entity kinds
    /// that reference is valid for. An entity referred to twice is listed
    /// twice.
    pub fn references(&self) -> Vec<(EntityId, &'static [EntityKind])> {
        let point = |p: &PointRef| (p.entity(), p.expected_kinds());
        match self {
            Constraint::Coincident(p1, p2) => vec![point(p1), point(p2)],
            Constraint::Fixed { point: p, .. } => vec![point(p)],
            Constraint::PointDistance { p1, p2, .. }
            | Constraint::DistanceX { p1, p2, .. }
            | Constraint::DistanceY { p1, p2, .. } => vec![point(p1), point(p2)],
            Constraint::LineLength { line, .. } => vec![(*line, LINE)],
            Constraint::Horizontal(line) | Constraint::Vertical(line) => vec![(*line, LINE)],
            Constraint::Parallel(l1, l2)
            | Constraint::Perpendicular(l1, l2)
            | Constraint::EqualLength(l1, l2)
            | Constraint::Angle { l1, l2, .. } => vec![(*l1, LINE), (*l2, LINE)],
            Constraint::EqualRadius(c1, c2)
            | Constraint::Concentric(c1, c2)
            | Constraint::CircleCircleTangent { c1, c2, .. } => {
                vec![(*c1, CIRCULAR), (*c2, CIRCULAR)]
            }
            Constraint::ArcArcTangent { a1, a2, .. } => vec![(*a1, CIRCULAR), (*a2, CIRCULAR)],
            Constraint::PointOnCircle { point: p, circle } => {
                vec![point(p), (*circle, CIRCULAR)]
            }
            Constraint::PointOnLine { point: p, line }
            | Constraint::Midpoint { point: p, line } => {
                vec![point(p), (*line, LINE)]
            }
            Constraint::Tangent { circle, line } => vec![(*circle, CIRCULAR), (*line, LINE)],
            Constraint::Radius { circle, .. } | Constraint::Diameter { circle, .. } => {
                vec![(*circle, CIRCULAR)]
            }
            Constraint::Symmetric { p1, p2, line } => vec![point(p1), point(p2), (*line, LINE)],
            Constraint::SymmetricPoint { p1, p2, about } => {
                vec![point(p1), point(p2), point(about)]
            }
        }
    }

    /// The driving value of a dimensional constraint: the distance, length,
    /// angle (radians), radius or diameter it asks for. `None` for
    /// constraints without a single scalar value (`Fixed` has two; the
    /// geometric relations have none).
    pub fn value(&self) -> Option<f64> {
        match self {
            Constraint::PointDistance { distance, .. }
            | Constraint::DistanceX { distance, .. }
            | Constraint::DistanceY { distance, .. } => Some(*distance),
            Constraint::LineLength { length, .. } => Some(*length),
            Constraint::Angle { angle_rad, .. } => Some(*angle_rad),
            Constraint::Radius { radius, .. } => Some(*radius),
            Constraint::Diameter { diameter, .. } => Some(*diameter),
            _ => None,
        }
    }

    /// Set the driving value of a dimensional constraint (see
    /// [`Constraint::value`]). Returns `false`, changing nothing, for a
    /// constraint that has none.
    pub fn set_value(&mut self, value: f64) -> bool {
        match self {
            Constraint::PointDistance { distance, .. }
            | Constraint::DistanceX { distance, .. }
            | Constraint::DistanceY { distance, .. } => *distance = value,
            Constraint::LineLength { length, .. } => *length = value,
            Constraint::Angle { angle_rad, .. } => *angle_rad = value,
            Constraint::Radius { radius, .. } => *radius = value,
            Constraint::Diameter { diameter, .. } => *diameter = value,
            _ => return false,
        }
        true
    }

    /// Whether every number the constraint carries is finite: its driving
    /// value, or the target position of a `Fixed`.
    pub fn has_finite_values(&self) -> bool {
        match self {
            Constraint::Fixed { x, y, .. } => x.is_finite() && y.is_finite(),
            other => other.value().is_none_or(f64::is_finite),
        }
    }

    /// What a dimensional constraint measures in the current geometry: the
    /// value that [`Constraint::value`] would have to be for the constraint
    /// to hold as things stand. This is how a dimension is read without
    /// driving anything. `None` for constraints without a scalar value.
    ///
    /// An angle is reported from `l1` to `l2` in `(-π, π]`.
    pub fn measure(&self, params: &[f64], entities: &[Entity]) -> Option<f64> {
        let direction = |line: EntityId| {
            let e = &entities[line];
            (
                params[e.param(2)] - params[e.param(0)],
                params[e.param(3)] - params[e.param(1)],
            )
        };
        match self {
            Constraint::PointDistance { p1, p2, .. } => {
                let (x1, y1) = p1.point_coords(params, entities);
                let (x2, y2) = p2.point_coords(params, entities);
                Some(((x1 - x2).powi(2) + (y1 - y2).powi(2)).sqrt())
            }
            Constraint::DistanceX { p1, p2, .. } => {
                Some(p1.point_coords(params, entities).0 - p2.point_coords(params, entities).0)
            }
            Constraint::DistanceY { p1, p2, .. } => {
                Some(p1.point_coords(params, entities).1 - p2.point_coords(params, entities).1)
            }
            Constraint::LineLength { line, .. } => {
                let (dx, dy) = direction(*line);
                Some((dx * dx + dy * dy).sqrt())
            }
            Constraint::Angle { l1, l2, .. } => {
                let (dx1, dy1) = direction(*l1);
                let (dx2, dy2) = direction(*l2);
                Some((dx1 * dy2 - dy1 * dx2).atan2(dx1 * dx2 + dy1 * dy2))
            }
            Constraint::Radius { circle, .. } => Some(params[entities[*circle].param(2)]),
            Constraint::Diameter { circle, .. } => Some(2.0 * params[entities[*circle].param(2)]),
            _ => None,
        }
    }
}

/// Reflexive AsRef so the solver's residual/Jacobian helpers accept both
/// `&[Constraint]` (full sketches) and `&[&Constraint]` (component
/// subsets) through one generic bound.
impl AsRef<Constraint> for Constraint {
    fn as_ref(&self) -> &Constraint {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::EntityKind;
    use crate::numerics::numerical_jacobian_into;

    /// The generic M1.2-style lockstep gate (solver roadmap M2.1): EVERY
    /// constraint variant's analytic Jacobian must match central finite
    /// differences of its residuals, across deterministic sample
    /// configurations. A variant added without derivatives — or with a
    /// derivative that drifts from its residual — fails here.
    #[test]
    fn analytic_jacobian_matches_finite_difference() {
        // Entity layout: points 0-2, lines 3-4, circles 5-6, arcs 7-8.
        let entities: Vec<Entity> = vec![
            Entity::new(EntityKind::Point, 0),   // params 0-1
            Entity::new(EntityKind::Point, 2),   // params 2-3
            Entity::new(EntityKind::Point, 4),   // params 4-5
            Entity::new(EntityKind::Line, 6),    // params 6-9
            Entity::new(EntityKind::Line, 10),   // params 10-13
            Entity::new(EntityKind::Circle, 14), // params 14-16
            Entity::new(EntityKind::Circle, 17), // params 17-19
            Entity::new(EntityKind::Arc, 20),    // params 20-24
            Entity::new(EntityKind::Arc, 25),    // params 25-29
        ];
        let n_params = 30;

        // Deterministic sample configurations (non-degenerate: no zero
        // lengths, no point-on-circle, no coincident points — the FD
        // comparison is meaningless at kinks).
        let samples: Vec<Vec<f64>> = (0..3)
            .map(|s| {
                let t = s as f64;
                vec![
                    1.0 + 0.3 * t,
                    2.0 - 0.2 * t, // p0
                    -1.0 + 0.4 * t,
                    3.0 + 0.1 * t, // p1
                    2.5 - 0.3 * t,
                    -2.0 + 0.5 * t, // p2
                    0.5 + 0.2 * t,
                    0.5,
                    4.0,
                    2.0 - 0.3 * t, // line 3
                    -3.0,
                    -1.0 + 0.4 * t,
                    1.0,
                    -2.5, // line 4
                    0.5 + 0.1 * t,
                    -1.0,
                    1.5, // circle 5
                    2.0,
                    2.0 + 0.2 * t,
                    2.5, // circle 6
                    1.0,
                    1.0,
                    2.0,
                    0.2,
                    2.4, // arc 7
                    -1.0 + 0.2 * t,
                    0.5,
                    3.0,
                    3.5,
                    5.9, // arc 8
                ]
            })
            .collect();

        let variants: Vec<(&str, Constraint)> = vec![
            (
                "Coincident",
                Constraint::Coincident(PointRef::Point(0), PointRef::Point(1)),
            ),
            (
                "Fixed",
                Constraint::Fixed {
                    point: PointRef::Point(2),
                    x: 1.0,
                    y: 1.0,
                },
            ),
            (
                "PointDistance",
                Constraint::PointDistance {
                    p1: PointRef::Point(0),
                    p2: PointRef::LineEnd(3),
                    distance: 4.0,
                },
            ),
            (
                "DistanceX",
                Constraint::DistanceX {
                    p1: PointRef::Center(5),
                    p2: PointRef::Point(1),
                    distance: 2.0,
                },
            ),
            (
                "DistanceY",
                Constraint::DistanceY {
                    p1: PointRef::LineStart(4),
                    p2: PointRef::Center(6),
                    distance: 3.0,
                },
            ),
            (
                "LineLength",
                Constraint::LineLength {
                    line: 3,
                    length: 4.0,
                },
            ),
            ("Horizontal", Constraint::Horizontal(3)),
            ("Vertical", Constraint::Vertical(4)),
            ("Parallel", Constraint::Parallel(3, 4)),
            ("Perpendicular", Constraint::Perpendicular(3, 4)),
            ("EqualLength", Constraint::EqualLength(3, 4)),
            (
                "Angle",
                Constraint::Angle {
                    l1: 3,
                    l2: 4,
                    angle_rad: 0.7,
                },
            ),
            ("EqualRadius", Constraint::EqualRadius(5, 6)),
            (
                "Radius",
                Constraint::Radius {
                    circle: 5,
                    radius: 1.5,
                },
            ),
            (
                "PointOnCircle",
                Constraint::PointOnCircle {
                    point: PointRef::Point(2),
                    circle: 5,
                },
            ),
            (
                "PointOnLine",
                Constraint::PointOnLine {
                    point: PointRef::Point(2),
                    line: 3,
                },
            ),
            ("Tangent", Constraint::Tangent { circle: 5, line: 3 }),
            (
                "CircleCircleTangent-ext",
                Constraint::CircleCircleTangent {
                    c1: 5,
                    c2: 6,
                    external: true,
                },
            ),
            (
                "CircleCircleTangent-int",
                Constraint::CircleCircleTangent {
                    c1: 5,
                    c2: 6,
                    external: false,
                },
            ),
            (
                "ArcArcTangent",
                Constraint::ArcArcTangent {
                    a1: 7,
                    a2: 8,
                    external: false,
                },
            ),
            (
                "Symmetric",
                Constraint::Symmetric {
                    p1: PointRef::Point(0),
                    p2: PointRef::Point(1),
                    line: 3,
                },
            ),
            (
                "SymmetricPoint",
                Constraint::SymmetricPoint {
                    p1: PointRef::Point(0),
                    p2: PointRef::Center(5),
                    about: PointRef::Point(1),
                },
            ),
            ("Concentric", Constraint::Concentric(5, 6)),
            (
                "Midpoint",
                Constraint::Midpoint {
                    point: PointRef::Point(2),
                    line: 4,
                },
            ),
            (
                "Diameter",
                Constraint::Diameter {
                    circle: 6,
                    diameter: 5.0,
                },
            ),
        ];

        for params in &samples {
            for (name, c) in &variants {
                let n_eq = c.equation_count();
                let mut analytic = vec![Vec::new(); n_eq];
                c.jacobian(params, &entities, &mut analytic);

                // Dense analytic for comparison.
                let mut ja = vec![vec![0.0_f64; n_params]; n_eq];
                for (r, entries) in analytic.iter().enumerate() {
                    assert!(
                        !entries.is_empty() || n_eq == 0,
                        "{name}: row {r} left empty by the analytic jacobian"
                    );
                    for &(p, v) in entries {
                        ja[r][p] += v;
                    }
                }

                // Finite differences of the residuals.
                let free: Vec<usize> = (0..n_params).collect();
                let mut eval = |p: &[f64], out: &mut [f64]| {
                    c.residuals(p, &entities, out);
                };
                let mut jfd = vec![vec![0.0_f64; n_params]; n_eq];
                numerical_jacobian_into(params, &free, n_eq, &mut eval, &mut jfd);

                for (r, (ra, rf)) in ja.iter().zip(jfd.iter()).enumerate() {
                    for (p, (va, vf)) in ra.iter().zip(rf.iter()).enumerate() {
                        assert!(
                            (va - vf).abs() <= 1.0e-5 * (1.0 + vf.abs()),
                            "{name}: analytic row {r} d/param {p} = {va} vs fd {vf} (sample params {params:?})"
                        );
                    }
                }
            }
        }
    }
}
