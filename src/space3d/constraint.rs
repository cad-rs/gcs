//! 3D geometric constraints and their residual equations.
//!
//! Each [`SpaceConstraint`] variant encodes one or more scalar equations f(x)=0.
//! The [`SpaceConstraint::residuals`] method evaluates those equations given the
//! current parameter vector, and [`SpaceConstraint::equation_count`] returns how
//! many equations the constraint contributes.

use crate::solver::SOLVER_NORM_TOL;
use crate::space3d::entity::{SpaceEntity, SpaceEntityId, SpacePointRef};

/// A geometric constraint between 3D space entities.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub enum SpaceConstraint {
    // ── Point constraints ────────────────────────────────────────────────────
    /// Two points occupy the same location in 3D. (3 equations)
    Coincident(SpacePointRef, SpacePointRef),

    /// A point is fixed at a specific (x, y, z) location. (3 equations)
    Fixed {
        point: SpacePointRef,
        x: f64,
        y: f64,
        z: f64,
    },

    // ── Distance ─────────────────────────────────────────────────────────────
    /// Euclidean distance between two points equals `distance`. (1 equation)
    PointDistance {
        p1: SpacePointRef,
        p2: SpacePointRef,
        distance: f64,
    },

    /// Directed distance between two points along a direction:
    /// `(p1 − p2)·dir̂ = distance`. (1 equation)
    ///
    /// The direction is normalized internally; the sign follows the
    /// argument order, matching the 2D `DistanceX`/`DistanceY` discipline.
    PointDistanceAlongDir {
        p1: SpacePointRef,
        p2: SpacePointRef,
        dx: f64,
        dy: f64,
        dz: f64,
        distance: f64,
    },

    /// Length of a SpaceLine entity equals `length`. (1 equation)
    LineLength { line: SpaceEntityId, length: f64 },

    // ── Point-on-entity constraints ──────────────────────────────────────────
    /// A point lies on a line (infinite extension). (2 equations)
    PointOnLine {
        point: SpacePointRef,
        line: SpaceEntityId,
    },

    /// A point lies on a plane. (1 equation)
    PointOnPlane {
        point: SpacePointRef,
        plane: SpaceEntityId,
    },

    /// A point lies on a sphere surface. (1 equation)
    PointOnSphere {
        point: SpacePointRef,
        sphere: SpaceEntityId,
    },

    // ── Line orientation ─────────────────────────────────────────────────────
    /// A SpaceLine is parallel to a given direction vector (dx, dy, dz). (2 equations)
    /// The direction is normalized internally.
    LineParallelToDir {
        line: SpaceEntityId,
        dx: f64,
        dy: f64,
        dz: f64,
    },

    /// A SpaceLine is parallel to another SpaceLine. (2 equations)
    LineParallelLine(SpaceEntityId, SpaceEntityId),

    /// A SpaceLine is perpendicular to another SpaceLine. (1 equation)
    LinePerpendicularLine(SpaceEntityId, SpaceEntityId),

    /// Two SpaceLines have equal length. (1 equation)
    EqualLineLength(SpaceEntityId, SpaceEntityId),

    // ── Plane constraints ────────────────────────────────────────────────────
    /// A plane has a fixed normal direction (nx, ny, nz). (2 equations)
    /// The normal is normalized internally.
    PlaneNormal {
        plane: SpaceEntityId,
        nx: f64,
        ny: f64,
        nz: f64,
    },

    /// Two planes are parallel (normals are parallel). (2 equations)
    PlaneParallel(SpaceEntityId, SpaceEntityId),

    /// Two planes are perpendicular (normals are orthogonal). (1 equation)
    PlanePerpendicular(SpaceEntityId, SpaceEntityId),

    /// A line lies in a plane (both endpoints on the plane). (2 equations)
    LineInPlane {
        line: SpaceEntityId,
        plane: SpaceEntityId,
    },

    /// A line is perpendicular to a plane (line direction parallel to plane normal). (2 equations)
    LinePerpendicularPlane {
        line: SpaceEntityId,
        plane: SpaceEntityId,
    },

    /// A line is parallel to a plane (line direction perpendicular to plane normal). (1 equation)
    LineParallelPlane {
        line: SpaceEntityId,
        plane: SpaceEntityId,
    },

    // ── Sphere constraints ───────────────────────────────────────────────────
    /// Sphere radius equals a fixed value. (1 equation)
    SphereRadius { sphere: SpaceEntityId, radius: f64 },

    /// Two spheres are tangent (distance between centers = sum/diff of radii). (1 equation)
    SphereTangent {
        s1: SpaceEntityId,
        s2: SpaceEntityId,
        external: bool,
    },

    /// Two spheres have equal radius. (1 equation)
    EqualSphereRadius(SpaceEntityId, SpaceEntityId),

    // ── Angle ────────────────────────────────────────────────────────────────
    /// Angle between two lines equals `angle_rad`. (1 equation)
    AngleLineLine {
        l1: SpaceEntityId,
        l2: SpaceEntityId,
        angle_rad: f64,
    },

    /// Angle between two planes equals `angle_rad`. (1 equation)
    AnglePlanePlane {
        p1: SpaceEntityId,
        p2: SpaceEntityId,
        angle_rad: f64,
    },

    /// Angle between a line and a plane equals `angle_rad`. (1 equation)
    AngleLinePlane {
        line: SpaceEntityId,
        plane: SpaceEntityId,
        angle_rad: f64,
    },

    // ── Distance (entity to entity) ──────────────────────────────────────────
    /// Distance from a point to a plane equals `distance`. (1 equation)
    PointPlaneDistance {
        point: SpacePointRef,
        plane: SpaceEntityId,
        distance: f64,
    },

    // ── Cylinder constraints ─────────────────────────────────────────────────
    /// Two Cylinder/Sphere entities share the same center. (3 equations)
    Concentric3D(SpaceEntityId, SpaceEntityId),

    /// Cylinder radius equals a fixed value. (1 equation)
    CylinderRadius {
        cylinder: SpaceEntityId,
        radius: f64,
    },

    /// A point lies on a cylinder surface. (1 equation)
    PointOnCylinder {
        point: SpacePointRef,
        cylinder: SpaceEntityId,
    },

    /// Two cylinders with parallel axes are tangent. (1 equation)
    CylinderTangent {
        c1: SpaceEntityId,
        c2: SpaceEntityId,
        external: bool,
    },

    /// A plane is tangent to a sphere. (1 equation)
    PlaneTangentToSphere {
        plane: SpaceEntityId,
        sphere: SpaceEntityId,
    },

    // ── Circle/arc constraints (3D-3) ────────────────────────────────────────
    /// A point lies on a circle: it is in the circle's plane and at the
    /// radius in that plane. (2 equations)
    ///
    /// Equation 1 is the signed plane distance; equation 2 is the squared
    /// in-plane radius `|u|² − (u·n)² − r²` (smooth everywhere, no square
    /// root, and exactly the squared 2D radial residual when the circle lies
    /// on a workplane).
    PointOnCircle {
        point: SpacePointRef,
        circle: SpaceEntityId,
    },

    /// A point lies on the supporting circle of an arc. (2 equations)
    ///
    /// The equations are [`SpaceConstraint::PointOnCircle`]'s; the arc span
    /// is enforced by the editable-model translation, which projects an
    /// out-of-span point onto the nearer arc endpoint before solving (the
    /// same stick-to-the-endpoint behavior a 2D point-on-arc constraint has
    /// when dragged past its span).
    PointOnArc {
        point: SpacePointRef,
        arc: SpaceEntityId,
    },

    /// A line is tangent to a circle/arc: the distance from the center to
    /// the infinite line equals the radius. (1 equation)
    ///
    /// In 3D the distance is the perpendicular distance to the line; when
    /// the line and the circle share the workplane this is exactly the 2D
    /// `dist(center, line) − r` equation.
    LineTangentCircle {
        line: SpaceEntityId,
        circle: SpaceEntityId,
    },

    /// Two circles/arcs are tangent (center distance = sum or absolute
    /// difference of radii). (1 equation)
    CircleCircleTangent {
        c1: SpaceEntityId,
        c2: SpaceEntityId,
        external: bool,
    },

    /// Two circles/arcs have equal radius. (1 equation)
    EqualRadius(SpaceEntityId, SpaceEntityId),

    /// Circle/arc radius equals a fixed value. (1 equation)
    CircleRadius { circle: SpaceEntityId, radius: f64 },

    /// A point is the midpoint of a line segment. (3 equations)
    Midpoint {
        point: SpacePointRef,
        line: SpaceEntityId,
    },

    /// Two points are mirror images about a line: their midpoint lies on the
    /// line (2 equations) and the segment joining them is perpendicular to
    /// the line (1 equation). (3 equations)
    Symmetric {
        p1: SpacePointRef,
        p2: SpacePointRef,
        line: SpaceEntityId,
    },

    /// A point coincides with one endpoint of an arc (the frame's start
    /// angle or its normalized end angle). (3 equations)
    PointAtArcEndpoint {
        point: SpacePointRef,
        arc: SpaceEntityId,
        at_start: bool,
    },

    // ── Spline constraints (3D-4) ──────────────────────────────────────────
    /// A point lies on a Bezier spline curve. (2 equations)
    ///
    /// The residual is the vector from the point to the nearest point on the
    /// curve, projected onto the plane perpendicular to the curve tangent at
    /// that point. The component along the tangent is free (the point can
    /// slide along the curve).
    PointOnSpline {
        point: SpacePointRef,
        spline: SpaceEntityId,
    },
    /// A point coincides with a spline endpoint (start or end control point).
    /// (3 equations)
    SplineEndpoint {
        point: SpacePointRef,
        spline: SpaceEntityId,
        at_start: bool,
    },
    /// A control point is the midpoint of two other control points on the
    /// same spline. (3 equations)
    SplineControlPointMidpoint {
        spline: SpaceEntityId,
        midpoint_index: usize,
        first_index: usize,
        second_index: usize,
    },
    /// Two control points on the same spline are mirror images about a line.
    /// (3 equations)
    SplineControlPointSymmetric {
        spline: SpaceEntityId,
        first_index: usize,
        second_index: usize,
        line: SpaceEntityId,
    },
}

impl SpaceConstraint {
    // ── Convenience constructors ─────────────────────────────────────────────

    pub fn fix_point(point: impl Into<SpacePointRef>, x: f64, y: f64, z: f64) -> Self {
        SpaceConstraint::Fixed {
            point: point.into(),
            x,
            y,
            z,
        }
    }

    pub fn point_distance(
        p1: impl Into<SpacePointRef>,
        p2: impl Into<SpacePointRef>,
        distance: f64,
    ) -> Self {
        SpaceConstraint::PointDistance {
            p1: p1.into(),
            p2: p2.into(),
            distance,
        }
    }

    pub fn coincident(p1: impl Into<SpacePointRef>, p2: impl Into<SpacePointRef>) -> Self {
        SpaceConstraint::Coincident(p1.into(), p2.into())
    }

    // ── Equation count ───────────────────────────────────────────────────────

    /// Number of scalar equations this constraint contributes.
    pub fn equation_count(&self) -> usize {
        match self {
            SpaceConstraint::Coincident(..)
            | SpaceConstraint::Fixed { .. }
            | SpaceConstraint::Concentric3D(..)
            | SpaceConstraint::Midpoint { .. }
            | SpaceConstraint::Symmetric { .. }
            | SpaceConstraint::PointAtArcEndpoint { .. }
            | SpaceConstraint::SplineEndpoint { .. }
            | SpaceConstraint::SplineControlPointMidpoint { .. }
            | SpaceConstraint::SplineControlPointSymmetric { .. } => 3,
            SpaceConstraint::PointOnCircle { .. } | SpaceConstraint::PointOnArc { .. } => 2,
            SpaceConstraint::PointOnSpline { .. } => 2,
            SpaceConstraint::PointOnLine { .. }
            | SpaceConstraint::PlaneNormal { .. }
            | SpaceConstraint::PlaneParallel(..)
            | SpaceConstraint::LineInPlane { .. }
            | SpaceConstraint::LinePerpendicularPlane { .. }
            | SpaceConstraint::LineParallelToDir { .. }
            | SpaceConstraint::LineParallelLine(..) => 2,
            _ => 1,
        }
    }

    // ── Residual evaluation ──────────────────────────────────────────────────

    /// Evaluate the constraint residuals f(x) into `out`.
    ///
    /// `out` must have length == `equation_count()`.
    /// Returns `false` if the constraint references an entity of the wrong kind.
    pub fn residuals(&self, params: &[f64], entities: &[SpaceEntity], out: &mut [f64]) -> bool {
        match self {
            // ── Coincident ───────────────────────────────────────────────────
            SpaceConstraint::Coincident(p1, p2) => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                out[0] = params[x1] - params[x2];
                out[1] = params[y1] - params[y2];
                out[2] = params[z1] - params[z2];
            }

            // ── Fixed ────────────────────────────────────────────────────────
            SpaceConstraint::Fixed { point, x, y, z } => {
                let (xi, yi, zi) = point.param_indices(entities);
                out[0] = params[xi] - x;
                out[1] = params[yi] - y;
                out[2] = params[zi] - z;
            }

            // ── PointDistance ────────────────────────────────────────────────
            SpaceConstraint::PointDistance { p1, p2, distance } => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                let dx = params[x1] - params[x2];
                let dy = params[y1] - params[y2];
                let dz = params[z1] - params[z2];
                out[0] = (dx * dx + dy * dy + dz * dz).sqrt() - distance;
            }

            // ── LineLength ───────────────────────────────────────────────────
            SpaceConstraint::LineLength { line, length } => {
                let e = &entities[*line];
                let dx = params[e.param(3)] - params[e.param(0)];
                let dy = params[e.param(4)] - params[e.param(1)];
                let dz = params[e.param(5)] - params[e.param(2)];
                out[0] = (dx * dx + dy * dy + dz * dz).sqrt() - length;
            }

            // ── PointDistanceAlongDir ────────────────────────────────────────
            SpaceConstraint::PointDistanceAlongDir {
                p1,
                p2,
                dx,
                dy,
                dz,
                distance,
            } => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                let len_d = (dx * dx + dy * dy + dz * dz).sqrt().max(SOLVER_NORM_TOL);
                let ux = dx / len_d;
                let uy = dy / len_d;
                let uz = dz / len_d;
                out[0] = (params[x1] - params[x2]) * ux
                    + (params[y1] - params[y2]) * uy
                    + (params[z1] - params[z2]) * uz
                    - distance;
            }

            // ── PointOnLine ──────────────────────────────────────────────────
            SpaceConstraint::PointOnLine { point, line } => {
                let (px, py, pz) = point.param_indices(entities);
                let el = &entities[*line];
                let x1 = params[el.param(0)];
                let y1 = params[el.param(1)];
                let z1 = params[el.param(2)];
                let x2 = params[el.param(3)];
                let y2 = params[el.param(4)];
                let z2 = params[el.param(5)];
                let ldx = x2 - x1;
                let ldy = y2 - y1;
                let ldz = z2 - z1;
                let len = (ldx * ldx + ldy * ldy + ldz * ldz)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                // Cross product of (point - p1) with line direction = 0.
                // The cross is orthogonal to the line direction, so the
                // component along that direction's dominant axis is implied
                // by the other two. Selecting the retained pair from the
                // direction (instead of a fixed y/z pair) keeps both rows
                // independent for every line orientation — a fixed pair
                // degenerates for an X-aligned line (row 0 vanishes and z
                // stays free).
                let dx = params[px] - x1;
                let dy = params[py] - y1;
                let dz = params[pz] - z1;
                let (i0, i1) = cross_pair(dominant_axis(ldx, ldy, ldz));
                let c = cross3((dx, dy, dz), (ldx, ldy, ldz));
                out[0] = comp(c, i0) / len;
                out[1] = comp(c, i1) / len;
            }

            // ── PointOnPlane ─────────────────────────────────────────────────
            SpaceConstraint::PointOnPlane { point, plane } => {
                let (px, py, pz) = point.param_indices(entities);
                let ep = &entities[*plane];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let d = params[ep.param(3)];
                out[0] = nx * params[px] + ny * params[py] + nz * params[pz] - d;
            }

            // ── PointOnSphere ────────────────────────────────────────────────
            SpaceConstraint::PointOnSphere { point, sphere } => {
                let (px, py, pz) = point.param_indices(entities);
                let es = &entities[*sphere];
                let cx = params[es.param(0)];
                let cy = params[es.param(1)];
                let cz = params[es.param(2)];
                let r = params[es.param(3)];
                let dx = params[px] - cx;
                let dy = params[py] - cy;
                let dz = params[pz] - cz;
                out[0] = (dx * dx + dy * dy + dz * dz).sqrt() - r;
            }

            // ── LineParallelToDir ────────────────────────────────────────────
            SpaceConstraint::LineParallelToDir { line, dx, dy, dz } => {
                let e = &entities[*line];
                let ldx = params[e.param(3)] - params[e.param(0)];
                let ldy = params[e.param(4)] - params[e.param(1)];
                let ldz = params[e.param(5)] - params[e.param(2)];
                let len_l = (ldx * ldx + ldy * ldy + ldz * ldz)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len_d = (dx * dx + dy * dy + dz * dz).sqrt().max(SOLVER_NORM_TOL);
                // Cross product = 0 (normalized). The cross is orthogonal to
                // the target direction, so the component along the target's
                // dominant axis is implied by the other two.
                let ndx = dx / len_d;
                let ndy = dy / len_d;
                let ndz = dz / len_d;
                let (i0, i1) = cross_pair(dominant_axis(*dx, *dy, *dz));
                let c = cross3((ldx, ldy, ldz), (ndx, ndy, ndz));
                out[0] = comp(c, i0) / len_l;
                out[1] = comp(c, i1) / len_l;
            }

            // ── LineParallelLine ─────────────────────────────────────────────
            SpaceConstraint::LineParallelLine(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(3)] - params[e1.param(0)];
                let dy1 = params[e1.param(4)] - params[e1.param(1)];
                let dz1 = params[e1.param(5)] - params[e1.param(2)];
                let dx2 = params[e2.param(3)] - params[e2.param(0)];
                let dy2 = params[e2.param(4)] - params[e2.param(1)];
                let dz2 = params[e2.param(5)] - params[e2.param(2)];
                let len1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                // Cross product = 0 (normalized). The cross is orthogonal to
                // the first line's direction: dropping the component along
                // that direction's dominant axis keeps two independent rows.
                let (i0, i1) = cross_pair(dominant_axis(dx1, dy1, dz1));
                let c = cross3((dx1, dy1, dz1), (dx2, dy2, dz2));
                out[0] = comp(c, i0) / (len1 * len2);
                out[1] = comp(c, i1) / (len1 * len2);
            }

            // ── LinePerpendicularLine ────────────────────────────────────────
            SpaceConstraint::LinePerpendicularLine(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(3)] - params[e1.param(0)];
                let dy1 = params[e1.param(4)] - params[e1.param(1)];
                let dz1 = params[e1.param(5)] - params[e1.param(2)];
                let dx2 = params[e2.param(3)] - params[e2.param(0)];
                let dy2 = params[e2.param(4)] - params[e2.param(1)];
                let dz2 = params[e2.param(5)] - params[e2.param(2)];
                let len1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                // Dot product = 0 (normalized)
                out[0] = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (len1 * len2);
            }

            // ── EqualLineLength ──────────────────────────────────────────────
            SpaceConstraint::EqualLineLength(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(3)] - params[e1.param(0)];
                let dy1 = params[e1.param(4)] - params[e1.param(1)];
                let dz1 = params[e1.param(5)] - params[e1.param(2)];
                let dx2 = params[e2.param(3)] - params[e2.param(0)];
                let dy2 = params[e2.param(4)] - params[e2.param(1)];
                let dz2 = params[e2.param(5)] - params[e2.param(2)];
                out[0] = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1) - (dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
            }

            // ── PlaneNormal ──────────────────────────────────────────────────
            SpaceConstraint::PlaneNormal { plane, nx, ny, nz } => {
                let e = &entities[*plane];
                let enx = params[e.param(0)];
                let eny = params[e.param(1)];
                let enz = params[e.param(2)];
                let len = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                let nnx = nx / len;
                let nny = ny / len;
                let nnz = nz / len;
                // Cross product of current normal and target normal = 0;
                // the target's dominant axis component is implied.
                let (i0, i1) = cross_pair(dominant_axis(*nx, *ny, *nz));
                let c = cross3((enx, eny, enz), (nnx, nny, nnz));
                out[0] = comp(c, i0);
                out[1] = comp(c, i1);
            }

            // ── PlaneParallel ────────────────────────────────────────────────
            SpaceConstraint::PlaneParallel(p1, p2) => {
                let e1 = &entities[*p1];
                let e2 = &entities[*p2];
                let nx1 = params[e1.param(0)];
                let ny1 = params[e1.param(1)];
                let nz1 = params[e1.param(2)];
                let nx2 = params[e2.param(0)];
                let ny2 = params[e2.param(1)];
                let nz2 = params[e2.param(2)];
                let len1 = (nx1 * nx1 + ny1 * ny1 + nz1 * nz1)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len2 = (nx2 * nx2 + ny2 * ny2 + nz2 * nz2)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                // Cross product = 0 (normalized); two independent components
                // are selected against the first normal's dominant axis.
                let (i0, i1) = cross_pair(dominant_axis(nx1, ny1, nz1));
                let c = cross3((nx1, ny1, nz1), (nx2, ny2, nz2));
                out[0] = comp(c, i0) / (len1 * len2);
                out[1] = comp(c, i1) / (len1 * len2);
            }

            // ── PlanePerpendicular ───────────────────────────────────────────
            SpaceConstraint::PlanePerpendicular(p1, p2) => {
                let e1 = &entities[*p1];
                let e2 = &entities[*p2];
                let nx1 = params[e1.param(0)];
                let ny1 = params[e1.param(1)];
                let nz1 = params[e1.param(2)];
                let nx2 = params[e2.param(0)];
                let ny2 = params[e2.param(1)];
                let nz2 = params[e2.param(2)];
                let len1 = (nx1 * nx1 + ny1 * ny1 + nz1 * nz1)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len2 = (nx2 * nx2 + ny2 * ny2 + nz2 * nz2)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                out[0] = (nx1 * nx2 + ny1 * ny2 + nz1 * nz2) / (len1 * len2);
            }

            // ── LineInPlane ──────────────────────────────────────────────────
            SpaceConstraint::LineInPlane { line, plane } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let d = params[ep.param(3)];
                // Both endpoints must satisfy plane equation
                let x1 = params[el.param(0)];
                let y1 = params[el.param(1)];
                let z1 = params[el.param(2)];
                let x2 = params[el.param(3)];
                let y2 = params[el.param(4)];
                let z2 = params[el.param(5)];
                out[0] = nx * x1 + ny * y1 + nz * z1 - d;
                out[1] = nx * x2 + ny * y2 + nz * z2 - d;
            }

            // ── LinePerpendicularPlane ───────────────────────────────────────
            SpaceConstraint::LinePerpendicularPlane { line, plane } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let ldx = params[el.param(3)] - params[el.param(0)];
                let ldy = params[el.param(4)] - params[el.param(1)];
                let ldz = params[el.param(5)] - params[el.param(2)];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let len_l = (ldx * ldx + ldy * ldy + ldz * ldz)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                // Cross product of line direction and plane normal = 0; the
                // normal's dominant axis component is implied.
                let (i0, i1) = cross_pair(dominant_axis(nx, ny, nz));
                let c = cross3((ldx, ldy, ldz), (nx, ny, nz));
                out[0] = comp(c, i0) / (len_l * len_n);
                out[1] = comp(c, i1) / (len_l * len_n);
            }

            // ── LineParallelPlane ────────────────────────────────────────────
            SpaceConstraint::LineParallelPlane { line, plane } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let ldx = params[el.param(3)] - params[el.param(0)];
                let ldy = params[el.param(4)] - params[el.param(1)];
                let ldz = params[el.param(5)] - params[el.param(2)];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let len_l = (ldx * ldx + ldy * ldy + ldz * ldz)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                // Dot product = 0 (line direction perpendicular to plane normal)
                out[0] = (ldx * nx + ldy * ny + ldz * nz) / (len_l * len_n);
            }

            // ── SphereRadius ─────────────────────────────────────────────────
            SpaceConstraint::SphereRadius { sphere, radius } => {
                let e = &entities[*sphere];
                out[0] = params[e.param(3)] - radius;
            }

            // ── SphereTangent ────────────────────────────────────────────────
            SpaceConstraint::SphereTangent { s1, s2, external } => {
                let e1 = &entities[*s1];
                let e2 = &entities[*s2];
                let cx1 = params[e1.param(0)];
                let cy1 = params[e1.param(1)];
                let cz1 = params[e1.param(2)];
                let r1 = params[e1.param(3)];
                let cx2 = params[e2.param(0)];
                let cy2 = params[e2.param(1)];
                let cz2 = params[e2.param(2)];
                let r2 = params[e2.param(3)];
                let dx = cx2 - cx1;
                let dy = cy2 - cy1;
                let dz = cz2 - cz1;
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if *external {
                    out[0] = dist - (r1 + r2);
                } else {
                    out[0] = dist - (r1 - r2).abs();
                }
            }

            // ── EqualSphereRadius ────────────────────────────────────────────
            SpaceConstraint::EqualSphereRadius(s1, s2) => {
                let e1 = &entities[*s1];
                let e2 = &entities[*s2];
                out[0] = params[e1.param(3)] - params[e2.param(3)];
            }

            // ── AngleLineLine ────────────────────────────────────────────────
            SpaceConstraint::AngleLineLine { l1, l2, angle_rad } => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(3)] - params[e1.param(0)];
                let dy1 = params[e1.param(4)] - params[e1.param(1)];
                let dz1 = params[e1.param(5)] - params[e1.param(2)];
                let dx2 = params[e2.param(3)] - params[e2.param(0)];
                let dy2 = params[e2.param(4)] - params[e2.param(1)];
                let dz2 = params[e2.param(5)] - params[e2.param(2)];
                let len1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let dot = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (len1 * len2);
                // cos(angle) - dot = 0
                out[0] = angle_rad.cos() - dot.clamp(-1.0, 1.0);
            }

            // ── AnglePlanePlane ──────────────────────────────────────────────
            SpaceConstraint::AnglePlanePlane { p1, p2, angle_rad } => {
                let e1 = &entities[*p1];
                let e2 = &entities[*p2];
                let nx1 = params[e1.param(0)];
                let ny1 = params[e1.param(1)];
                let nz1 = params[e1.param(2)];
                let nx2 = params[e2.param(0)];
                let ny2 = params[e2.param(1)];
                let nz2 = params[e2.param(2)];
                let len1 = (nx1 * nx1 + ny1 * ny1 + nz1 * nz1)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len2 = (nx2 * nx2 + ny2 * ny2 + nz2 * nz2)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let dot = (nx1 * nx2 + ny1 * ny2 + nz1 * nz2) / (len1 * len2);
                out[0] = angle_rad.cos() - dot.clamp(-1.0, 1.0);
            }

            // ── AngleLinePlane ───────────────────────────────────────────────
            SpaceConstraint::AngleLinePlane {
                line,
                plane,
                angle_rad,
            } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let ldx = params[el.param(3)] - params[el.param(0)];
                let ldy = params[el.param(4)] - params[el.param(1)];
                let ldz = params[el.param(5)] - params[el.param(2)];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let len_l = (ldx * ldx + ldy * ldy + ldz * ldz)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                // Angle between line and plane: sin(angle) = |dir·n| / (|dir|*|n|)
                let sin_val = (ldx * nx + ldy * ny + ldz * nz).abs() / (len_l * len_n);
                out[0] = angle_rad.sin() - sin_val.clamp(0.0, 1.0);
            }

            // ── PointPlaneDistance ───────────────────────────────────────────
            SpaceConstraint::PointPlaneDistance {
                point,
                plane,
                distance,
            } => {
                let (px, py, pz) = point.param_indices(entities);
                let ep = &entities[*plane];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let d = params[ep.param(3)];
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                let num = (nx * params[px] + ny * params[py] + nz * params[pz] - d).abs();
                out[0] = num / len_n - distance;
            }

            // ── Concentric3D ─────────────────────────────────────────────────
            SpaceConstraint::Concentric3D(e1, e2) => {
                let a = &entities[*e1];
                let b = &entities[*e2];
                out[0] = params[a.param(0)] - params[b.param(0)];
                out[1] = params[a.param(1)] - params[b.param(1)];
                out[2] = params[a.param(2)] - params[b.param(2)];
            }

            // ── CylinderRadius ───────────────────────────────────────────────
            SpaceConstraint::CylinderRadius { cylinder, radius } => {
                let e = &entities[*cylinder];
                out[0] = params[e.param(6)] - radius;
            }

            // ── PointOnCylinder ──────────────────────────────────────────────
            SpaceConstraint::PointOnCylinder { point, cylinder } => {
                let (px, py, pz) = point.param_indices(entities);
                let ec = &entities[*cylinder];
                let cx = params[ec.param(0)];
                let cy = params[ec.param(1)];
                let cz = params[ec.param(2)];
                let ax = params[ec.param(3)];
                let ay = params[ec.param(4)];
                let az = params[ec.param(5)];
                let r = params[ec.param(6)];
                let len_a = (ax * ax + ay * ay + az * az).sqrt().max(SOLVER_NORM_TOL);
                let nax = ax / len_a;
                let nay = ay / len_a;
                let naz = az / len_a;
                // Vector from axis point to point
                let vx = params[px] - cx;
                let vy = params[py] - cy;
                let vz = params[pz] - cz;
                // Project onto axis
                let proj = vx * nax + vy * nay + vz * naz;
                // Perpendicular component
                let perp_x = vx - proj * nax;
                let perp_y = vy - proj * nay;
                let perp_z = vz - proj * naz;
                let dist = (perp_x * perp_x + perp_y * perp_y + perp_z * perp_z).sqrt();
                out[0] = dist - r;
            }

            // ── CylinderTangent ──────────────────────────────────────────────
            SpaceConstraint::CylinderTangent { c1, c2, external } => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                let cx1 = params[e1.param(0)];
                let cy1 = params[e1.param(1)];
                let cz1 = params[e1.param(2)];
                let ax1 = params[e1.param(3)];
                let ay1 = params[e1.param(4)];
                let az1 = params[e1.param(5)];
                let r1 = params[e1.param(6)];
                let cx2 = params[e2.param(0)];
                let cy2 = params[e2.param(1)];
                let cz2 = params[e2.param(2)];
                let _ax2 = params[e2.param(3)];
                let _ay2 = params[e2.param(4)];
                let _az2 = params[e2.param(5)];
                let r2 = params[e2.param(6)];
                // Normalize axis 1
                let len_a1 = (ax1 * ax1 + ay1 * ay1 + az1 * az1)
                    .sqrt()
                    .max(SOLVER_NORM_TOL);
                let nax1 = ax1 / len_a1;
                let nay1 = ay1 / len_a1;
                let naz1 = az1 / len_a1;
                // Vector between axis points
                let dx = cx2 - cx1;
                let dy = cy2 - cy1;
                let dz = cz2 - cz1;
                // Perpendicular distance from c2 axis point to c1 axis line
                let proj = dx * nax1 + dy * nay1 + dz * naz1;
                let perp_x = dx - proj * nax1;
                let perp_y = dy - proj * nay1;
                let perp_z = dz - proj * naz1;
                let perp_dist = (perp_x * perp_x + perp_y * perp_y + perp_z * perp_z).sqrt();
                if *external {
                    out[0] = perp_dist - (r1 + r2);
                } else {
                    out[0] = perp_dist - (r1 - r2).abs();
                }
            }

            // ── PlaneTangentToSphere ─────────────────────────────────────────
            SpaceConstraint::PlaneTangentToSphere { plane, sphere } => {
                let ep = &entities[*plane];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let d = params[ep.param(3)];
                let es = &entities[*sphere];
                let cx = params[es.param(0)];
                let cy = params[es.param(1)];
                let cz = params[es.param(2)];
                let r = params[es.param(3)];
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                let dist = (nx * cx + ny * cy + nz * cz - d).abs() / len_n;
                out[0] = dist - r;
            }

            // ── PointOnCircle / PointOnArc ───────────────────────────────────
            SpaceConstraint::PointOnCircle { point, circle }
            | SpaceConstraint::PointOnArc { point, arc: circle } => {
                let (px, py, pz) = point.param_indices(entities);
                let e = &entities[*circle];
                let frame = e.curve_frame();
                let cx = params[e.param(0)];
                let cy = params[e.param(1)];
                let cz = params[e.param(2)];
                let r = params[e.param(3)];
                let ux = params[px] - cx;
                let uy = params[py] - cy;
                let uz = params[pz] - cz;
                let h = ux * frame.normal.x + uy * frame.normal.y + uz * frame.normal.z;
                // Squared in-plane radius: |u|^2 - (u.n)^2 - r^2. Smooth,
                // no square root, and exactly the squared 2D radial residual
                // when the point and the curve share the workplane.
                out[0] = h;
                out[1] = ux * ux + uy * uy + uz * uz - h * h - r * r;
            }

            // ── LineTangentCircle ────────────────────────────────────────────
            SpaceConstraint::LineTangentCircle { line, circle } => {
                let el = &entities[*line];
                let ec = &entities[*circle];
                let (dx, dy, dz) = line_dir(params, el);
                let qx = params[ec.param(0)] - params[el.param(0)];
                let qy = params[ec.param(1)] - params[el.param(1)];
                let qz = params[ec.param(2)] - params[el.param(2)];
                let wx = qy * dz - qz * dy;
                let wy = qz * dx - qx * dz;
                let wz = qx * dy - qy * dx;
                let (len, _) = guarded_norm(dx * dx + dy * dy + dz * dz);
                let (cross_len, _) = guarded_norm(wx * wx + wy * wy + wz * wz);
                out[0] = cross_len / len - params[ec.param(3)];
            }

            // ── CircleCircleTangent ──────────────────────────────────────────
            SpaceConstraint::CircleCircleTangent { c1, c2, external } => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                let dx = params[e2.param(0)] - params[e1.param(0)];
                let dy = params[e2.param(1)] - params[e1.param(1)];
                let dz = params[e2.param(2)] - params[e1.param(2)];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                let r1 = params[e1.param(3)];
                let r2 = params[e2.param(3)];
                if *external {
                    out[0] = dist - (r1 + r2);
                } else {
                    out[0] = dist - (r1 - r2).abs();
                }
            }

            // ── EqualRadius ──────────────────────────────────────────────────
            SpaceConstraint::EqualRadius(c1, c2) => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                out[0] = params[e1.param(3)] - params[e2.param(3)];
            }

            // ── CircleRadius ─────────────────────────────────────────────────
            SpaceConstraint::CircleRadius { circle, radius } => {
                let e = &entities[*circle];
                out[0] = params[e.param(3)] - radius;
            }

            // ── Midpoint ─────────────────────────────────────────────────────
            SpaceConstraint::Midpoint { point, line } => {
                let (px, py, pz) = point.param_indices(entities);
                let el = &entities[*line];
                out[0] = params[px] - 0.5 * (params[el.param(0)] + params[el.param(3)]);
                out[1] = params[py] - 0.5 * (params[el.param(1)] + params[el.param(4)]);
                out[2] = params[pz] - 0.5 * (params[el.param(2)] + params[el.param(5)]);
            }

            // ── Symmetric ────────────────────────────────────────────────────
            SpaceConstraint::Symmetric { p1, p2, line } => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                let el = &entities[*line];
                let lx1 = params[el.param(0)];
                let ly1 = params[el.param(1)];
                let lz1 = params[el.param(2)];
                let d = line_dir(params, el);
                let (len, _) = guarded_norm(d.0 * d.0 + d.1 * d.1 + d.2 * d.2);
                // Midpoint lies on the line: two independent cross rows.
                let x = (
                    0.5 * (params[x1] + params[x2]) - lx1,
                    0.5 * (params[y1] + params[y2]) - ly1,
                    0.5 * (params[z1] + params[z2]) - lz1,
                );
                let c = cross3(x, d);
                let (i0, i1) = cross_pair(dominant_axis(d.0, d.1, d.2));
                out[0] = comp(c, i0) / len;
                out[1] = comp(c, i1) / len;
                // The joining segment is perpendicular to the line.
                let sx = params[x2] - params[x1];
                let sy = params[y2] - params[y1];
                let sz = params[z2] - params[z1];
                out[2] = (sx * d.0 + sy * d.1 + sz * d.2) / len;
            }

            // ── PointAtArcEndpoint ───────────────────────────────────────────
            SpaceConstraint::PointAtArcEndpoint {
                point,
                arc,
                at_start,
            } => {
                let (px, py, pz) = point.param_indices(entities);
                let e = &entities[*arc];
                let frame = e.curve_frame();
                let angle = if *at_start {
                    frame.start_angle
                } else {
                    frame.normalized_end_angle()
                };
                let radial = frame.radial(angle);
                let r = params[e.param(3)];
                let ex = params[e.param(0)] + r * radial.x;
                let ey = params[e.param(1)] + r * radial.y;
                let ez = params[e.param(2)] + r * radial.z;
                out[0] = params[px] - ex;
                out[1] = params[py] - ey;
                out[2] = params[pz] - ez;
            }

            // ── Spline constraints (3D-4) ──────────────────────────────────
            SpaceConstraint::PointOnSpline { point, spline } => {
                let (px, py, pz) = point.param_indices(entities);
                let e = &entities[*spline];
                let n = e.control_point_count;
                if n < 2 {
                    out[0] = 0.0;
                    out[1] = 0.0;
                    return true;
                }
                let cp: Vec<(f64, f64, f64)> = (0..n)
                    .map(|i| {
                        (
                            params[e.param(i * 3)],
                            params[e.param(i * 3 + 1)],
                            params[e.param(i * 3 + 2)],
                        )
                    })
                    .collect();
                let p = (params[px], params[py], params[pz]);
                let (t, _dist) = nearest_on_bezier(&cp, p);
                let (cx, cy, cz) = bezier_eval(&cp, t);
                let (tx, ty, tz) = bezier_deriv(&cp, t);
                // Vector from curve point to the constraint point.
                let vx = p.0 - cx;
                let vy = p.1 - cy;
                let vz = p.2 - cz;
                // Project onto the plane perpendicular to the tangent.
                let t_len_sq = tx * tx + ty * ty + tz * tz;
                let t_len = t_len_sq.sqrt().max(SOLVER_NORM_TOL);
                let ux = tx / t_len;
                let uy = ty / t_len;
                let uz = tz / t_len;
                let dot = vx * ux + vy * uy + vz * uz;
                let rx = vx - dot * ux;
                let ry = vy - dot * uy;
                let rz = vz - dot * uz;
                // Select 2 independent components using the tangent's
                // dominant axis (the component along the dominant axis is
                // implied by the other two).
                let (i0, i1) = cross_pair(dominant_axis(tx, ty, tz));
                out[0] = comp((rx, ry, rz), i0);
                out[1] = comp((rx, ry, rz), i1);
            }

            SpaceConstraint::SplineEndpoint {
                point,
                spline,
                at_start,
            } => {
                let (px, py, pz) = point.param_indices(entities);
                let e = &entities[*spline];
                let idx = if *at_start {
                    0
                } else {
                    e.control_point_count.saturating_sub(1)
                };
                out[0] = params[px] - params[e.param(idx * 3)];
                out[1] = params[py] - params[e.param(idx * 3 + 1)];
                out[2] = params[pz] - params[e.param(idx * 3 + 2)];
            }

            SpaceConstraint::SplineControlPointMidpoint {
                spline,
                midpoint_index,
                first_index,
                second_index,
            } => {
                let e = &entities[*spline];
                let mi = *midpoint_index * 3;
                let fi = *first_index * 3;
                let si = *second_index * 3;
                out[0] = params[e.param(mi)] - 0.5 * (params[e.param(fi)] + params[e.param(si)]);
                out[1] = params[e.param(mi + 1)]
                    - 0.5 * (params[e.param(fi + 1)] + params[e.param(si + 1)]);
                out[2] = params[e.param(mi + 2)]
                    - 0.5 * (params[e.param(fi + 2)] + params[e.param(si + 2)]);
            }

            SpaceConstraint::SplineControlPointSymmetric {
                spline,
                first_index,
                second_index,
                line,
            } => {
                let e = &entities[*spline];
                let fi = *first_index * 3;
                let si = *second_index * 3;
                let el = &entities[*line];
                let lx1 = params[el.param(0)];
                let ly1 = params[el.param(1)];
                let lz1 = params[el.param(2)];
                let d = line_dir(params, el);
                let (len, _) = guarded_norm(d.0 * d.0 + d.1 * d.1 + d.2 * d.2);
                // Midpoint lies on the line: two independent cross rows.
                let x = (
                    0.5 * (params[e.param(fi)] + params[e.param(si)]) - lx1,
                    0.5 * (params[e.param(fi + 1)] + params[e.param(si + 1)]) - ly1,
                    0.5 * (params[e.param(fi + 2)] + params[e.param(si + 2)]) - lz1,
                );
                let c = cross3(x, d);
                let (i0, i1) = cross_pair(dominant_axis(d.0, d.1, d.2));
                out[0] = comp(c, i0) / len;
                out[1] = comp(c, i1) / len;
                // The joining segment is perpendicular to the line.
                let sx = params[e.param(si)] - params[e.param(fi)];
                let sy = params[e.param(si + 1)] - params[e.param(fi + 1)];
                let sz = params[e.param(si + 2)] - params[e.param(fi + 2)];
                out[2] = (sx * d.0 + sy * d.1 + sz * d.2) / len;
            }
        }
        true
    }

    // ── Analytic Jacobian ────────────────────────────────────────────────────

    /// Analytic Jacobian of this constraint's equations — the companion of
    /// [`SpaceConstraint::residuals`], kept in lockstep with it by the
    /// generic analytic-versus-finite-difference test in this module.
    ///
    /// `out` must have length == [`SpaceConstraint::equation_count`];
    /// `out[r]` receives the `(absolute parameter index, dr/dparam)` pairs
    /// for equation row `r`. Entries are pushed for every parameter the
    /// constraint structurally touches — even where the derivative is zero —
    /// so the constraint-graph decomposition sees the complete touch set.
    /// Parameter indices held by the solver's fixed mask are filtered by the
    /// caller.
    pub fn jacobian(
        &self,
        params: &[f64],
        entities: &[SpaceEntity],
        out: &mut [Vec<(usize, f64)>],
    ) {
        // Push (param, derivative) entries; local shorthand for readability.
        macro_rules! row {
            ($r:expr, $($p:expr => $d:expr),+ $(,)?) => {{
                $( $r.push(($p, $d)); )+
            }};
        }
        match self {
            // ── Coincident ───────────────────────────────────────────────────
            SpaceConstraint::Coincident(p1, p2) => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                row!(out[0], x1 => 1.0, x2 => -1.0);
                row!(out[1], y1 => 1.0, y2 => -1.0);
                row!(out[2], z1 => 1.0, z2 => -1.0);
            }

            // ── Fixed ────────────────────────────────────────────────────────
            SpaceConstraint::Fixed { point, .. } => {
                let (xi, yi, zi) = point.param_indices(entities);
                row!(out[0], xi => 1.0);
                row!(out[1], yi => 1.0);
                row!(out[2], zi => 1.0);
            }

            // ── PointDistance ────────────────────────────────────────────────
            SpaceConstraint::PointDistance { p1, p2, .. } => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                let dx = params[x1] - params[x2];
                let dy = params[y1] - params[y2];
                let dz = params[z1] - params[z2];
                let (_, dl) = guarded_norm(dx * dx + dy * dy + dz * dz);
                let (ix, iy, iz) = (2.0 * dx * dl, 2.0 * dy * dl, 2.0 * dz * dl);
                row!(out[0],
                    x1 => ix, y1 => iy, z1 => iz,
                    x2 => -ix, y2 => -iy, z2 => -iz,
                );
            }

            // ── LineLength ───────────────────────────────────────────────────
            SpaceConstraint::LineLength { line, .. } => {
                let e = &entities[*line];
                let (p0, p1, p2) = (e.param(0), e.param(1), e.param(2));
                let (p3, p4, p5) = (e.param(3), e.param(4), e.param(5));
                let dx = params[p3] - params[p0];
                let dy = params[p4] - params[p1];
                let dz = params[p5] - params[p2];
                let (_, dl) = guarded_norm(dx * dx + dy * dy + dz * dz);
                let (ix, iy, iz) = (2.0 * dx * dl, 2.0 * dy * dl, 2.0 * dz * dl);
                row!(out[0],
                    p0 => -ix, p1 => -iy, p2 => -iz,
                    p3 => ix, p4 => iy, p5 => iz,
                );
            }

            // ── PointDistanceAlongDir ────────────────────────────────────────
            SpaceConstraint::PointDistanceAlongDir {
                p1, p2, dx, dy, dz, ..
            } => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                let len_d = (dx * dx + dy * dy + dz * dz).sqrt().max(SOLVER_NORM_TOL);
                let ux = dx / len_d;
                let uy = dy / len_d;
                let uz = dz / len_d;
                row!(out[0],
                    x1 => ux, y1 => uy, z1 => uz,
                    x2 => -ux, y2 => -uy, z2 => -uz,
                );
            }

            // ── PointOnLine ──────────────────────────────────────────────────
            SpaceConstraint::PointOnLine { point, line } => {
                let (px, py, pz) = point.param_indices(entities);
                let el = &entities[*line];
                let (lx1, ly1, lz1) = (el.param(0), el.param(1), el.param(2));
                let (lx2, ly2, lz2) = (el.param(3), el.param(4), el.param(5));
                let d = (
                    params[lx2] - params[lx1],
                    params[ly2] - params[ly1],
                    params[lz2] - params[lz1],
                );
                let (l, dl) = guarded_norm(d.0 * d.0 + d.1 * d.1 + d.2 * d.2);
                let x = (
                    params[px] - params[lx1],
                    params[py] - params[ly1],
                    params[pz] - params[lz1],
                );
                // r_i = (x x d)_i / |d|. The retained component pair follows
                // the residual selection: the cross is orthogonal to d, so
                // the dominant-axis component is implied (see `residuals`).
                let (i0, i1) = cross_pair(dominant_axis(d.0, d.1, d.2));
                let c = cross3(x, d);
                for (row_idx, i) in [(0usize, i0), (1usize, i1)] {
                    // d(c_i)/d(point) = e_k x d; d(c_i)/d(p2) = x x e_k.
                    let dc_point = cross_u_deriv(i, d);
                    let dc_end = cross_v_deriv(i, x);
                    let r = comp(c, i) / l;
                    let q = |dc: f64, dk: f64| dc / l - r * (2.0 * dl * dk / l);
                    row!(out[row_idx],
                        px => q(dc_point.0, 0.0),
                        py => q(dc_point.1, 0.0),
                        pz => q(dc_point.2, 0.0),
                        lx1 => q(-dc_point.0 - dc_end.0, -d.0),
                        ly1 => q(-dc_point.1 - dc_end.1, -d.1),
                        lz1 => q(-dc_point.2 - dc_end.2, -d.2),
                        lx2 => q(dc_end.0, d.0),
                        ly2 => q(dc_end.1, d.1),
                        lz2 => q(dc_end.2, d.2),
                    );
                }
            }

            // ── PointOnPlane ─────────────────────────────────────────────────
            SpaceConstraint::PointOnPlane { point, plane } => {
                let (px, py, pz) = point.param_indices(entities);
                let ep = &entities[*plane];
                let (nx, ny, nz) = (
                    params[ep.param(0)],
                    params[ep.param(1)],
                    params[ep.param(2)],
                );
                row!(out[0],
                    px => nx,
                    py => ny,
                    pz => nz,
                    ep.param(0) => params[px],
                    ep.param(1) => params[py],
                    ep.param(2) => params[pz],
                    ep.param(3) => -1.0,
                );
            }

            // ── PointOnSphere ────────────────────────────────────────────────
            SpaceConstraint::PointOnSphere { point, sphere } => {
                let (px, py, pz) = point.param_indices(entities);
                let es = &entities[*sphere];
                let cx = params[es.param(0)];
                let cy = params[es.param(1)];
                let cz = params[es.param(2)];
                let dx = params[px] - cx;
                let dy = params[py] - cy;
                let dz = params[pz] - cz;
                let (_, dl) = guarded_norm(dx * dx + dy * dy + dz * dz);
                let (ix, iy, iz) = (2.0 * dx * dl, 2.0 * dy * dl, 2.0 * dz * dl);
                row!(out[0],
                    px => ix, py => iy, pz => iz,
                    es.param(0) => -ix, es.param(1) => -iy, es.param(2) => -iz,
                    es.param(3) => -1.0,
                );
            }

            // ── LineParallelToDir ────────────────────────────────────────────
            SpaceConstraint::LineParallelToDir { line, dx, dy, dz } => {
                let e = &entities[*line];
                let (p0, p1, p2) = (e.param(0), e.param(1), e.param(2));
                let (p3, p4, p5) = (e.param(3), e.param(4), e.param(5));
                let u = line_dir(params, e);
                let (l, dl) = guarded_norm(u.0 * u.0 + u.1 * u.1 + u.2 * u.2);
                let len_d = (dx * dx + dy * dy + dz * dz).sqrt().max(SOLVER_NORM_TOL);
                let v = (dx / len_d, dy / len_d, dz / len_d);
                let (i0, i1) = cross_pair(dominant_axis(*dx, *dy, *dz));
                let c = cross3(u, v);
                for (row_idx, i) in [(0usize, i0), (1usize, i1)] {
                    let du = cross_u_deriv(i, v);
                    let r = comp(c, i) / l;
                    let q = |dc: f64, uk: f64| dc / l - r * (2.0 * dl * uk / l);
                    row!(out[row_idx],
                        p0 => -q(du.0, u.0),
                        p1 => -q(du.1, u.1),
                        p2 => -q(du.2, u.2),
                        p3 => q(du.0, u.0),
                        p4 => q(du.1, u.1),
                        p5 => q(du.2, u.2),
                    );
                }
            }

            // ── LineParallelLine ─────────────────────────────────────────────
            SpaceConstraint::LineParallelLine(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let u = line_dir(params, e1);
                let v = line_dir(params, e2);
                let (len1, dl1) = guarded_norm(u.0 * u.0 + u.1 * u.1 + u.2 * u.2);
                let (len2, dl2) = guarded_norm(v.0 * v.0 + v.1 * v.1 + v.2 * v.2);
                let d = len1 * len2;
                let (i0, i1) = cross_pair(dominant_axis(u.0, u.1, u.2));
                let c = cross3(u, v);
                for (row_idx, i) in [(0usize, i0), (1usize, i1)] {
                    let du = cross_u_deriv(i, v);
                    let dv = cross_v_deriv(i, u);
                    let r = comp(c, i) / d;
                    let qu = |dc: f64, uk: f64| dc / d - r * (2.0 * dl1 * uk / len1);
                    let qv = |dc: f64, vk: f64| dc / d - r * (2.0 * dl2 * vk / len2);
                    row!(out[row_idx],
                        e1.param(0) => -qu(du.0, u.0),
                        e1.param(1) => -qu(du.1, u.1),
                        e1.param(2) => -qu(du.2, u.2),
                        e1.param(3) => qu(du.0, u.0),
                        e1.param(4) => qu(du.1, u.1),
                        e1.param(5) => qu(du.2, u.2),
                        e2.param(0) => -qv(dv.0, v.0),
                        e2.param(1) => -qv(dv.1, v.1),
                        e2.param(2) => -qv(dv.2, v.2),
                        e2.param(3) => qv(dv.0, v.0),
                        e2.param(4) => qv(dv.1, v.1),
                        e2.param(5) => qv(dv.2, v.2),
                    );
                }
            }

            // ── LinePerpendicularLine ────────────────────────────────────────
            SpaceConstraint::LinePerpendicularLine(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let (dx1, dy1, dz1) = line_dir(params, e1);
                let (dx2, dy2, dz2) = line_dir(params, e2);
                let (l1, dl1) = guarded_norm(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
                let (l2, dl2) = guarded_norm(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                let d = l1 * l2;
                let r = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / d;
                let q = |dk_other: f64, dk: f64, l: f64, dl: f64| {
                    dk_other / d - r * (2.0 * dl * dk / l)
                };
                row!(out[0],
                    e1.param(0) => -q(dx2, dx1, l1, dl1),
                    e1.param(1) => -q(dy2, dy1, l1, dl1),
                    e1.param(2) => -q(dz2, dz1, l1, dl1),
                    e1.param(3) => q(dx2, dx1, l1, dl1),
                    e1.param(4) => q(dy2, dy1, l1, dl1),
                    e1.param(5) => q(dz2, dz1, l1, dl1),
                    e2.param(0) => -q(dx1, dx2, l2, dl2),
                    e2.param(1) => -q(dy1, dy2, l2, dl2),
                    e2.param(2) => -q(dz1, dz2, l2, dl2),
                    e2.param(3) => q(dx1, dx2, l2, dl2),
                    e2.param(4) => q(dy1, dy2, l2, dl2),
                    e2.param(5) => q(dz1, dz2, l2, dl2),
                );
            }

            // ── EqualLineLength ──────────────────────────────────────────────
            SpaceConstraint::EqualLineLength(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let (dx1, dy1, dz1) = line_dir(params, e1);
                let (dx2, dy2, dz2) = line_dir(params, e2);
                row!(out[0],
                    e1.param(0) => -2.0 * dx1,
                    e1.param(1) => -2.0 * dy1,
                    e1.param(2) => -2.0 * dz1,
                    e1.param(3) => 2.0 * dx1,
                    e1.param(4) => 2.0 * dy1,
                    e1.param(5) => 2.0 * dz1,
                    e2.param(0) => 2.0 * dx2,
                    e2.param(1) => 2.0 * dy2,
                    e2.param(2) => 2.0 * dz2,
                    e2.param(3) => -2.0 * dx2,
                    e2.param(4) => -2.0 * dy2,
                    e2.param(5) => -2.0 * dz2,
                );
            }

            // ── PlaneNormal ──────────────────────────────────────────────────
            SpaceConstraint::PlaneNormal { plane, nx, ny, nz } => {
                let e = &entities[*plane];
                let len = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                let v = (nx / len, ny / len, nz / len);
                let (i0, i1) = cross_pair(dominant_axis(*nx, *ny, *nz));
                for (row_idx, i) in [(0usize, i0), (1usize, i1)] {
                    let du = cross_u_deriv(i, v);
                    row!(out[row_idx],
                        e.param(0) => du.0,
                        e.param(1) => du.1,
                        e.param(2) => du.2,
                    );
                }
            }

            // ── PlaneParallel ────────────────────────────────────────────────
            SpaceConstraint::PlaneParallel(p1, p2) => {
                let e1 = &entities[*p1];
                let e2 = &entities[*p2];
                let u = plane_normal(params, e1);
                let v = plane_normal(params, e2);
                let (l1, dl1) = guarded_norm(u.0 * u.0 + u.1 * u.1 + u.2 * u.2);
                let (l2, dl2) = guarded_norm(v.0 * v.0 + v.1 * v.1 + v.2 * v.2);
                let d = l1 * l2;
                let (i0, i1) = cross_pair(dominant_axis(u.0, u.1, u.2));
                let c = cross3(u, v);
                for (row_idx, i) in [(0usize, i0), (1usize, i1)] {
                    let du = cross_u_deriv(i, v);
                    let dv = cross_v_deriv(i, u);
                    let r = comp(c, i) / d;
                    let qu = |dc: f64, uk: f64| dc / d - r * (2.0 * dl1 * uk / l1);
                    let qv = |dc: f64, vk: f64| dc / d - r * (2.0 * dl2 * vk / l2);
                    row!(out[row_idx],
                        e1.param(0) => qu(du.0, u.0),
                        e1.param(1) => qu(du.1, u.1),
                        e1.param(2) => qu(du.2, u.2),
                        e2.param(0) => qv(dv.0, v.0),
                        e2.param(1) => qv(dv.1, v.1),
                        e2.param(2) => qv(dv.2, v.2),
                    );
                }
            }

            // ── PlanePerpendicular ───────────────────────────────────────────
            SpaceConstraint::PlanePerpendicular(p1, p2) => {
                let e1 = &entities[*p1];
                let e2 = &entities[*p2];
                let (nx1, ny1, nz1) = plane_normal(params, e1);
                let (nx2, ny2, nz2) = plane_normal(params, e2);
                let (l1, dl1) = guarded_norm(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);
                let (l2, dl2) = guarded_norm(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                let d = l1 * l2;
                let r = (nx1 * nx2 + ny1 * ny2 + nz1 * nz2) / d;
                let q = |dk_other: f64, dk: f64, l: f64, dl: f64| {
                    dk_other / d - r * (2.0 * dl * dk / l)
                };
                row!(out[0],
                    e1.param(0) => q(nx2, nx1, l1, dl1),
                    e1.param(1) => q(ny2, ny1, l1, dl1),
                    e1.param(2) => q(nz2, nz1, l1, dl1),
                    e2.param(0) => q(nx1, nx2, l2, dl2),
                    e2.param(1) => q(ny1, ny2, l2, dl2),
                    e2.param(2) => q(nz1, nz2, l2, dl2),
                );
            }

            // ── LineInPlane ──────────────────────────────────────────────────
            SpaceConstraint::LineInPlane { line, plane } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let (nx, ny, nz) = plane_normal(params, ep);
                let (x1, y1, z1) = (
                    params[el.param(0)],
                    params[el.param(1)],
                    params[el.param(2)],
                );
                let (x2, y2, z2) = (
                    params[el.param(3)],
                    params[el.param(4)],
                    params[el.param(5)],
                );
                row!(out[0],
                    el.param(0) => nx,
                    el.param(1) => ny,
                    el.param(2) => nz,
                    ep.param(0) => x1,
                    ep.param(1) => y1,
                    ep.param(2) => z1,
                    ep.param(3) => -1.0,
                );
                row!(out[1],
                    el.param(3) => nx,
                    el.param(4) => ny,
                    el.param(5) => nz,
                    ep.param(0) => x2,
                    ep.param(1) => y2,
                    ep.param(2) => z2,
                    ep.param(3) => -1.0,
                );
            }

            // ── LinePerpendicularPlane ───────────────────────────────────────
            SpaceConstraint::LinePerpendicularPlane { line, plane } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let u = line_dir(params, el);
                let v = plane_normal(params, ep);
                let (ll, dll) = guarded_norm(u.0 * u.0 + u.1 * u.1 + u.2 * u.2);
                let (ln, dln) = guarded_norm(v.0 * v.0 + v.1 * v.1 + v.2 * v.2);
                let d = ll * ln;
                let (i0, i1) = cross_pair(dominant_axis(v.0, v.1, v.2));
                let c = cross3(u, v);
                for (row_idx, i) in [(0usize, i0), (1usize, i1)] {
                    let du = cross_u_deriv(i, v);
                    let dv = cross_v_deriv(i, u);
                    let r = comp(c, i) / d;
                    let ql = |dc: f64, uk: f64| dc / d - r * (2.0 * dll * uk / ll);
                    let qn = |dc: f64, nk: f64| dc / d - r * (2.0 * dln * nk / ln);
                    row!(out[row_idx],
                        el.param(0) => -ql(du.0, u.0),
                        el.param(1) => -ql(du.1, u.1),
                        el.param(2) => -ql(du.2, u.2),
                        el.param(3) => ql(du.0, u.0),
                        el.param(4) => ql(du.1, u.1),
                        el.param(5) => ql(du.2, u.2),
                        ep.param(0) => qn(dv.0, v.0),
                        ep.param(1) => qn(dv.1, v.1),
                        ep.param(2) => qn(dv.2, v.2),
                    );
                }
            }

            // ── LineParallelPlane ────────────────────────────────────────────
            SpaceConstraint::LineParallelPlane { line, plane } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let (ldx, ldy, ldz) = line_dir(params, el);
                let (nx, ny, nz) = plane_normal(params, ep);
                let (ll, dll) = guarded_norm(ldx * ldx + ldy * ldy + ldz * ldz);
                let (ln, dln) = guarded_norm(nx * nx + ny * ny + nz * nz);
                let d = ll * ln;
                let r = (ldx * nx + ldy * ny + ldz * nz) / d;
                let ql = |nk: f64, dk: f64| nk / d - r * (2.0 * dll * dk / ll);
                let qn = |ldk: f64, nk: f64| ldk / d - r * (2.0 * dln * nk / ln);
                row!(out[0],
                    el.param(0) => -ql(nx, ldx),
                    el.param(1) => -ql(ny, ldy),
                    el.param(2) => -ql(nz, ldz),
                    el.param(3) => ql(nx, ldx),
                    el.param(4) => ql(ny, ldy),
                    el.param(5) => ql(nz, ldz),
                    ep.param(0) => qn(ldx, nx),
                    ep.param(1) => qn(ldy, ny),
                    ep.param(2) => qn(ldz, nz),
                );
            }

            // ── SphereRadius ─────────────────────────────────────────────────
            SpaceConstraint::SphereRadius { sphere, .. } => {
                let e = &entities[*sphere];
                row!(out[0], e.param(3) => 1.0);
            }

            // ── SphereTangent ────────────────────────────────────────────────
            SpaceConstraint::SphereTangent { s1, s2, external } => {
                let e1 = &entities[*s1];
                let e2 = &entities[*s2];
                let dx = params[e2.param(0)] - params[e1.param(0)];
                let dy = params[e2.param(1)] - params[e1.param(1)];
                let dz = params[e2.param(2)] - params[e1.param(2)];
                let (_, dd) = guarded_norm(dx * dx + dy * dy + dz * dz);
                let (ix, iy, iz) = (2.0 * dx * dd, 2.0 * dy * dd, 2.0 * dz * dd);
                let (r1, r2) = (
                    if *external {
                        -1.0
                    } else {
                        -(params[e1.param(3)] - params[e2.param(3)]).signum()
                    },
                    if *external {
                        -1.0
                    } else {
                        (params[e1.param(3)] - params[e2.param(3)]).signum()
                    },
                );
                row!(out[0],
                    e1.param(0) => -ix,
                    e1.param(1) => -iy,
                    e1.param(2) => -iz,
                    e1.param(3) => r1,
                    e2.param(0) => ix,
                    e2.param(1) => iy,
                    e2.param(2) => iz,
                    e2.param(3) => r2,
                );
            }

            // ── EqualSphereRadius ────────────────────────────────────────────
            SpaceConstraint::EqualSphereRadius(s1, s2) => {
                let e1 = &entities[*s1];
                let e2 = &entities[*s2];
                row!(out[0], e1.param(3) => 1.0, e2.param(3) => -1.0);
            }

            // ── AngleLineLine ────────────────────────────────────────────────
            SpaceConstraint::AngleLineLine { l1, l2, angle_rad } => {
                let _ = angle_rad;
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let (dx1, dy1, dz1) = line_dir(params, e1);
                let (dx2, dy2, dz2) = line_dir(params, e2);
                let (ll1, dl1) = guarded_norm(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
                let (ll2, dl2) = guarded_norm(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
                let d = ll1 * ll2;
                let v = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / d;
                if v.abs() >= 1.0 {
                    // clamp derivative is zero (the residual is saturated).
                    row!(out[0],
                        e1.param(0) => 0.0, e1.param(1) => 0.0, e1.param(2) => 0.0,
                        e1.param(3) => 0.0, e1.param(4) => 0.0, e1.param(5) => 0.0,
                        e2.param(0) => 0.0, e2.param(1) => 0.0, e2.param(2) => 0.0,
                        e2.param(3) => 0.0, e2.param(4) => 0.0, e2.param(5) => 0.0,
                    );
                } else {
                    let q = |dk_other: f64, dk: f64, l: f64, dl: f64| {
                        dk_other / d - v * (2.0 * dl * dk / l)
                    };
                    row!(out[0],
                        e1.param(0) => q(dx2, dx1, ll1, dl1),
                        e1.param(1) => q(dy2, dy1, ll1, dl1),
                        e1.param(2) => q(dz2, dz1, ll1, dl1),
                        e1.param(3) => -q(dx2, dx1, ll1, dl1),
                        e1.param(4) => -q(dy2, dy1, ll1, dl1),
                        e1.param(5) => -q(dz2, dz1, ll1, dl1),
                        e2.param(0) => q(dx1, dx2, ll2, dl2),
                        e2.param(1) => q(dy1, dy2, ll2, dl2),
                        e2.param(2) => q(dz1, dz2, ll2, dl2),
                        e2.param(3) => -q(dx1, dx2, ll2, dl2),
                        e2.param(4) => -q(dy1, dy2, ll2, dl2),
                        e2.param(5) => -q(dz1, dz2, ll2, dl2),
                    );
                }
            }

            // ── AnglePlanePlane ──────────────────────────────────────────────
            SpaceConstraint::AnglePlanePlane { p1, p2, angle_rad } => {
                let _ = angle_rad;
                let e1 = &entities[*p1];
                let e2 = &entities[*p2];
                let (nx1, ny1, nz1) = plane_normal(params, e1);
                let (nx2, ny2, nz2) = plane_normal(params, e2);
                let (l1, dl1) = guarded_norm(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);
                let (l2, dl2) = guarded_norm(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                let d = l1 * l2;
                let v = (nx1 * nx2 + ny1 * ny2 + nz1 * nz2) / d;
                if v.abs() >= 1.0 {
                    row!(out[0],
                        e1.param(0) => 0.0, e1.param(1) => 0.0, e1.param(2) => 0.0,
                        e2.param(0) => 0.0, e2.param(1) => 0.0, e2.param(2) => 0.0,
                    );
                } else {
                    let q = |dk_other: f64, dk: f64, l: f64, dl: f64| {
                        dk_other / d - v * (2.0 * dl * dk / l)
                    };
                    row!(out[0],
                        e1.param(0) => -q(nx2, nx1, l1, dl1),
                        e1.param(1) => -q(ny2, ny1, l1, dl1),
                        e1.param(2) => -q(nz2, nz1, l1, dl1),
                        e2.param(0) => -q(nx1, nx2, l2, dl2),
                        e2.param(1) => -q(ny1, ny2, l2, dl2),
                        e2.param(2) => -q(nz1, nz2, l2, dl2),
                    );
                }
            }

            // ── AngleLinePlane ───────────────────────────────────────────────
            SpaceConstraint::AngleLinePlane {
                line,
                plane,
                angle_rad,
            } => {
                let _ = angle_rad;
                let el = &entities[*line];
                let ep = &entities[*plane];
                let (ldx, ldy, ldz) = line_dir(params, el);
                let (nx, ny, nz) = plane_normal(params, ep);
                let (ll, dll) = guarded_norm(ldx * ldx + ldy * ldy + ldz * ldz);
                let (ln, dln) = guarded_norm(nx * nx + ny * ny + nz * nz);
                let d = ll * ln;
                let v = (ldx * nx + ldy * ny + ldz * nz) / d;
                if v.abs() >= 1.0 {
                    row!(out[0],
                        el.param(0) => 0.0, el.param(1) => 0.0, el.param(2) => 0.0,
                        el.param(3) => 0.0, el.param(4) => 0.0, el.param(5) => 0.0,
                        ep.param(0) => 0.0, ep.param(1) => 0.0, ep.param(2) => 0.0,
                    );
                } else {
                    let s = v.signum();
                    let ql = |nk: f64, dk: f64| -s * (nk / d - v * (2.0 * dll * dk / ll));
                    let qn = |ldk: f64, nk: f64| -s * (ldk / d - v * (2.0 * dln * nk / ln));
                    row!(out[0],
                        el.param(0) => -ql(nx, ldx),
                        el.param(1) => -ql(ny, ldy),
                        el.param(2) => -ql(nz, ldz),
                        el.param(3) => ql(nx, ldx),
                        el.param(4) => ql(ny, ldy),
                        el.param(5) => ql(nz, ldz),
                        ep.param(0) => qn(ldx, nx),
                        ep.param(1) => qn(ldy, ny),
                        ep.param(2) => qn(ldz, nz),
                    );
                }
            }

            // ── PointPlaneDistance ───────────────────────────────────────────
            SpaceConstraint::PointPlaneDistance { point, plane, .. } => {
                let (px, py, pz) = point.param_indices(entities);
                let ep = &entities[*plane];
                let (nx, ny, nz) = plane_normal(params, ep);
                let d = params[ep.param(3)];
                let (ln, dln) = guarded_norm(nx * nx + ny * ny + nz * nz);
                let n_val = nx * params[px] + ny * params[py] + nz * params[pz] - d;
                let s = n_val / ln;
                let sign = if s > 0.0 {
                    1.0
                } else if s < 0.0 {
                    -1.0
                } else {
                    0.0
                };
                row!(out[0],
                    px => sign * nx / ln,
                    py => sign * ny / ln,
                    pz => sign * nz / ln,
                    ep.param(0) => sign * (params[px] / ln - s * 2.0 * nx * dln / ln),
                    ep.param(1) => sign * (params[py] / ln - s * 2.0 * ny * dln / ln),
                    ep.param(2) => sign * (params[pz] / ln - s * 2.0 * nz * dln / ln),
                    ep.param(3) => -sign / ln,
                );
            }

            // ── Concentric3D ─────────────────────────────────────────────────
            SpaceConstraint::Concentric3D(e1, e2) => {
                let a = &entities[*e1];
                let b = &entities[*e2];
                row!(out[0], a.param(0) => 1.0, b.param(0) => -1.0);
                row!(out[1], a.param(1) => 1.0, b.param(1) => -1.0);
                row!(out[2], a.param(2) => 1.0, b.param(2) => -1.0);
            }

            // ── CylinderRadius ───────────────────────────────────────────────
            SpaceConstraint::CylinderRadius { cylinder, .. } => {
                let e = &entities[*cylinder];
                row!(out[0], e.param(6) => 1.0);
            }

            // ── PointOnCylinder ──────────────────────────────────────────────
            SpaceConstraint::PointOnCylinder { point, cylinder } => {
                let (px, py, pz) = point.param_indices(entities);
                let ec = &entities[*cylinder];
                let (cx, cy, cz) = (
                    params[ec.param(0)],
                    params[ec.param(1)],
                    params[ec.param(2)],
                );
                let (ax, ay, az) = (
                    params[ec.param(3)],
                    params[ec.param(4)],
                    params[ec.param(5)],
                );
                let (l, _) = guarded_norm(ax * ax + ay * ay + az * az);
                let (nax, nay, naz) = (ax / l, ay / l, az / l);
                let vx = params[px] - cx;
                let vy = params[py] - cy;
                let vz = params[pz] - cz;
                let proj = vx * nax + vy * nay + vz * naz;
                let perp_x = vx - proj * nax;
                let perp_y = vy - proj * nay;
                let perp_z = vz - proj * naz;
                let (_, dpd) = guarded_norm(perp_x * perp_x + perp_y * perp_y + perp_z * perp_z);
                let c = 2.0 * dpd;
                let ka = -proj * c / l;
                row!(out[0],
                    px => c * perp_x,
                    py => c * perp_y,
                    pz => c * perp_z,
                    ec.param(0) => -c * perp_x,
                    ec.param(1) => -c * perp_y,
                    ec.param(2) => -c * perp_z,
                    ec.param(3) => ka * perp_x,
                    ec.param(4) => ka * perp_y,
                    ec.param(5) => ka * perp_z,
                    ec.param(6) => -1.0,
                );
            }

            // ── CylinderTangent ──────────────────────────────────────────────
            SpaceConstraint::CylinderTangent { c1, c2, external } => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                let (cx1, cy1, cz1) = (
                    params[e1.param(0)],
                    params[e1.param(1)],
                    params[e1.param(2)],
                );
                let (ax1, ay1, az1) = (
                    params[e1.param(3)],
                    params[e1.param(4)],
                    params[e1.param(5)],
                );
                let (l, _) = guarded_norm(ax1 * ax1 + ay1 * ay1 + az1 * az1);
                let (nax, nay, naz) = (ax1 / l, ay1 / l, az1 / l);
                let dx = params[e2.param(0)] - cx1;
                let dy = params[e2.param(1)] - cy1;
                let dz = params[e2.param(2)] - cz1;
                let proj = dx * nax + dy * nay + dz * naz;
                let perp_x = dx - proj * nax;
                let perp_y = dy - proj * nay;
                let perp_z = dz - proj * naz;
                let (_, dpd) = guarded_norm(perp_x * perp_x + perp_y * perp_y + perp_z * perp_z);
                let c = 2.0 * dpd;
                let ka = -proj * c / l;
                // The axis of the second cylinder does not enter the residual.
                let (r1, r2) = (
                    if *external {
                        -1.0
                    } else {
                        -(params[e1.param(6)] - params[e2.param(6)]).signum()
                    },
                    if *external {
                        -1.0
                    } else {
                        (params[e1.param(6)] - params[e2.param(6)]).signum()
                    },
                );
                row!(out[0],
                    e1.param(0) => -c * perp_x,
                    e1.param(1) => -c * perp_y,
                    e1.param(2) => -c * perp_z,
                    e1.param(3) => ka * perp_x,
                    e1.param(4) => ka * perp_y,
                    e1.param(5) => ka * perp_z,
                    e1.param(6) => r1,
                    e2.param(0) => c * perp_x,
                    e2.param(1) => c * perp_y,
                    e2.param(2) => c * perp_z,
                    e2.param(3) => 0.0,
                    e2.param(4) => 0.0,
                    e2.param(5) => 0.0,
                    e2.param(6) => r2,
                );
            }

            // ── PlaneTangentToSphere ─────────────────────────────────────────
            SpaceConstraint::PlaneTangentToSphere { plane, sphere } => {
                let ep = &entities[*plane];
                let es = &entities[*sphere];
                let (nx, ny, nz) = plane_normal(params, ep);
                let d = params[ep.param(3)];
                let (cx, cy, cz) = (
                    params[es.param(0)],
                    params[es.param(1)],
                    params[es.param(2)],
                );
                let (ln, dln) = guarded_norm(nx * nx + ny * ny + nz * nz);
                let n_val = nx * cx + ny * cy + nz * cz - d;
                let s = n_val / ln;
                let sign = if s > 0.0 {
                    1.0
                } else if s < 0.0 {
                    -1.0
                } else {
                    0.0
                };
                row!(out[0],
                    es.param(0) => sign * nx / ln,
                    es.param(1) => sign * ny / ln,
                    es.param(2) => sign * nz / ln,
                    es.param(3) => -1.0,
                    ep.param(0) => sign * (cx / ln - s * 2.0 * nx * dln / ln),
                    ep.param(1) => sign * (cy / ln - s * 2.0 * ny * dln / ln),
                    ep.param(2) => sign * (cz / ln - s * 2.0 * nz * dln / ln),
                    ep.param(3) => -sign / ln,
                );
            }

            // ── PointOnCircle / PointOnArc ───────────────────────────────────
            SpaceConstraint::PointOnCircle { point, circle }
            | SpaceConstraint::PointOnArc { point, arc: circle } => {
                let (px, py, pz) = point.param_indices(entities);
                let e = &entities[*circle];
                let frame = e.curve_frame();
                let (cx, cy, cz) = (params[e.param(0)], params[e.param(1)], params[e.param(2)]);
                let r = params[e.param(3)];
                let ux = params[px] - cx;
                let uy = params[py] - cy;
                let uz = params[pz] - cz;
                let n = frame.normal;
                let h = ux * n.x + uy * n.y + uz * n.z;
                // Row 0: u.n. Row 1: |u|^2 - (u.n)^2 - r^2; its gradient is
                // 2(u - (u.n)n), the in-plane radial direction scaled by 2r.
                let gx = 2.0 * (ux - h * n.x);
                let gy = 2.0 * (uy - h * n.y);
                let gz = 2.0 * (uz - h * n.z);
                row!(out[0],
                    px => n.x, py => n.y, pz => n.z,
                    e.param(0) => -n.x, e.param(1) => -n.y, e.param(2) => -n.z,
                    e.param(3) => 0.0,
                );
                row!(out[1],
                    px => gx, py => gy, pz => gz,
                    e.param(0) => -gx, e.param(1) => -gy, e.param(2) => -gz,
                    e.param(3) => -2.0 * r,
                );
            }

            // ── LineTangentCircle ────────────────────────────────────────────
            SpaceConstraint::LineTangentCircle { line, circle } => {
                let el = &entities[*line];
                let ec = &entities[*circle];
                let d = line_dir(params, el);
                let q = (
                    params[ec.param(0)] - params[el.param(0)],
                    params[ec.param(1)] - params[el.param(1)],
                    params[ec.param(2)] - params[el.param(2)],
                );
                let w = cross3(q, d);
                let (len, _) = guarded_norm(d.0 * d.0 + d.1 * d.1 + d.2 * d.2);
                let (s, ds) = guarded_norm(w.0 * w.0 + w.1 * w.1 + w.2 * w.2);
                let r = s / len;
                // grad_q s = (d x w)/s, grad_d s = (w x q)/s; scaled by 2*ds.
                let gq = cross3(d, w);
                let gd = cross3(w, q);
                let (gq, gd) = if ds > 0.0 {
                    (
                        (gq.0 * 2.0 * ds, gq.1 * 2.0 * ds, gq.2 * 2.0 * ds),
                        (gd.0 * 2.0 * ds, gd.1 * 2.0 * ds, gd.2 * 2.0 * ds),
                    )
                } else {
                    ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
                };
                row!(out[0],
                    ec.param(0) => gq.0 / len,
                    ec.param(1) => gq.1 / len,
                    ec.param(2) => gq.2 / len,
                    ec.param(3) => -1.0,
                    el.param(0) => -(gq.0 + gd.0) / len + r * d.0 / (len * len),
                    el.param(1) => -(gq.1 + gd.1) / len + r * d.1 / (len * len),
                    el.param(2) => -(gq.2 + gd.2) / len + r * d.2 / (len * len),
                    el.param(3) => gd.0 / len - r * d.0 / (len * len),
                    el.param(4) => gd.1 / len - r * d.1 / (len * len),
                    el.param(5) => gd.2 / len - r * d.2 / (len * len),
                );
            }

            // ── CircleCircleTangent ──────────────────────────────────────────
            SpaceConstraint::CircleCircleTangent { c1, c2, external } => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                let dx = params[e2.param(0)] - params[e1.param(0)];
                let dy = params[e2.param(1)] - params[e1.param(1)];
                let dz = params[e2.param(2)] - params[e1.param(2)];
                let (_, dd) = guarded_norm(dx * dx + dy * dy + dz * dz);
                let (ix, iy, iz) = (2.0 * dx * dd, 2.0 * dy * dd, 2.0 * dz * dd);
                let (r1, r2) = (
                    if *external {
                        -1.0
                    } else {
                        -(params[e1.param(3)] - params[e2.param(3)]).signum()
                    },
                    if *external {
                        -1.0
                    } else {
                        (params[e1.param(3)] - params[e2.param(3)]).signum()
                    },
                );
                row!(out[0],
                    e1.param(0) => -ix,
                    e1.param(1) => -iy,
                    e1.param(2) => -iz,
                    e1.param(3) => r1,
                    e2.param(0) => ix,
                    e2.param(1) => iy,
                    e2.param(2) => iz,
                    e2.param(3) => r2,
                );
            }

            // ── EqualRadius ──────────────────────────────────────────────────
            SpaceConstraint::EqualRadius(c1, c2) => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                row!(out[0], e1.param(3) => 1.0, e2.param(3) => -1.0);
            }

            // ── CircleRadius ─────────────────────────────────────────────────
            SpaceConstraint::CircleRadius { circle, .. } => {
                let e = &entities[*circle];
                row!(out[0], e.param(3) => 1.0);
            }

            // ── Midpoint ─────────────────────────────────────────────────────
            SpaceConstraint::Midpoint { point, line } => {
                let (px, py, pz) = point.param_indices(entities);
                let el = &entities[*line];
                row!(out[0],
                    px => 1.0,
                    el.param(0) => -0.5, el.param(3) => -0.5,
                );
                row!(out[1],
                    py => 1.0,
                    el.param(1) => -0.5, el.param(4) => -0.5,
                );
                row!(out[2],
                    pz => 1.0,
                    el.param(2) => -0.5, el.param(5) => -0.5,
                );
            }

            // ── Symmetric ────────────────────────────────────────────────────
            SpaceConstraint::Symmetric { p1, p2, line } => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                let el = &entities[*line];
                let (lx1, ly1, lz1) = (el.param(0), el.param(1), el.param(2));
                let (lx2, ly2, lz2) = (el.param(3), el.param(4), el.param(5));
                let d = (
                    params[lx2] - params[lx1],
                    params[ly2] - params[ly1],
                    params[lz2] - params[lz1],
                );
                let (len, dl) = guarded_norm(d.0 * d.0 + d.1 * d.1 + d.2 * d.2);
                let x = (
                    0.5 * (params[x1] + params[x2]) - params[lx1],
                    0.5 * (params[y1] + params[y2]) - params[ly1],
                    0.5 * (params[z1] + params[z2]) - params[lz1],
                );
                let c = cross3(x, d);
                let (i0, i1) = cross_pair(dominant_axis(d.0, d.1, d.2));
                for (row_idx, i) in [(0usize, i0), (1usize, i1)] {
                    // Same chain rule as PointOnLine, with mid = (p1+p2)/2.
                    let dc_point = cross_u_deriv(i, d);
                    let dc_end = cross_v_deriv(i, x);
                    let r = comp(c, i) / len;
                    let q = |dc: f64, dk: f64| dc / len - r * (2.0 * dl * dk / len);
                    row!(out[row_idx],
                        x1 => 0.5 * q(dc_point.0, 0.0),
                        y1 => 0.5 * q(dc_point.1, 0.0),
                        z1 => 0.5 * q(dc_point.2, 0.0),
                        x2 => 0.5 * q(dc_point.0, 0.0),
                        y2 => 0.5 * q(dc_point.1, 0.0),
                        z2 => 0.5 * q(dc_point.2, 0.0),
                        lx1 => q(-dc_point.0 - dc_end.0, -d.0),
                        ly1 => q(-dc_point.1 - dc_end.1, -d.1),
                        lz1 => q(-dc_point.2 - dc_end.2, -d.2),
                        lx2 => q(dc_end.0, d.0),
                        ly2 => q(dc_end.1, d.1),
                        lz2 => q(dc_end.2, d.2),
                    );
                }
                // Perpendicular joining segment: (p2-p1).d/|d|.
                let s = (
                    params[x2] - params[x1],
                    params[y2] - params[y1],
                    params[z2] - params[z1],
                );
                let e = (s.0 * d.0 + s.1 * d.1 + s.2 * d.2) / len;
                // d/d(line_k) of s/|d|: the sign splits between the two
                // endpoints; |d| also depends on both endpoints.
                let line_deriv = |k: usize| {
                    let (s_k, d_k) = (comp(s, k), comp(d, k));
                    s_k / len - e * d_k / (len * len)
                };
                row!(out[2],
                    x1 => -d.0 / len,
                    y1 => -d.1 / len,
                    z1 => -d.2 / len,
                    x2 => d.0 / len,
                    y2 => d.1 / len,
                    z2 => d.2 / len,
                    lx1 => -line_deriv(0),
                    ly1 => -line_deriv(1),
                    lz1 => -line_deriv(2),
                    lx2 => line_deriv(0),
                    ly2 => line_deriv(1),
                    lz2 => line_deriv(2),
                );
            }

            // ── PointAtArcEndpoint ───────────────────────────────────────────
            SpaceConstraint::PointAtArcEndpoint {
                point,
                arc,
                at_start,
            } => {
                let (px, py, pz) = point.param_indices(entities);
                let e = &entities[*arc];
                let frame = e.curve_frame();
                let angle = if *at_start {
                    frame.start_angle
                } else {
                    frame.normalized_end_angle()
                };
                let radial = frame.radial(angle);
                row!(out[0],
                    px => 1.0,
                    e.param(0) => -1.0,
                    e.param(3) => -radial.x,
                );
                row!(out[1],
                    py => 1.0,
                    e.param(1) => -1.0,
                    e.param(3) => -radial.y,
                );
                row!(out[2],
                    pz => 1.0,
                    e.param(2) => -1.0,
                    e.param(3) => -radial.z,
                );
            }

            // ── Spline constraints (3D-4) ──────────────────────────────────
            // The envelope theorem at the nearest point: the derivative of
            // the perpendicular residual through t* vanishes (the tangent is
            // orthogonal to the perpendicular plane), so d r/d p is the
            // perpendicular projector M = I - uu^T and d r/d cp_i is
            // -B_i(t*) * M — exact wherever the nearest point is regular,
            // matching the finite differences the residual is checked
            // against.
            SpaceConstraint::PointOnSpline { point, spline } => {
                let (px, py, pz) = point.param_indices(entities);
                let e = &entities[*spline];
                let n = e.control_point_count;
                if n < 2 {
                    // A degenerate spline carries no row (the residual arm
                    // reports zeros for it).
                }
                let cp: Vec<(f64, f64, f64)> = (0..n)
                    .map(|i| {
                        (
                            params[e.param(i * 3)],
                            params[e.param(i * 3 + 1)],
                            params[e.param(i * 3 + 2)],
                        )
                    })
                    .collect();
                let p = (params[px], params[py], params[pz]);
                let (t, _dist) = nearest_on_bezier(&cp, p);
                let (tx, ty, tz) = bezier_deriv(&cp, t);
                let t_len = (tx * tx + ty * ty + tz * tz).sqrt().max(SOLVER_NORM_TOL);
                let (ux, uy, uz) = (tx / t_len, ty / t_len, tz / t_len);
                // The projector's rows: d comp_k(M v) / d v = row k of M.
                let ex = (1.0 - ux * ux, -ux * uy, -ux * uz);
                let ey = (-uy * ux, 1.0 - uy * uy, -uy * uz);
                let ez = (-uz * ux, -uz * uy, 1.0 - uz * uz);
                let (i0, i1) = cross_pair(dominant_axis(tx, ty, tz));
                let c0 = [comp(ex, i0), comp(ey, i0), comp(ez, i0)];
                let c1 = [comp(ex, i1), comp(ey, i1), comp(ez, i1)];
                out[0].push((px, c0[0]));
                out[0].push((py, c0[1]));
                out[0].push((pz, c0[2]));
                out[1].push((px, c1[0]));
                out[1].push((py, c1[1]));
                out[1].push((pz, c1[2]));
                let k = n - 1;
                for i in 0..n {
                    // The Bernstein basis value B_i(t*), negated for -c'.
                    let mut binom = 1.0;
                    for j in 0..i {
                        binom *= (k - j) as f64 / (j + 1) as f64;
                    }
                    let b = -binom * t.powi(i as i32) * (1.0 - t).powi((k - i) as i32);
                    for a in 0..3 {
                        let param = e.param(i * 3 + a);
                        out[0].push((param, b * c0[a]));
                        out[1].push((param, b * c1[a]));
                    }
                }
            }

            SpaceConstraint::SplineEndpoint {
                point,
                spline,
                at_start,
            } => {
                let (px, py, pz) = point.param_indices(entities);
                let e = &entities[*spline];
                let idx = if *at_start {
                    0
                } else {
                    e.control_point_count.saturating_sub(1)
                };
                row!(out[0],
                    px => 1.0,
                    e.param(idx * 3) => -1.0,
                );
                row!(out[1],
                    py => 1.0,
                    e.param(idx * 3 + 1) => -1.0,
                );
                row!(out[2],
                    pz => 1.0,
                    e.param(idx * 3 + 2) => -1.0,
                );
            }

            SpaceConstraint::SplineControlPointMidpoint {
                spline,
                midpoint_index,
                first_index,
                second_index,
            } => {
                let e = &entities[*spline];
                let mi = *midpoint_index * 3;
                let fi = *first_index * 3;
                let si = *second_index * 3;
                row!(out[0],
                    e.param(mi) => 1.0,
                    e.param(fi) => -0.5,
                    e.param(si) => -0.5,
                );
                row!(out[1],
                    e.param(mi + 1) => 1.0,
                    e.param(fi + 1) => -0.5,
                    e.param(si + 1) => -0.5,
                );
                row!(out[2],
                    e.param(mi + 2) => 1.0,
                    e.param(fi + 2) => -0.5,
                    e.param(si + 2) => -0.5,
                );
            }

            SpaceConstraint::SplineControlPointSymmetric {
                spline,
                first_index,
                second_index,
                line,
            } => {
                let e = &entities[*spline];
                let fi = *first_index * 3;
                let si = *second_index * 3;
                let el = &entities[*line];
                let lx1 = params[el.param(0)];
                let ly1 = params[el.param(1)];
                let lz1 = params[el.param(2)];
                let d = line_dir(params, el);
                let (len, _) = guarded_norm(d.0 * d.0 + d.1 * d.1 + d.2 * d.2);
                // Midpoint lies on the line: two independent cross rows.
                let x = (
                    0.5 * (params[e.param(fi)] + params[e.param(si)]) - lx1,
                    0.5 * (params[e.param(fi + 1)] + params[e.param(si + 1)]) - ly1,
                    0.5 * (params[e.param(fi + 2)] + params[e.param(si + 2)]) - lz1,
                );
                let c = cross3(x, d);
                let (i0, i1) = cross_pair(dominant_axis(d.0, d.1, d.2));
                // d(cross)/d(midpoint) = cross3(0.5, d) for each control point.
                let half_d = (0.5 * d.0, 0.5 * d.1, 0.5 * d.2);
                row!(out[0],
                    e.param(fi) => comp(cross3((1.0, 0.0, 0.0), d), i0) / len,
                    e.param(si) => comp(cross3((1.0, 0.0, 0.0), d), i0) / len,
                );
                row!(out[1],
                    e.param(fi + 1) => comp(cross3((0.0, 1.0, 0.0), d), i1) / len,
                    e.param(si + 1) => comp(cross3((0.0, 1.0, 0.0), d), i1) / len,
                );
                // The joining segment is perpendicular to the line.
                row!(out[2],
                    e.param(fi) => -d.0 / len,
                    e.param(si) => d.0 / len,
                    e.param(fi + 1) => -d.1 / len,
                    e.param(si + 1) => d.1 / len,
                    e.param(fi + 2) => -d.2 / len,
                    e.param(si + 2) => d.2 / len,
                );
                let _ = half_d;
                let _ = c;
            }
        }
    }
}

/// Reflexive `AsRef` so the solver's residual/Jacobian helpers accept both
/// owned and borrowed constraint slices.
impl AsRef<SpaceConstraint> for SpaceConstraint {
    fn as_ref(&self) -> &SpaceConstraint {
        self
    }
}

/// Index of the coordinate with the largest magnitude, ties broken toward
/// the lowest index (deterministic).
#[inline]
fn dominant_axis(x: f64, y: f64, z: f64) -> usize {
    let (ax, ay, az) = (x.abs(), y.abs(), z.abs());
    if ax >= ay && ax >= az {
        0
    } else if ay >= az {
        1
    } else {
        2
    }
}

/// The two component indices retained for a cross-product congruence test,
/// omitting the dominant axis `k`: the cross is orthogonal to the reference
/// direction, so the omitted component is a linear combination of the other
/// two (and the two retained rows are independent whenever the reference's
/// `k` coordinate is nonzero).
#[inline]
fn cross_pair(k: usize) -> (usize, usize) {
    match k {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    }
}

#[inline]
fn cross3(u: (f64, f64, f64), v: (f64, f64, f64)) -> (f64, f64, f64) {
    (
        u.1 * v.2 - u.2 * v.1,
        u.2 * v.0 - u.0 * v.2,
        u.0 * v.1 - u.1 * v.0,
    )
}

#[inline]
fn comp(c: (f64, f64, f64), i: usize) -> f64 {
    match i {
        0 => c.0,
        1 => c.1,
        _ => c.2,
    }
}

/// Derivative of `(u x v)_i` with respect to the components of `u`.
#[inline]
fn cross_u_deriv(i: usize, v: (f64, f64, f64)) -> (f64, f64, f64) {
    match i {
        0 => (0.0, v.2, -v.1),
        1 => (-v.2, 0.0, v.0),
        _ => (v.1, -v.0, 0.0),
    }
}

/// Derivative of `(u x v)_i` with respect to the components of `v`.
#[inline]
fn cross_v_deriv(i: usize, u: (f64, f64, f64)) -> (f64, f64, f64) {
    match i {
        0 => (0.0, -u.2, u.1),
        1 => (u.2, 0.0, -u.0),
        _ => (-u.1, u.0, 0.0),
    }
}

/// Guarded Euclidean norm matching the residual formulas' `.max(SOLVER_NORM_TOL)`:
/// returns `(guarded norm, d norm / d squared-norm)`. The derivative is zero
/// while the floor is active (the residual then treats the norm as a constant).
#[inline]
fn guarded_norm(sq: f64) -> (f64, f64) {
    let l = sq.sqrt();
    if l > SOLVER_NORM_TOL {
        (l, 0.5 / l)
    } else {
        (SOLVER_NORM_TOL, 0.0)
    }
}

/// Direction vector of a line entity, `end - start`.
#[inline]
fn line_dir(params: &[f64], e: &SpaceEntity) -> (f64, f64, f64) {
    (
        params[e.param(3)] - params[e.param(0)],
        params[e.param(4)] - params[e.param(1)],
        params[e.param(5)] - params[e.param(2)],
    )
}

/// Normal vector of a plane entity.
#[inline]
fn plane_normal(params: &[f64], e: &SpaceEntity) -> (f64, f64, f64) {
    (params[e.param(0)], params[e.param(1)], params[e.param(2)])
}

// ── Bezier curve evaluation (3D-4) ─────────────────────────────────────────

/// Binomial coefficient C(n, k).
fn binom(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut result = 1.0;
    for i in 0..k {
        result = result * (n - i) as f64 / (i + 1) as f64;
    }
    result
}

/// Evaluate a Bezier curve at parameter `t`.
///
/// `cp` are the control points; the curve degree is `cp.len() - 1`.
fn bezier_eval(cp: &[(f64, f64, f64)], t: f64) -> (f64, f64, f64) {
    let n = cp.len() - 1;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for (i, &(px, py, pz)) in cp.iter().enumerate() {
        let basis = binom(n, i) * t.powi(i as i32) * (1.0 - t).powi((n - i) as i32);
        x += basis * px;
        y += basis * py;
        z += basis * pz;
    }
    (x, y, z)
}

/// Evaluate the first derivative of a Bezier curve at parameter `t`.
fn bezier_deriv(cp: &[(f64, f64, f64)], t: f64) -> (f64, f64, f64) {
    let n = cp.len() - 1;
    if n == 0 {
        return (0.0, 0.0, 0.0);
    }
    // The derivative of a degree-n Bezier is a degree-(n-1) Bezier with
    // control points n * (P_{i+1} - P_i).
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for i in 0..n {
        let basis = binom(n - 1, i) * t.powi(i as i32) * (1.0 - t).powi((n - 1 - i) as i32);
        let dx = n as f64 * (cp[i + 1].0 - cp[i].0);
        let dy = n as f64 * (cp[i + 1].1 - cp[i].1);
        let dz = n as f64 * (cp[i + 1].2 - cp[i].2);
        x += basis * dx;
        y += basis * dy;
        z += basis * dz;
    }
    (x, y, z)
}

/// Find the parameter `t` of the point on the Bezier curve nearest to `p`.
///
/// Uses a coarse sampling followed by Newton refinement on the distance
/// function. Returns `(t, distance)`.
fn nearest_on_bezier(cp: &[(f64, f64, f64)], p: (f64, f64, f64)) -> (f64, f64) {
    const SAMPLES: usize = 20;
    let mut best_t = 0.0;
    let mut best_dist = f64::INFINITY;
    for i in 0..=SAMPLES {
        let t = i as f64 / SAMPLES as f64;
        let (cx, cy, cz) = bezier_eval(cp, t);
        let dx = p.0 - cx;
        let dy = p.1 - cy;
        let dz = p.2 - cz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist < best_dist {
            best_dist = dist;
            best_t = t;
        }
    }
    // Newton refinement on f(t) = (p - C(t)) . C'(t) = 0.
    let mut t = best_t;
    for _ in 0..10 {
        let (cx, cy, cz) = bezier_eval(cp, t);
        let (tx, ty, tz) = bezier_deriv(cp, t);
        let vx = p.0 - cx;
        let vy = p.1 - cy;
        let vz = p.2 - cz;
        let f = vx * tx + vy * ty + vz * tz;
        // f'(t) = -|C'(t)|^2 + (p - C(t)) . C''(t)
        // Approximate C''(t) numerically.
        let h = 1.0e-6;
        let t0 = (t - h).max(0.0);
        let t1 = (t + h).min(1.0);
        let (_, _, _) = (t0, t1, h);
        let (tx0, ty0, tz0) = bezier_deriv(cp, t0);
        let (tx1, ty1, tz1) = bezier_deriv(cp, t1);
        let dtx = (tx1 - tx0) / (t1 - t0);
        let dty = (ty1 - ty0) / (t1 - t0);
        let dtz = (tz1 - tz0) / (t1 - t0);
        let df = -(tx * tx + ty * ty + tz * tz) + vx * dtx + vy * dty + vz * dtz;
        if df.abs() < 1.0e-12 {
            break;
        }
        let dt = f / df;
        t = (t - dt).clamp(0.0, 1.0);
    }
    let (cx, cy, cz) = bezier_eval(cp, t);
    let dx = p.0 - cx;
    let dy = p.1 - cy;
    let dz = p.2 - cz;
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
    (t, dist)
}

#[cfg(test)]
mod jacobian_tests {
    use super::*;
    use crate::space3d::entity::{SpaceEntityKind, SpaceFrame3};

    /// Fixture entities: index 0/1 points, 2/3 lines, 4/5 planes,
    /// 6/7 spheres, 8/9 cylinders, 10/11 circles, 12/13 arcs.
    fn fixture() -> (Vec<f64>, Vec<SpaceEntity>) {
        let params = vec![
            0.3, -0.7, 1.1, // point 0
            2.2, 0.4, -0.9, // point 1
            1.3, -0.2, 0.5, 2.7, 1.4, -1.3, // line 2
            0.2, 1.1, -0.4, -1.6, 2.3, 0.8, // line 3
            0.6, -0.48, 0.64, 1.7, // plane 4
            0.0, 0.8, 0.6, -0.4, // plane 5
            0.5, 1.5, -2.0, 1.2, // sphere 6
            -1.1, 0.3, 0.9, 2.3, // sphere 7
            0.4, -0.6, 1.2, 0.3, 0.9, -0.2, 0.8, // cylinder 8
            1.7, 0.2, -0.5, -0.5, 0.4, 0.7, 1.4, // cylinder 9
            1.2, -0.8, 0.4, 1.5, // circle 10
            -0.6, 1.4, 0.9, 2.1, // circle 11
            0.7, 0.2, -1.1, 1.3, // arc 12
            -1.3, -0.5, 0.8, 0.9, // arc 13
        ];
        let frame = |normal: glam::DVec3, x: glam::DVec3, a0: f64, a1: f64| {
            SpaceFrame3::sanitized(normal, x, a0, a1)
        };
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 3),
            SpaceEntity::new(SpaceEntityKind::SpaceLine, 6),
            SpaceEntity::new(SpaceEntityKind::SpaceLine, 12),
            SpaceEntity::new(SpaceEntityKind::Plane, 18),
            SpaceEntity::new(SpaceEntityKind::Plane, 22),
            SpaceEntity::new(SpaceEntityKind::Sphere, 26),
            SpaceEntity::new(SpaceEntityKind::Sphere, 30),
            SpaceEntity::new(SpaceEntityKind::Cylinder, 34),
            SpaceEntity::new(SpaceEntityKind::Cylinder, 41),
            SpaceEntity::with_frame(
                SpaceEntityKind::SpaceCircle,
                48,
                frame(glam::DVec3::Z, glam::DVec3::X, 0.0, 0.0),
            ),
            SpaceEntity::with_frame(
                SpaceEntityKind::SpaceCircle,
                52,
                frame(
                    glam::DVec3::new(0.6, -0.48, 0.64),
                    glam::DVec3::new(0.0, 0.8, 0.6),
                    0.0,
                    0.0,
                ),
            ),
            SpaceEntity::with_frame(
                SpaceEntityKind::SpaceArc,
                56,
                frame(
                    glam::DVec3::new(0.6, -0.48, 0.64),
                    glam::DVec3::new(0.8, 0.6, 0.0),
                    -1.2,
                    1.1,
                ),
            ),
            SpaceEntity::with_frame(
                SpaceEntityKind::SpaceArc,
                60,
                frame(
                    glam::DVec3::new(0.0, 0.6, 0.8),
                    glam::DVec3::new(0.0, -0.8, 0.6),
                    2.0,
                    5.5,
                ),
            ),
        ];
        (params, entities)
    }

    /// Every variant of the 3D constraint vocabulary at two parameter
    /// points. A variant added to `residuals` without a matching
    /// `jacobian` arm cannot compile; a wrong derivative fails here.
    fn all_variants() -> Vec<SpaceConstraint> {
        use SpacePointRef::{ArcCenter, CircleCenter, LineEnd, LineStart, Point, SphereCenter};
        vec![
            SpaceConstraint::Coincident(Point(0), Point(1)),
            SpaceConstraint::Fixed {
                point: Point(0),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::Fixed {
                point: LineEnd(2),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::Fixed {
                point: SphereCenter(7),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::PointDistance {
                p1: Point(0),
                p2: LineStart(2),
                distance: 2.5,
            },
            SpaceConstraint::PointDistanceAlongDir {
                p1: Point(0),
                p2: LineEnd(3),
                dx: 0.2,
                dy: -0.7,
                dz: 0.4,
                distance: 1.5,
            },
            SpaceConstraint::LineLength {
                line: 2,
                length: 3.0,
            },
            SpaceConstraint::PointOnLine {
                point: Point(0),
                line: 2,
            },
            SpaceConstraint::PointOnPlane {
                point: Point(0),
                plane: 4,
            },
            SpaceConstraint::PointOnSphere {
                point: Point(0),
                sphere: 6,
            },
            SpaceConstraint::LineParallelToDir {
                line: 2,
                dx: 0.2,
                dy: 0.9,
                dz: -0.3,
            },
            SpaceConstraint::LineParallelLine(2, 3),
            SpaceConstraint::LinePerpendicularLine(2, 3),
            SpaceConstraint::EqualLineLength(2, 3),
            SpaceConstraint::PlaneNormal {
                plane: 4,
                nx: 0.1,
                ny: 0.7,
                nz: 0.3,
            },
            SpaceConstraint::PlaneParallel(4, 5),
            SpaceConstraint::PlanePerpendicular(4, 5),
            SpaceConstraint::LineInPlane { line: 2, plane: 4 },
            SpaceConstraint::LinePerpendicularPlane { line: 2, plane: 4 },
            SpaceConstraint::LineParallelPlane { line: 2, plane: 4 },
            SpaceConstraint::SphereRadius {
                sphere: 6,
                radius: 1.0,
            },
            SpaceConstraint::SphereTangent {
                s1: 6,
                s2: 7,
                external: true,
            },
            SpaceConstraint::SphereTangent {
                s1: 6,
                s2: 7,
                external: false,
            },
            SpaceConstraint::EqualSphereRadius(6, 7),
            SpaceConstraint::AngleLineLine {
                l1: 2,
                l2: 3,
                angle_rad: 0.7,
            },
            SpaceConstraint::AnglePlanePlane {
                p1: 4,
                p2: 5,
                angle_rad: 0.6,
            },
            SpaceConstraint::AngleLinePlane {
                line: 2,
                plane: 4,
                angle_rad: 0.4,
            },
            SpaceConstraint::PointPlaneDistance {
                point: Point(0),
                plane: 4,
                distance: 0.5,
            },
            SpaceConstraint::Concentric3D(6, 7),
            SpaceConstraint::CylinderRadius {
                cylinder: 8,
                radius: 0.8,
            },
            SpaceConstraint::PointOnCylinder {
                point: Point(0),
                cylinder: 8,
            },
            SpaceConstraint::CylinderTangent {
                c1: 8,
                c2: 9,
                external: true,
            },
            SpaceConstraint::CylinderTangent {
                c1: 8,
                c2: 9,
                external: false,
            },
            SpaceConstraint::PlaneTangentToSphere {
                plane: 4,
                sphere: 6,
            },
            // ── 3D-3 circle/arc family ───────────────────────────────────────
            SpaceConstraint::PointOnCircle {
                point: Point(0),
                circle: 10,
            },
            SpaceConstraint::PointOnArc {
                point: Point(1),
                arc: 12,
            },
            SpaceConstraint::LineTangentCircle {
                line: 2,
                circle: 10,
            },
            SpaceConstraint::LineTangentCircle {
                line: 3,
                circle: 12,
            },
            SpaceConstraint::CircleCircleTangent {
                c1: 10,
                c2: 11,
                external: true,
            },
            SpaceConstraint::CircleCircleTangent {
                c1: 10,
                c2: 11,
                external: false,
            },
            SpaceConstraint::EqualRadius(10, 11),
            SpaceConstraint::CircleRadius {
                circle: 10,
                radius: 1.5,
            },
            SpaceConstraint::Midpoint {
                point: Point(0),
                line: 2,
            },
            SpaceConstraint::Symmetric {
                p1: Point(0),
                p2: Point(1),
                line: 3,
            },
            SpaceConstraint::PointAtArcEndpoint {
                point: Point(0),
                arc: 12,
                at_start: true,
            },
            SpaceConstraint::PointAtArcEndpoint {
                point: LineEnd(3),
                arc: 13,
                at_start: false,
            },
            SpaceConstraint::Concentric3D(10, 12),
            SpaceConstraint::Fixed {
                point: CircleCenter(10),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::Fixed {
                point: ArcCenter(13),
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            SpaceConstraint::PointDistance {
                p1: CircleCenter(11),
                p2: ArcCenter(12),
                distance: 2.0,
            },
        ]
    }

    fn check_one(c: &SpaceConstraint, params: &[f64], entities: &[SpaceEntity]) {
        let n_eq = c.equation_count();
        let mut f = vec![0.0_f64; n_eq];
        assert!(c.residuals(params, entities, &mut f));
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_eq];
        c.jacobian(params, entities, &mut rows);

        let h = 1.0e-6;
        for col in 0..params.len() {
            let mut p_plus = params.to_vec();
            p_plus[col] += h;
            let mut f_plus = vec![0.0_f64; n_eq];
            assert!(c.residuals(&p_plus, entities, &mut f_plus));
            let mut p_minus = params.to_vec();
            p_minus[col] -= h;
            let mut f_minus = vec![0.0_f64; n_eq];
            assert!(c.residuals(&p_minus, entities, &mut f_minus));

            for r in 0..n_eq {
                let fd = (f_plus[r] - f_minus[r]) / (2.0 * h);
                let mut analytic = 0.0_f64;
                for &(p, v) in &rows[r] {
                    if p == col {
                        analytic += v;
                    }
                }
                assert!(
                    (analytic - fd).abs() <= 1.0e-5 * (1.0 + fd.abs()),
                    "analytic Jacobian drifts from finite differences at row {r} col {col}: \
                     analytic {analytic} vs finite difference {fd} for {c:?}"
                );
            }
        }
    }

    #[test]
    fn analytic_jacobian_matches_finite_difference() {
        let (params, entities) = fixture();
        // Two evaluation points: the fixture and a deterministic perturbation,
        // so clamp/guard branches are exercised away from a single state.
        let mut perturbed = params.clone();
        for (i, v) in perturbed.iter_mut().enumerate() {
            *v += 0.037 * ((i % 7) as f64 - 3.0);
        }
        for c in all_variants() {
            check_one(&c, &params, &entities);
            check_one(&c, &perturbed, &entities);
        }
    }

    /// Structural completeness: every variant must report at least one
    /// touched parameter on every equation row (otherwise the constraint
    /// graph would drop its entities from the component it belongs to).
    #[test]
    fn every_row_touches_parameters() {
        let (params, entities) = fixture();
        for c in all_variants() {
            let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); c.equation_count()];
            c.jacobian(&params, &entities, &mut rows);
            for (r, row) in rows.iter().enumerate() {
                assert!(!row.is_empty(), "row {r} of {c:?} touches no parameter");
            }
        }
    }

    /// The PointOnLine rows select their cross-product pair from the line
    /// direction's dominant axis. The Jacobian must stay exact on every
    /// world-axis-aligned support line (the fixed y/z row pair used to
    /// degenerate on X/Y-aligned lines).
    #[test]
    fn point_on_line_jacobian_is_exact_for_axis_aligned_lines() {
        for direction in [(10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0)] {
            let params = vec![
                1.0,
                -2.0,
                3.0, // point
                4.0,
                5.0,
                6.0, // line start
                4.0 + direction.0,
                5.0 + direction.1,
                6.0 + direction.2, // line end
            ];
            let entities = vec![
                SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
                SpaceEntity::new(SpaceEntityKind::SpaceLine, 3),
            ];
            let c = SpaceConstraint::PointOnLine {
                point: SpacePointRef::Point(0),
                line: 1,
            };
            check_one(&c, &params, &entities);
        }
    }
}

#[cfg(test)]
mod point_on_line_tests {
    use super::*;
    use crate::space3d::sketch::SpaceSketch;

    /// Regression for the axis-aligned degeneracy: the point must reach the
    /// support line on X-, Y- and Z-aligned lines alike. The old fixed
    /// y/z cross-row pair left one coordinate free for X- and Y-aligned
    /// lines.
    #[test]
    fn point_reaches_every_axis_aligned_line() {
        for (start, end) in [
            ((1.0, -2.0, 3.0), (11.0, -2.0, 3.0)),
            ((1.0, -2.0, 3.0), (1.0, 8.0, 3.0)),
            ((1.0, -2.0, 3.0), (1.0, -2.0, 13.0)),
        ] {
            let mut sk = SpaceSketch::new();
            let line = sk.add_line(start.0, start.1, start.2, end.0, end.1, end.2);
            sk.fix_entity(line);
            // Offset in both transverse coordinates (and beyond the end
            // point, exercising the infinite support line).
            let point = sk.add_point(20.0, 6.0, 9.0);
            sk.add_constraint(SpaceConstraint::PointOnLine {
                point: SpacePointRef::Point(point),
                line,
            });

            let diagnosis = sk.diagnose();
            assert_eq!(
                diagnosis.dof, 1,
                "slide along the line remains: {diagnosis:?}"
            );
            assert!(diagnosis.conflicting.is_empty(), "{diagnosis:?}");

            let result = sk.solve();
            assert!(
                result.converged,
                "residual {} for {start:?}->{end:?}",
                result.residual
            );

            let (x, y, z) = sk.point_coords(point).unwrap();
            let p = glam::DVec3::new(x, y, z);
            let a = glam::DVec3::new(start.0, start.1, start.2);
            let b = glam::DVec3::new(end.0, end.1, end.2);
            let dir = (b - a).normalize();
            let transverse = (p - a) - (p - a).dot(dir) * dir;
            assert!(
                transverse.length() < 1.0e-9,
                "point {p:?} is {} off the line {a:?}->{b:?}",
                transverse.length()
            );
        }
    }

    /// A point on a diagonal line keeps both transverse rows independent.
    #[test]
    fn point_reaches_a_diagonal_line() {
        let mut sk = SpaceSketch::new();
        let line = sk.add_line(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        sk.fix_entity(line);
        let point = sk.add_point(4.0, -2.0, 1.0);
        sk.add_constraint(SpaceConstraint::PointOnLine {
            point: SpacePointRef::Point(point),
            line,
        });
        let result = sk.solve();
        assert!(result.converged, "residual {}", result.residual);
        let (x, y, z) = sk.point_coords(point).unwrap();
        assert!((x - y).abs() < 1.0e-9, "({x}, {y}, {z})");
        assert!((y - z).abs() < 1.0e-9, "({x}, {y}, {z})");
    }
}
