//! 3D geometric entity types and their parameter layout.
//!
//! All entity parameters are stored in a flat `Vec<f64>` owned by [`crate::space3d::SpaceSketch`].
//! Each entity records its `param_start` offset into that vector.
//!
//! | Kind        | Params (relative offsets)                              |
//! |-------------|--------------------------------------------------------|
//! | SpacePoint  | `[0]`=x, `[1]`=y, `[2]`=z                                    |
//! | SpaceLine   | `[0]`=x1, `[1]`=y1, `[2]`=z1, `[3]`=x2, `[4]`=y2, `[5]`=z2       |
//! | Plane       | `[0]`=nx, `[1]`=ny, `[2]`=nz, `[3]`=d (distance from origin)  |
//! | Sphere      | `[0]`=cx, `[1]`=cy, `[2]`=cz, `[3]`=r                          |
//! | Cylinder    | `[0..3]`=axis point, `[3..6]`=axis direction, `[6]`=radius  |
//! | SpaceCircle | `[0]`=cx, `[1]`=cy, `[2]`=cz, `[3]`=r, plus a fixed frame      |
//! | SpaceArc    | `[0]`=cx, `[1]`=cy, `[2]`=cz, `[3]`=r, plus a fixed frame      |
//!
//! Circles and arcs carry their plane orientation as *frame metadata*
//! ([`SpaceFrame3`]), not as parameters: a 3D-3 curve always lies on a
//! workplane (the sketch workplane or an explicit per-element frame), so the
//! orientation is datum geometry held constant. Their free parameters are the
//! center (3) and the radius (1) — exactly the degrees of freedom the solver
//! may move, matching the plane-local reduction of the editable model.

use glam::DVec3;

/// Index into the entity list of a [`crate::space3d::SpaceSketch`].
pub type SpaceEntityId = usize;

/// The geometric kind of a 3D entity.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpaceEntityKind {
    SpacePoint,
    SpaceLine,
    Plane,
    Sphere,
    Cylinder,
    /// A full circle: center + radius parameters over a fixed frame.
    SpaceCircle,
    /// A circular arc: center + radius parameters over a fixed frame plus a
    /// start/end angle span in the frame.
    SpaceArc,
    /// A Bezier spline: `3 * num_control_points` parameters (each control
    /// point is a world-space (x, y, z) triple).
    SpaceSpline,
}

impl SpaceEntityKind {
    /// Number of scalar parameters this entity type uses.
    ///
    /// For [`SpaceEntityKind::SpaceSpline`] this returns 0 because the count
    /// depends on the number of control points (stored per-entity in
    /// [`SpaceEntity::control_point_count`]); callers must use
    /// [`SpaceEntity::param_count()` for the per-entity count.
    pub fn param_count(self) -> usize {
        match self {
            SpaceEntityKind::SpacePoint => 3,
            SpaceEntityKind::SpaceLine => 6,
            SpaceEntityKind::Plane => 4,
            SpaceEntityKind::Sphere => 4,
            SpaceEntityKind::Cylinder => 7,
            SpaceEntityKind::SpaceCircle => 4,
            SpaceEntityKind::SpaceArc => 4,
            SpaceEntityKind::SpaceSpline => 0,
        }
    }

    /// Whether this kind is a circular curve (circle or arc).
    pub fn is_curve(self) -> bool {
        matches!(
            self,
            SpaceEntityKind::SpaceCircle | SpaceEntityKind::SpaceArc
        )
    }
}

/// The fixed in-plane frame of a circle/arc entity.
///
/// `normal` and `x_dir` are unit length and orthogonal; the in-plane Y axis
/// is `normal × x_dir`. Angles are measured in that frame (radians, right
/// handed around `normal`). The frame origin is the entity center itself.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpaceFrame3 {
    /// Unit plane normal.
    pub normal: DVec3,
    /// Unit in-plane X direction (start-angle reference).
    pub x_dir: DVec3,
    /// Arc start angle in radians (ignored by full circles).
    pub start_angle: f64,
    /// Arc end angle in radians (ignored by full circles).
    pub end_angle: f64,
}

impl SpaceFrame3 {
    /// The canonical world XY frame at angle zero (the fallback frame).
    pub const XY: SpaceFrame3 = SpaceFrame3 {
        normal: DVec3::Z,
        x_dir: DVec3::X,
        start_angle: 0.0,
        end_angle: 0.0,
    };

    /// Builds a sanitized frame deterministically.
    ///
    /// Degenerate input never panics and never produces NaN:
    /// - a non-finite or zero normal falls back to `+Z`;
    /// - the X direction is projected into the plane; a non-finite,
    ///   zero-length or parallel X direction falls back to a deterministic
    ///   axis helper (`gp_Ax2`-style: the coordinate axis least aligned with
    ///   the normal);
    /// - non-finite angles fall back to zero.
    pub fn sanitized(normal: DVec3, x_dir: DVec3, start_angle: f64, end_angle: f64) -> Self {
        let normal = if normal.is_finite() {
            normal.try_normalize().unwrap_or(DVec3::Z)
        } else {
            DVec3::Z
        };
        let projected = if x_dir.is_finite() {
            x_dir - normal * x_dir.dot(normal)
        } else {
            DVec3::ZERO
        };
        let x_dir = projected.try_normalize().unwrap_or_else(|| {
            let (ax, ay, az) = (normal.x.abs(), normal.y.abs(), normal.z.abs());
            let helper = if ax <= ay && ax <= az {
                DVec3::X
            } else if ay <= az {
                DVec3::Y
            } else {
                DVec3::Z
            };
            normal.cross(helper).normalize()
        });
        Self {
            normal,
            x_dir,
            start_angle: if start_angle.is_finite() {
                start_angle
            } else {
                0.0
            },
            end_angle: if end_angle.is_finite() {
                end_angle
            } else {
                0.0
            },
        }
    }

    /// The derived in-plane Y direction.
    #[inline]
    pub fn y_dir(&self) -> DVec3 {
        self.normal.cross(self.x_dir)
    }

    /// Unit radial direction at frame angle `angle`.
    #[inline]
    pub fn radial(&self, angle: f64) -> DVec3 {
        self.x_dir * angle.cos() + self.y_dir() * angle.sin()
    }

    /// The arc span in radians, normalized into `(0, 2π]`.
    ///
    /// The end angle is re-expressed as `start + span`, so the arc always
    /// travels counter-clockwise (around `normal`) from the start angle.
    pub fn arc_span(&self) -> f64 {
        let tau = std::f64::consts::TAU;
        let mut span = (self.end_angle - self.start_angle) % tau;
        if span <= 0.0 {
            span += tau;
        }
        span
    }

    /// The normalized end angle (`start_angle + arc_span()`).
    pub fn normalized_end_angle(&self) -> f64 {
        self.start_angle + self.arc_span()
    }

    /// Re-expresses `angle` as `start_angle + t` with `t` in `[0, 2π)`,
    /// i.e. the counter-clockwise offset from the arc start.
    pub fn span_offset(&self, angle: f64) -> f64 {
        let tau = std::f64::consts::TAU;
        let mut t = (angle - self.start_angle) % tau;
        if t < 0.0 {
            t += tau;
        }
        t
    }

    /// Whether `angle` falls inside the arc span.
    pub fn span_contains(&self, angle: f64) -> bool {
        self.span_offset(angle) <= self.arc_span()
    }

    /// The angle of the arc point nearest to `angle`: the angle itself when
    /// it is inside the span, otherwise the nearer span endpoint. Ties go to
    /// the start endpoint (deterministic).
    pub fn nearest_span_angle(&self, angle: f64) -> f64 {
        let span = self.arc_span();
        let t = self.span_offset(angle);
        if t <= span {
            return self.start_angle + t;
        }
        let tau = std::f64::consts::TAU;
        let to_end = t - span;
        let to_start = tau - t;
        if to_end < to_start {
            self.start_angle + span
        } else {
            self.start_angle
        }
    }

    /// The arc point nearest to a world point, given the arc's center and
    /// radius: in-plane when the projection falls inside the span, otherwise
    /// the nearer endpoint.
    pub fn nearest_arc_point(&self, center: DVec3, radius: f64, point: DVec3) -> DVec3 {
        let d = point - center;
        let u = d.dot(self.x_dir);
        let v = d.dot(self.y_dir());
        let angle = self.nearest_span_angle(v.atan2(u));
        center + self.radial(angle) * radius
    }
}

/// Metadata for a single 3D entity stored in a [`crate::space3d::SpaceSketch`].
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub struct SpaceEntity {
    pub kind: SpaceEntityKind,
    /// Index of the first parameter in the sketch's flat `params` vector.
    pub param_start: usize,
    /// Fixed in-plane frame for [`SpaceEntityKind::SpaceCircle`] and
    /// [`SpaceEntityKind::SpaceArc`]; `None` for every other kind.
    #[cfg_attr(feature = "serde", serde(default))]
    pub frame: Option<SpaceFrame3>,
    /// Number of control points for [`SpaceEntityKind::SpaceSpline`]; 0 for
    /// every other kind.
    #[cfg_attr(feature = "serde", serde(default))]
    pub control_point_count: usize,
}

impl SpaceEntity {
    pub fn new(kind: SpaceEntityKind, param_start: usize) -> Self {
        Self {
            kind,
            param_start,
            frame: None,
            control_point_count: 0,
        }
    }

    /// Constructs a circle/arc entity carrying its fixed frame.
    pub fn with_frame(kind: SpaceEntityKind, param_start: usize, frame: SpaceFrame3) -> Self {
        Self {
            kind,
            param_start,
            frame: Some(frame),
            control_point_count: 0,
        }
    }

    /// Constructs a spline entity with a given control point count.
    pub fn with_control_points(kind: SpaceEntityKind, param_start: usize, count: usize) -> Self {
        Self {
            kind,
            param_start,
            frame: None,
            control_point_count: count,
        }
    }

    /// Per-entity parameter count (for [`SpaceEntityKind::SpaceSpline`])
    /// which depends on the control point count).
    #[inline]
    pub fn param_count(&self) -> usize {
        match self.kind {
            SpaceEntityKind::SpaceSpline => self.control_point_count * 3,
            other => other.param_count(),
        }
    }

    /// Absolute parameter index for a relative offset within this entity.
    #[inline]
    pub fn param(&self, offset: usize) -> usize {
        self.param_start + offset
    }

    /// The fixed frame of a circle/arc entity (the canonical XY frame when a
    /// hand-built entity omits it).
    #[inline]
    pub fn curve_frame(&self) -> SpaceFrame3 {
        self.frame.unwrap_or(SpaceFrame3::XY)
    }
}

/// A reference to a specific 3D point within an entity.
///
/// Used by constraints that operate on points.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpacePointRef {
    /// A `SpacePoint` entity — the entity itself is the point.
    Point(SpaceEntityId),
    /// The start endpoint of a `SpaceLine` entity (params `[0,1,2]`).
    LineStart(SpaceEntityId),
    /// The end endpoint of a `SpaceLine` entity (params `[3,4,5]`).
    LineEnd(SpaceEntityId),
    /// The center of a `Sphere` entity (params `[0,1,2]`).
    SphereCenter(SpaceEntityId),
    /// The center of a `SpaceCircle` entity (params `[0,1,2]`).
    CircleCenter(SpaceEntityId),
    /// The center of a `SpaceArc` entity (params `[0,1,2]`).
    ArcCenter(SpaceEntityId),
}

impl SpacePointRef {
    /// Return the (x_param_idx, y_param_idx, z_param_idx) absolute indices
    /// into the sketch parameter vector.
    pub fn param_indices(&self, entities: &[SpaceEntity]) -> (usize, usize, usize) {
        match *self {
            SpacePointRef::Point(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1), e.param(2))
            }
            SpacePointRef::LineStart(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1), e.param(2))
            }
            SpacePointRef::LineEnd(id) => {
                let e = &entities[id];
                (e.param(3), e.param(4), e.param(5))
            }
            SpacePointRef::SphereCenter(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1), e.param(2))
            }
            SpacePointRef::CircleCenter(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1), e.param(2))
            }
            SpacePointRef::ArcCenter(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1), e.param(2))
            }
        }
    }
}

impl From<SpaceEntityId> for SpacePointRef {
    fn from(id: SpaceEntityId) -> Self {
        SpacePointRef::Point(id)
    }
}
