//! Geometric entity types and their parameter layout.
//!
//! All entity parameters are stored in a flat `Vec<f64>` owned by [`crate::Sketch`].
//! Each entity records its `param_start` offset into that vector.
//!
//! | Kind   | Params (relative offsets)                          |
//! |--------|----------------------------------------------------|
//! | Point  | `[0]`=x, `[1]`=y                                       |
//! | Line   | `[0]`=x1, `[1]`=y1, `[2]`=x2, `[3]`=y2                    |
//! | Circle | `[0]`=cx, `[1]`=cy, `[2]`=r                              |
//! | Arc    | `[0]`=cx, `[1]`=cy, `[2]`=r, `[3]`=start_angle, `[4]`=end_angle |

/// Index into the entity list of a [`crate::Sketch`].
pub type EntityId = usize;

/// The geometric kind of an entity.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntityKind {
    Point,
    Line,
    Circle,
    Arc,
}

impl EntityKind {
    /// Number of scalar parameters this entity type uses.
    pub fn param_count(self) -> usize {
        match self {
            EntityKind::Point => 2,
            EntityKind::Line => 4,
            EntityKind::Circle => 3,
            EntityKind::Arc => 5,
        }
    }
}

/// Metadata for a single entity stored in a [`crate::Sketch`].
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone)]
pub struct Entity {
    pub kind: EntityKind,
    /// Index of the first parameter in the sketch's flat `params` vector.
    pub param_start: usize,
}

impl Entity {
    pub fn new(kind: EntityKind, param_start: usize) -> Self {
        Self { kind, param_start }
    }

    /// Absolute parameter index for a relative offset within this entity.
    #[inline]
    pub fn param(&self, offset: usize) -> usize {
        self.param_start + offset
    }
}

/// A reference to a specific 2D point within an entity.
///
/// Used by constraints that operate on points (e.g. Coincident, Fixed).
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointRef {
    /// A `Point` entity — the entity itself is the point.
    Point(EntityId),
    /// The start endpoint of a `Line` entity (params `[0,1]`).
    LineStart(EntityId),
    /// The end endpoint of a `Line` entity (params `[2,3]`).
    LineEnd(EntityId),
    /// The center of a `Circle` or `Arc` entity (params `[0,1]`).
    Center(EntityId),
    /// The start endpoint of an `Arc` entity (derived from center + radius + start_angle).
    ArcStart(EntityId),
    /// The end endpoint of an `Arc` entity (derived from center + radius + end_angle).
    ArcEnd(EntityId),
}

impl PointRef {
    /// The entity this point belongs to.
    pub fn entity(&self) -> EntityId {
        match *self {
            PointRef::Point(id)
            | PointRef::LineStart(id)
            | PointRef::LineEnd(id)
            | PointRef::Center(id)
            | PointRef::ArcStart(id)
            | PointRef::ArcEnd(id) => id,
        }
    }

    /// The entity kinds this kind of reference may point into.
    pub fn expected_kinds(&self) -> &'static [EntityKind] {
        match self {
            PointRef::Point(_) => &[EntityKind::Point],
            PointRef::LineStart(_) | PointRef::LineEnd(_) => &[EntityKind::Line],
            PointRef::Center(_) => &[EntityKind::Circle, EntityKind::Arc],
            PointRef::ArcStart(_) | PointRef::ArcEnd(_) => &[EntityKind::Arc],
        }
    }

    /// Return the (x_param_idx, y_param_idx) absolute indices into the sketch
    /// parameter vector. For derived points (arc endpoints), returns the center
    /// indices as a fallback; use `point_coords()` for the actual coordinates.
    pub fn param_indices(&self, entities: &[Entity]) -> (usize, usize) {
        match *self {
            PointRef::Point(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1))
            }
            PointRef::LineStart(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1))
            }
            PointRef::LineEnd(id) => {
                let e = &entities[id];
                (e.param(2), e.param(3))
            }
            PointRef::Center(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1))
            }
            PointRef::ArcStart(id) | PointRef::ArcEnd(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1))
            }
        }
    }

    /// Compute the actual (x, y) coordinates of this point reference.
    /// For stored points (Point, LineStart, LineEnd, Center), reads directly
    /// from the parameter vector. For derived points (ArcStart, ArcEnd),
    /// computes center + radius * (cos(angle), sin(angle)).
    pub fn point_coords(&self, params: &[f64], entities: &[Entity]) -> (f64, f64) {
        match *self {
            PointRef::Point(id) => {
                let e = &entities[id];
                (params[e.param(0)], params[e.param(1)])
            }
            PointRef::LineStart(id) => {
                let e = &entities[id];
                (params[e.param(0)], params[e.param(1)])
            }
            PointRef::LineEnd(id) => {
                let e = &entities[id];
                (params[e.param(2)], params[e.param(3)])
            }
            PointRef::Center(id) => {
                let e = &entities[id];
                (params[e.param(0)], params[e.param(1)])
            }
            PointRef::ArcStart(id) => {
                let e = &entities[id];
                let cx = params[e.param(0)];
                let cy = params[e.param(1)];
                let r = params[e.param(2)];
                let angle = params[e.param(3)];
                (cx + r * angle.cos(), cy + r * angle.sin())
            }
            PointRef::ArcEnd(id) => {
                let e = &entities[id];
                let cx = params[e.param(0)];
                let cy = params[e.param(1)];
                let r = params[e.param(2)];
                let angle = params[e.param(4)];
                (cx + r * angle.cos(), cy + r * angle.sin())
            }
        }
    }
}
