//! Editing vocabulary shared by the sketch containers: the state of a
//! constraint slot and the errors a model check can report.
//!
//! Entities and constraints are addressed by their position in the lists
//! they were added to. Those positions are stable: removing something
//! leaves its slot in place (marked removed) instead of shifting everything
//! after it, so an index a host holds keeps meaning the same thing.

use crate::entity::EntityKind;

/// Whether a constraint takes part in solving.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ConstraintState {
    /// Enforced by the solver.
    #[default]
    Active,
    /// Kept but not enforced: a reference ("driven") dimension that is only
    /// measured, or a constraint switched off for the moment. It can be
    /// made active again.
    Suppressed,
    /// Deleted. The slot stays so the indices of other constraints do not
    /// shift; it cannot be brought back.
    Removed,
}

/// A defect found by a model check: [`Sketch::validate`][crate::Sketch::validate]
/// or [`SpaceSketch::validate`][crate::space3d::SpaceSketch::validate].
///
/// `K` is the entity-kind type of the family: [`EntityKind`] for the 2D
/// sketch (the default), `SpaceEntityKind` for the 3D one (see
/// [`SpaceSketchError`][crate::space3d::SpaceSketchError]).
#[derive(Debug, Clone, PartialEq)]
pub enum SketchError<K: 'static = EntityKind> {
    /// A constraint refers to an entity index beyond the entity list.
    EntityOutOfRange {
        /// Index of the constraint.
        constraint: usize,
        /// The entity index it refers to.
        entity: usize,
    },
    /// A constraint refers to an entity that has been removed.
    EntityRemoved {
        /// Index of the constraint.
        constraint: usize,
        /// The removed entity.
        entity: usize,
    },
    /// A constraint refers to an entity of a kind it cannot apply to (a
    /// radius on a line, the start point of a circle).
    WrongEntityKind {
        /// Index of the constraint.
        constraint: usize,
        /// The entity referred to.
        entity: usize,
        /// The kind that entity has.
        found: K,
        /// The kinds the reference accepts.
        expected: &'static [K],
    },
    /// A constraint that carries a number which is NaN or infinite: a
    /// dimension value, a target position or a direction.
    NonFiniteValue {
        /// Index of the constraint.
        constraint: usize,
    },
    /// A parameter that is NaN or infinite.
    NonFiniteParameter {
        /// Index into the parameter vector.
        index: usize,
    },
    /// A constraint addresses a control point its spline does not have.
    /// Only the 3D family has splines.
    ControlPointOutOfRange {
        /// Index of the constraint.
        constraint: usize,
        /// The spline entity.
        entity: usize,
        /// The control point asked for.
        index: usize,
        /// How many control points the spline has.
        count: usize,
    },
}

impl<K: std::fmt::Debug> std::fmt::Display for SketchError<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SketchError::EntityOutOfRange { constraint, entity } => write!(
                f,
                "constraint {constraint} refers to entity {entity}, which does not exist"
            ),
            SketchError::EntityRemoved { constraint, entity } => write!(
                f,
                "constraint {constraint} refers to entity {entity}, which has been removed"
            ),
            SketchError::WrongEntityKind {
                constraint,
                entity,
                found,
                expected,
            } => write!(
                f,
                "constraint {constraint} refers to entity {entity}, a {found:?}, where it needs one of {expected:?}"
            ),
            SketchError::NonFiniteValue { constraint } => {
                write!(f, "constraint {constraint} has a value that is not finite")
            }
            SketchError::NonFiniteParameter { index } => {
                write!(f, "parameter {index} is not finite")
            }
            SketchError::ControlPointOutOfRange {
                constraint,
                entity,
                index,
                count,
            } => write!(
                f,
                "constraint {constraint} addresses control point {index} of spline {entity}, which has {count}"
            ),
        }
    }
}

impl<K: std::fmt::Debug> std::error::Error for SketchError<K> {}
