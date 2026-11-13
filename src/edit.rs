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

/// A defect found by [`Sketch::validate`][crate::Sketch::validate].
#[derive(Debug, Clone, PartialEq)]
pub enum SketchError {
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
        found: EntityKind,
        /// The kinds the reference accepts.
        expected: &'static [EntityKind],
    },
    /// A dimensional constraint whose value is NaN or infinite.
    NonFiniteValue {
        /// Index of the constraint.
        constraint: usize,
    },
    /// A parameter that is NaN or infinite.
    NonFiniteParameter {
        /// Index into the parameter vector.
        index: usize,
    },
}

impl std::fmt::Display for SketchError {
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
        }
    }
}

impl std::error::Error for SketchError {}
