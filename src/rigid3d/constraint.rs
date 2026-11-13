//! Constraint definitions for the rigid-pose solver family.

use glam::{DQuat, DVec3};

use super::Point3;
use crate::edit::ConstraintState;

/// Unique identifier for a component in the assembly.
pub type ComponentId = usize;

/// The kind of a constraint, without its data: the tag of a
/// [`ConstraintKind`]. Every kind has solver equations — placeholder kinds
/// without equations are deliberately not modeled (a host implements them
/// as combinations of these, or waits until they get real equations).
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ConstraintType {
    /// Fix component in place (3 DOF removed, 6 with an orientation)
    Fixed,
    /// Coincident points (3 DOF removed)
    Coincident,
    /// Distance between points (1 DOF removed)
    Distance,
    /// Angle between vectors (1 DOF removed)
    Angle,
    /// Parallel vectors (2 DOF removed)
    Parallel,
    /// Perpendicular vectors (1 DOF removed)
    Perpendicular,
    /// Concentric axes: two body axes collinear (4 DOF removed)
    Concentric,
    /// A point lies on a plane (1 DOF removed)
    PointOnPlane,
    /// Lock the relative pose of two bodies (6 DOF removed)
    Lock,
}

/// What a constraint asks for, with exactly the data that takes.
///
/// Points and directions are in the local frame of the body they belong to;
/// the solver carries them to world space with the body's pose. Directions
/// need not be unit length: they are normalized where they are used.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConstraintKind {
    /// The body's origin is at `position`, and, when `orientation` is
    /// given, the body is turned to it. With an orientation the body is
    /// pinned in full and the solver treats it as ground.
    Fixed {
        body: ComponentId,
        /// World position of the body origin.
        position: Point3,
        /// World orientation of the body (`None`: only the position is
        /// held).
        orientation: Option<DQuat>,
    },
    /// `point1` on `body1` and `point2` on `body2` are the same world point.
    Coincident {
        body1: ComponentId,
        point1: Point3,
        body2: ComponentId,
        point2: Point3,
    },
    /// `point1` on `body1` and `point2` on `body2` are `distance` apart.
    Distance {
        body1: ComponentId,
        point1: Point3,
        body2: ComponentId,
        point2: Point3,
        distance: f64,
    },
    /// The directions `dir1` on `body1` and `dir2` on `body2` enclose
    /// `angle` (radians).
    Angle {
        body1: ComponentId,
        dir1: DVec3,
        body2: ComponentId,
        dir2: DVec3,
        angle: f64,
    },
    /// The directions `dir1` on `body1` and `dir2` on `body2` are parallel
    /// (either way round).
    Parallel {
        body1: ComponentId,
        dir1: DVec3,
        body2: ComponentId,
        dir2: DVec3,
    },
    /// The directions `dir1` on `body1` and `dir2` on `body2` are at right
    /// angles.
    Perpendicular {
        body1: ComponentId,
        dir1: DVec3,
        body2: ComponentId,
        dir2: DVec3,
    },
    /// The axis through `point1` along `dir1` on `body1` and the axis
    /// through `point2` along `dir2` on `body2` are one line — the
    /// shaft-in-bore mate.
    Concentric {
        body1: ComponentId,
        point1: Point3,
        dir1: DVec3,
        body2: ComponentId,
        point2: Point3,
        dir2: DVec3,
    },
    /// `point` on `body2` lies on the plane through `plane_point` with
    /// normal `plane_normal`, which rides on `body1`.
    PointOnPlane {
        body1: ComponentId,
        plane_point: Point3,
        plane_normal: DVec3,
        body2: ComponentId,
        point: Point3,
    },
    /// The pose of `body2` is the pose of `body1` composed with a fixed
    /// relative transform, `translation` and `rotation` both expressed in
    /// the frame of `body1`.
    Lock {
        body1: ComponentId,
        body2: ComponentId,
        translation: DVec3,
        rotation: DQuat,
    },
}

impl ConstraintKind {
    /// The kind without its data.
    pub fn constraint_type(&self) -> ConstraintType {
        match self {
            ConstraintKind::Fixed { .. } => ConstraintType::Fixed,
            ConstraintKind::Coincident { .. } => ConstraintType::Coincident,
            ConstraintKind::Distance { .. } => ConstraintType::Distance,
            ConstraintKind::Angle { .. } => ConstraintType::Angle,
            ConstraintKind::Parallel { .. } => ConstraintType::Parallel,
            ConstraintKind::Perpendicular { .. } => ConstraintType::Perpendicular,
            ConstraintKind::Concentric { .. } => ConstraintType::Concentric,
            ConstraintKind::PointOnPlane { .. } => ConstraintType::PointOnPlane,
            ConstraintKind::Lock { .. } => ConstraintType::Lock,
        }
    }

    /// The bodies the constraint acts on: the first, and the second if it
    /// has one (`Fixed` acts on a single body).
    pub fn bodies(&self) -> (ComponentId, Option<ComponentId>) {
        match *self {
            ConstraintKind::Fixed { body, .. } => (body, None),
            ConstraintKind::Coincident { body1, body2, .. }
            | ConstraintKind::Distance { body1, body2, .. }
            | ConstraintKind::Angle { body1, body2, .. }
            | ConstraintKind::Parallel { body1, body2, .. }
            | ConstraintKind::Perpendicular { body1, body2, .. }
            | ConstraintKind::Concentric { body1, body2, .. }
            | ConstraintKind::PointOnPlane { body1, body2, .. }
            | ConstraintKind::Lock { body1, body2, .. } => (body1, Some(body2)),
        }
    }

    /// The same constraint with every body id passed through `f`.
    pub fn map_bodies(mut self, mut f: impl FnMut(ComponentId) -> ComponentId) -> Self {
        match &mut self {
            ConstraintKind::Fixed { body, .. } => *body = f(*body),
            ConstraintKind::Coincident { body1, body2, .. }
            | ConstraintKind::Distance { body1, body2, .. }
            | ConstraintKind::Angle { body1, body2, .. }
            | ConstraintKind::Parallel { body1, body2, .. }
            | ConstraintKind::Perpendicular { body1, body2, .. }
            | ConstraintKind::Concentric { body1, body2, .. }
            | ConstraintKind::PointOnPlane { body1, body2, .. }
            | ConstraintKind::Lock { body1, body2, .. } => {
                *body1 = f(*body1);
                *body2 = f(*body2);
            }
        }
        self
    }

    /// Number of degrees of freedom the constraint removes when nothing
    /// else interferes.
    pub fn dof_removed(&self) -> usize {
        match self {
            ConstraintKind::Fixed {
                orientation: Some(_),
                ..
            } => 6,
            ConstraintKind::Fixed { .. } => 3,
            ConstraintKind::Coincident { .. } => 3,
            ConstraintKind::Distance { .. } => 1,
            ConstraintKind::Angle { .. } => 1,
            ConstraintKind::Parallel { .. } => 2,
            ConstraintKind::Perpendicular { .. } => 1,
            ConstraintKind::Concentric { .. } => 4,
            ConstraintKind::PointOnPlane { .. } => 1,
            ConstraintKind::Lock { .. } => 6,
        }
    }

    /// The driving value of a dimensional constraint: the distance of a
    /// `Distance`, the angle (radians) of an `Angle`. `None` for the kinds
    /// without a single scalar value.
    pub fn value(&self) -> Option<f64> {
        match *self {
            ConstraintKind::Distance { distance, .. } => Some(distance),
            ConstraintKind::Angle { angle, .. } => Some(angle),
            _ => None,
        }
    }

    /// Set the driving value of a dimensional constraint (see
    /// [`ConstraintKind::value`]). Returns `false`, changing nothing, for a
    /// kind that has none.
    pub fn set_value(&mut self, value: f64) -> bool {
        match self {
            ConstraintKind::Distance { distance, .. } => *distance = value,
            ConstraintKind::Angle { angle, .. } => *angle = value,
            _ => return false,
        }
        true
    }

    /// What a dimensional constraint measures at the given poses: the value
    /// [`ConstraintKind::value`] would have to be for the constraint to
    /// hold as the bodies stand. `None` for the kinds without a value.
    ///
    /// An angle is reported in `[0, π]`: the angle equation compares
    /// cosines, so that is all it distinguishes.
    pub fn measure(&self, positions: &[Point3], rotations: &[DQuat]) -> Option<f64> {
        match *self {
            ConstraintKind::Distance {
                body1,
                point1,
                body2,
                point2,
                ..
            } => {
                let w1 = positions[body1] + rotations[body1] * point1;
                let w2 = positions[body2] + rotations[body2] * point2;
                Some((w1 - w2).length())
            }
            ConstraintKind::Angle {
                body1,
                dir1,
                body2,
                dir2,
                ..
            } => {
                let d1 = (rotations[body1] * dir1).normalize_or_zero();
                let d2 = (rotations[body2] * dir2).normalize_or_zero();
                Some(d1.dot(d2).clamp(-1.0, 1.0).acos())
            }
            _ => None,
        }
    }

    /// Whether every number the constraint carries is finite: points,
    /// directions, the distance or angle, the target pose.
    pub fn has_finite_values(&self) -> bool {
        match *self {
            ConstraintKind::Fixed {
                position,
                orientation,
                ..
            } => position.is_finite() && orientation.is_none_or(|q| q.is_finite()),
            ConstraintKind::Coincident { point1, point2, .. } => {
                point1.is_finite() && point2.is_finite()
            }
            ConstraintKind::Distance {
                point1,
                point2,
                distance,
                ..
            } => point1.is_finite() && point2.is_finite() && distance.is_finite(),
            ConstraintKind::Angle {
                dir1, dir2, angle, ..
            } => dir1.is_finite() && dir2.is_finite() && angle.is_finite(),
            ConstraintKind::Parallel { dir1, dir2, .. }
            | ConstraintKind::Perpendicular { dir1, dir2, .. } => {
                dir1.is_finite() && dir2.is_finite()
            }
            ConstraintKind::Concentric {
                point1,
                dir1,
                point2,
                dir2,
                ..
            } => point1.is_finite() && dir1.is_finite() && point2.is_finite() && dir2.is_finite(),
            ConstraintKind::PointOnPlane {
                plane_point,
                plane_normal,
                point,
                ..
            } => plane_point.is_finite() && plane_normal.is_finite() && point.is_finite(),
            ConstraintKind::Lock {
                translation,
                rotation,
                ..
            } => translation.is_finite() && rotation.is_finite(),
        }
    }

    /// Whether every direction and rotation the constraint carries has a
    /// length: a zero direction does not say which way, and a zero
    /// quaternion is no rotation.
    pub fn has_usable_directions(&self) -> bool {
        const TINY: f64 = 1.0e-14;
        let direction = |v: DVec3| v.length() > TINY;
        let rotation = |q: DQuat| q.length() > TINY;
        match *self {
            ConstraintKind::Fixed { orientation, .. } => orientation.is_none_or(rotation),
            ConstraintKind::Angle { dir1, dir2, .. }
            | ConstraintKind::Parallel { dir1, dir2, .. }
            | ConstraintKind::Perpendicular { dir1, dir2, .. }
            | ConstraintKind::Concentric { dir1, dir2, .. } => direction(dir1) && direction(dir2),
            ConstraintKind::PointOnPlane { plane_normal, .. } => direction(plane_normal),
            ConstraintKind::Lock { rotation: q, .. } => rotation(q),
            ConstraintKind::Coincident { .. } | ConstraintKind::Distance { .. } => true,
        }
    }
}

/// A constraint of an assembly: what it asks for, under the id the host
/// knows it by.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Constraint {
    /// The host's identifier for this constraint (a mate id). The diagnosis
    /// reports by it. Ids need not be unique or dense.
    pub id: usize,
    /// What the constraint asks for.
    pub kind: ConstraintKind,
}

impl Constraint {
    /// A constraint of the given kind.
    pub fn new(id: usize, kind: ConstraintKind) -> Self {
        Self { id, kind }
    }

    /// Hold a component's origin at `point` (its orientation stays free).
    pub fn fixed(id: usize, component: ComponentId, point: Point3) -> Self {
        Self::new(
            id,
            ConstraintKind::Fixed {
                body: component,
                position: point,
                orientation: None,
            },
        )
    }

    /// Ground a component at a pose (position + world orientation quaternion).
    pub fn fixed_with_orientation(
        id: usize,
        component: ComponentId,
        position: Point3,
        orientation: DQuat,
    ) -> Self {
        Self::new(
            id,
            ConstraintKind::Fixed {
                body: component,
                position,
                orientation: Some(orientation.normalize()),
            },
        )
    }

    /// Create a coincident constraint.
    pub fn coincident(
        id: usize,
        comp1: ComponentId,
        pt1: Point3,
        comp2: ComponentId,
        pt2: Point3,
    ) -> Self {
        Self::new(
            id,
            ConstraintKind::Coincident {
                body1: comp1,
                point1: pt1,
                body2: comp2,
                point2: pt2,
            },
        )
    }

    /// Create a distance constraint.
    pub fn distance(
        id: usize,
        comp1: ComponentId,
        pt1: Point3,
        comp2: ComponentId,
        pt2: Point3,
        dist: f64,
    ) -> Self {
        Self::new(
            id,
            ConstraintKind::Distance {
                body1: comp1,
                point1: pt1,
                body2: comp2,
                point2: pt2,
                distance: dist,
            },
        )
    }

    /// Create an angle constraint between two directions (in component
    /// space).
    pub fn angle(
        id: usize,
        comp1: ComponentId,
        dir1: DVec3,
        comp2: ComponentId,
        dir2: DVec3,
        ang: f64,
    ) -> Self {
        Self::new(
            id,
            ConstraintKind::Angle {
                body1: comp1,
                dir1,
                body2: comp2,
                dir2,
                angle: ang,
            },
        )
    }

    /// Create a parallel constraint between two directions (in component
    /// space).
    pub fn parallel(
        id: usize,
        comp1: ComponentId,
        dir1: DVec3,
        comp2: ComponentId,
        dir2: DVec3,
    ) -> Self {
        Self::new(
            id,
            ConstraintKind::Parallel {
                body1: comp1,
                dir1,
                body2: comp2,
                dir2,
            },
        )
    }

    /// Create a perpendicular constraint between two directions (in
    /// component space).
    pub fn perpendicular(
        id: usize,
        comp1: ComponentId,
        dir1: DVec3,
        comp2: ComponentId,
        dir2: DVec3,
    ) -> Self {
        Self::new(
            id,
            ConstraintKind::Perpendicular {
                body1: comp1,
                dir1,
                body2: comp2,
                dir2,
            },
        )
    }

    /// Create a concentric-axes constraint: the two body axes (a point on
    /// the axis plus the axis direction, both in component space) must be
    /// collinear — the shaft/bore mate (4 DOF removed).
    pub fn concentric(
        id: usize,
        comp1: ComponentId,
        axis_point1: Point3,
        axis_dir1: DVec3,
        comp2: ComponentId,
        axis_point2: Point3,
        axis_dir2: DVec3,
    ) -> Self {
        Self::new(
            id,
            ConstraintKind::Concentric {
                body1: comp1,
                point1: axis_point1,
                dir1: axis_dir1,
                body2: comp2,
                point2: axis_point2,
                dir2: axis_dir2,
            },
        )
    }

    /// Create a point-on-plane constraint: the plane is carried by the
    /// FIRST body (`plane_point` + `plane_normal`), the point is `point` on
    /// the SECOND body. One DOF removed.
    pub fn point_on_plane(
        id: usize,
        comp1: ComponentId,
        plane_point: Point3,
        plane_normal: DVec3,
        comp2: ComponentId,
        point: Point3,
    ) -> Self {
        Self::new(
            id,
            ConstraintKind::PointOnPlane {
                body1: comp1,
                plane_point,
                plane_normal,
                body2: comp2,
                point,
            },
        )
    }

    /// Create a relative-pose lock: the second body's pose equals the first
    /// body's pose composed with `lock_translation` / `lock_rotation`
    /// (both expressed in the FIRST body's frame). Six DOF removed.
    pub fn lock(
        id: usize,
        comp1: ComponentId,
        comp2: ComponentId,
        lock_translation: DVec3,
        lock_rotation: DQuat,
    ) -> Self {
        Self::new(
            id,
            ConstraintKind::Lock {
                body1: comp1,
                body2: comp2,
                translation: lock_translation,
                rotation: lock_rotation.normalize(),
            },
        )
    }

    /// The kind of the constraint, without its data.
    pub fn constraint_type(&self) -> ConstraintType {
        self.kind.constraint_type()
    }

    /// The first (or only) component the constraint acts on.
    pub fn component1(&self) -> ComponentId {
        self.kind.bodies().0
    }

    /// The second component the constraint acts on (`None` for `Fixed`).
    pub fn component2(&self) -> Option<ComponentId> {
        self.kind.bodies().1
    }

    /// Get the number of DOF removed by this constraint.
    pub fn dof_removed(&self) -> usize {
        self.kind.dof_removed()
    }
}

/// A defect found by [`ConstraintSystem::validate`]. `constraint` is the
/// index of the constraint in [`ConstraintSystem::constraints`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyError {
    /// A constraint refers to a component the system does not have.
    ComponentOutOfRange {
        constraint: usize,
        component: ComponentId,
    },
    /// A constraint refers to a component that has been removed.
    ComponentRemoved {
        constraint: usize,
        component: ComponentId,
    },
    /// A point, direction, distance, angle or rotation is NaN or infinite.
    NonFiniteValue { constraint: usize },
    /// A direction or a rotation of zero length.
    ZeroLength { constraint: usize },
    /// Both sides of the constraint are on one component, so moving the
    /// component cannot change whether it holds.
    SameComponent {
        constraint: usize,
        component: ComponentId,
    },
}

impl std::fmt::Display for AssemblyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AssemblyError::ComponentOutOfRange {
                constraint,
                component,
            } => write!(
                f,
                "constraint {constraint} refers to component {component}, which does not exist"
            ),
            AssemblyError::ComponentRemoved {
                constraint,
                component,
            } => write!(
                f,
                "constraint {constraint} refers to component {component}, which has been removed"
            ),
            AssemblyError::NonFiniteValue { constraint } => {
                write!(f, "constraint {constraint} has a value that is not finite")
            }
            AssemblyError::ZeroLength { constraint } => write!(
                f,
                "constraint {constraint} has a direction or rotation of zero length"
            ),
            AssemblyError::SameComponent {
                constraint,
                component,
            } => write!(
                f,
                "constraint {constraint} has component {component} on both sides"
            ),
        }
    }
}

impl std::error::Error for AssemblyError {}

/// A system of constraints.
///
/// Constraints are addressed by their position in
/// [`ConstraintSystem::constraints`], and components by their id; both stay
/// put when something is removed (the slot is left behind, marked). Edit
/// through the methods: pushing into or deleting from the public list
/// directly goes unnoticed by [`ConstraintSystem::topology_revision`], and
/// deleting from it shifts indices.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ConstraintSystem {
    /// Constraint slots, including suppressed and removed ones (see
    /// [`ConstraintSystem::constraint_state`]).
    pub constraints: Vec<Constraint>,
    /// Number of component slots, removed ones included.
    pub num_components: usize,
    /// The state of each constraint slot (a missing entry means active).
    #[cfg_attr(feature = "serde", serde(default))]
    states: Vec<ConstraintState>,
    /// Per component: removed (a missing entry means present).
    #[cfg_attr(feature = "serde", serde(default))]
    removed: Vec<bool>,
    /// See [`ConstraintSystem::topology_revision`].
    #[cfg_attr(feature = "serde", serde(default))]
    revision: u64,
}

impl ConstraintSystem {
    /// Create an empty system.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a system with given number of components.
    pub fn with_components(num_components: usize) -> Self {
        Self {
            num_components,
            ..Self::default()
        }
    }

    /// Add a component. Returns its id.
    pub fn add_component(&mut self) -> ComponentId {
        self.num_components += 1;
        self.revision += 1;
        self.num_components - 1
    }

    /// Add a constraint. Returns its index, which stays valid for the life
    /// of the system: removing other constraints does not shift it.
    pub fn add_constraint(&mut self, constraint: Constraint) -> usize {
        self.constraints.push(constraint);
        self.states
            .resize(self.constraints.len(), ConstraintState::Active);
        self.revision += 1;
        self.constraints.len() - 1
    }

    /// Whether constraint `index` is enforced, suppressed or removed.
    pub fn constraint_state(&self, index: usize) -> ConstraintState {
        self.states.get(index).copied().unwrap_or_default()
    }

    /// Whether constraint `index` is enforced.
    pub fn is_active(&self, index: usize) -> bool {
        self.constraint_state(index) == ConstraintState::Active
    }

    /// The enforced constraints, each with its index.
    pub fn active_constraints(&self) -> impl Iterator<Item = (usize, &Constraint)> {
        self.constraints
            .iter()
            .enumerate()
            .filter(|(index, _)| self.is_active(*index))
    }

    /// Delete constraint `index`. Its slot stays (so other indices do not
    /// shift) but it is never enforced again. Returns `false` if it was
    /// already removed.
    pub fn remove_constraint(&mut self, index: usize) -> bool {
        self.set_state(index, ConstraintState::Removed)
    }

    /// Enforce constraint `index` (`true`) or keep it without enforcing it
    /// (`false`). A suppressed dimension is a reference dimension: nothing
    /// drives it, and [`ConstraintSystem::measure`] reads what it currently
    /// is. Returns `false`, changing nothing, for a removed constraint.
    pub fn set_constraint_enabled(&mut self, index: usize, enabled: bool) -> bool {
        if self.constraint_state(index) == ConstraintState::Removed {
            return false;
        }
        self.set_state(
            index,
            if enabled {
                ConstraintState::Active
            } else {
                ConstraintState::Suppressed
            },
        );
        true
    }

    /// Put a different constraint in slot `index` and make it active.
    pub fn replace_constraint(&mut self, index: usize, constraint: Constraint) {
        self.constraints[index] = constraint;
        self.states
            .resize(self.constraints.len(), ConstraintState::Active);
        self.states[index] = ConstraintState::Active;
        self.revision += 1;
    }

    /// The driving value of dimensional constraint `index` (see
    /// [`ConstraintKind::value`]).
    pub fn constraint_value(&self, index: usize) -> Option<f64> {
        self.constraints[index].kind.value()
    }

    /// Change the driving value of dimensional constraint `index` in place.
    /// Nothing moves until the next solve; solving from the poses of the
    /// last one ([`ConstraintSolver::solve_from`][crate::rigid3d::ConstraintSolver::solve_from])
    /// continues from where the assembly stands. Returns `false` for a
    /// constraint without a scalar value.
    ///
    /// This does not change which bodies any constraint joins, so
    /// [`ConstraintSystem::topology_revision`] stays the same.
    pub fn set_constraint_value(&mut self, index: usize, value: f64) -> bool {
        self.constraints[index].kind.set_value(value)
    }

    /// What dimensional constraint `index` measures at the given poses,
    /// whatever its state (see [`ConstraintKind::measure`]).
    pub fn measure(&self, index: usize, positions: &[Point3], rotations: &[DQuat]) -> Option<f64> {
        self.constraints[index].kind.measure(positions, rotations)
    }

    fn set_state(&mut self, index: usize, state: ConstraintState) -> bool {
        assert!(
            index < self.constraints.len(),
            "constraint index {index} out of range"
        );
        self.states
            .resize(self.constraints.len(), ConstraintState::Active);
        if self.states[index] == state {
            return false;
        }
        self.states[index] = state;
        self.revision += 1;
        true
    }

    /// Delete a component together with every constraint that refers to
    /// it. Returns the indices of the constraints removed with it.
    ///
    /// The component's slot stays (so other ids do not shift, and the pose
    /// arrays keep their length); it is no longer solved for and no longer
    /// counts as six freedoms.
    pub fn remove_component(&mut self, id: ComponentId) -> Vec<usize> {
        assert!(id < self.num_components, "component {id} out of range");
        self.removed.resize(self.num_components, false);
        if self.removed[id] {
            return Vec::new();
        }
        self.removed[id] = true;
        self.revision += 1;

        let dependents: Vec<usize> = (0..self.constraints.len())
            .filter(|&index| self.constraint_state(index) != ConstraintState::Removed)
            .filter(|&index| {
                let (first, second) = self.constraints[index].kind.bodies();
                first == id || second == Some(id)
            })
            .collect();
        for &index in &dependents {
            self.set_state(index, ConstraintState::Removed);
        }
        dependents
    }

    /// Whether component `id` has been removed.
    pub fn is_component_removed(&self, id: ComponentId) -> bool {
        self.removed.get(id).copied().unwrap_or(false)
    }

    /// A counter that changes whenever the structure of the system does: a
    /// component or constraint added or removed, a constraint enabled,
    /// suppressed or replaced. Changing a dimension value leaves it alone,
    /// so anything derived from the structure alone (the clusters, for one)
    /// can be kept while it holds.
    pub fn topology_revision(&self) -> u64 {
        self.revision
    }

    /// Check the model for defects that would make a solve meaningless or
    /// panic: constraints referring to components that do not exist or
    /// were removed, non-finite numbers, directions or rotations of zero
    /// length, and constraints with one component on both sides. Removed
    /// constraints are not checked.
    pub fn validate(&self) -> Result<(), Vec<AssemblyError>> {
        let mut errors = Vec::new();
        for (constraint, c) in self.constraints.iter().enumerate() {
            if self.constraint_state(constraint) == ConstraintState::Removed {
                continue;
            }
            let (first, second) = c.kind.bodies();
            for component in std::iter::once(first).chain(second) {
                if component >= self.num_components {
                    errors.push(AssemblyError::ComponentOutOfRange {
                        constraint,
                        component,
                    });
                } else if self.is_component_removed(component) {
                    errors.push(AssemblyError::ComponentRemoved {
                        constraint,
                        component,
                    });
                }
            }
            if second == Some(first) {
                errors.push(AssemblyError::SameComponent {
                    constraint,
                    component: first,
                });
            }
            if !c.kind.has_finite_values() {
                errors.push(AssemblyError::NonFiniteValue { constraint });
            } else if !c.kind.has_usable_directions() {
                errors.push(AssemblyError::ZeroLength { constraint });
            }
        }
        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }

    /// The enforced constraints involving a specific component.
    pub fn constraints_for_component(&self, component: ComponentId) -> Vec<&Constraint> {
        self.active_constraints()
            .map(|(_, c)| c)
            .filter(|c| c.component1() == component || c.component2() == Some(component))
            .collect()
    }

    /// The freedoms the enforced constraints would remove if none of them
    /// overlapped: a count, not a rank. The diagnosis
    /// ([`diagnose_assembly`][crate::rigid3d::diagnose_assembly]) gives the
    /// number actually removed.
    pub fn total_dof_removed(&self) -> usize {
        self.active_constraints()
            .map(|(_, c)| c.dof_removed())
            .sum()
    }

    /// Get the maximum possible DOF (6 per component that has not been
    /// removed).
    pub fn max_dof(&self) -> usize {
        self.live_components() * 6
    }

    /// Number of components that have not been removed.
    pub fn live_components(&self) -> usize {
        (0..self.num_components)
            .filter(|&id| !self.is_component_removed(id))
            .count()
    }

    /// Per column of the pose parameterization (six per component slot):
    /// held, because its component has been removed.
    pub(crate) fn held_columns(&self) -> Vec<bool> {
        (0..self.num_components)
            .flat_map(|id| [self.is_component_removed(id); 6])
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fixed_constraint() {
        let c = Constraint::fixed(0, 1, DVec3::ZERO);
        assert_eq!(c.constraint_type(), ConstraintType::Fixed);
        assert_eq!(c.dof_removed(), 3);
        assert_eq!((c.component1(), c.component2()), (1, None));
    }

    #[test]
    fn test_fixed_with_orientation_dof() {
        let c = Constraint::fixed_with_orientation(0, 0, DVec3::ZERO, DQuat::IDENTITY);
        assert_eq!(c.dof_removed(), 6);
        assert!(matches!(
            c.kind,
            ConstraintKind::Fixed {
                orientation: Some(_),
                ..
            }
        ));
    }

    #[test]
    fn test_coincident_constraint() {
        let c = Constraint::coincident(0, 1, DVec3::ZERO, 2, DVec3::new(1.0, 0.0, 0.0));
        assert_eq!(c.constraint_type(), ConstraintType::Coincident);
        assert_eq!(c.dof_removed(), 3);
        assert_eq!((c.component1(), c.component2()), (1, Some(2)));
    }

    #[test]
    fn test_constraint_system() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 1, DVec3::ZERO));
        system.add_constraint(Constraint::coincident(
            1,
            1,
            DVec3::ZERO,
            2,
            DVec3::new(1.0, 0.0, 0.0),
        ));

        assert_eq!(system.constraints.len(), 2);
        assert_eq!(system.max_dof(), 12);
        assert_eq!(system.constraints_for_component(2).len(), 1);
    }

    /// Every kind names its bodies, and renumbering reaches all of them.
    #[test]
    fn bodies_are_named_and_renumbered() {
        let d = DVec3::X;
        let kinds = [
            Constraint::fixed(0, 3, d).kind,
            Constraint::fixed_with_orientation(0, 3, d, DQuat::IDENTITY).kind,
            Constraint::coincident(0, 3, d, 5, d).kind,
            Constraint::distance(0, 3, d, 5, d, 1.0).kind,
            Constraint::angle(0, 3, d, 5, d, 1.0).kind,
            Constraint::parallel(0, 3, d, 5, d).kind,
            Constraint::perpendicular(0, 3, d, 5, d).kind,
            Constraint::concentric(0, 3, d, d, 5, d, d).kind,
            Constraint::point_on_plane(0, 3, d, d, 5, d).kind,
            Constraint::lock(0, 3, 5, d, DQuat::IDENTITY).kind,
        ];
        for kind in kinds {
            let (first, second) = kind.bodies();
            assert_eq!(first, 3, "{kind:?}");
            assert!(second.is_none() || second == Some(5), "{kind:?}");
            let moved = kind.map_bodies(|body| body + 10);
            assert_eq!(moved.bodies(), (13, second.map(|_| 15)), "{kind:?}");
            assert_eq!(moved.constraint_type(), kind.constraint_type());
            assert_eq!(moved.dof_removed(), kind.dof_removed());
        }
    }
}
