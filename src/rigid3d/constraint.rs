//! Constraint definitions for the rigid-pose solver family.

use glam::{DQuat, DVec3};

use super::Point3;

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

/// A system of constraints.
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ConstraintSystem {
    /// All constraints in the system
    pub constraints: Vec<Constraint>,
    /// Number of components
    pub num_components: usize,
}

impl ConstraintSystem {
    /// Create an empty system.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a system with given number of components.
    pub fn with_components(num_components: usize) -> Self {
        Self {
            constraints: Vec::new(),
            num_components,
        }
    }

    /// Add a constraint.
    pub fn add_constraint(&mut self, constraint: Constraint) {
        self.constraints.push(constraint);
    }

    /// Get constraints involving a specific component.
    pub fn constraints_for_component(&self, component: ComponentId) -> Vec<&Constraint> {
        self.constraints
            .iter()
            .filter(|c| c.component1() == component || c.component2() == Some(component))
            .collect()
    }

    /// Get the total DOF removed by all constraints.
    pub fn total_dof_removed(&self) -> usize {
        self.constraints.iter().map(|c| c.dof_removed()).sum()
    }

    /// Get the maximum possible DOF (6 per component).
    pub fn max_dof(&self) -> usize {
        self.num_components * 6
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
