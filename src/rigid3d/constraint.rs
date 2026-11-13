//! Constraint definitions for the rigid-pose solver family.

use glam::{DQuat, DVec3};

use super::Point3;

/// Unique identifier for a component in the assembly.
pub type ComponentId = usize;

/// Types of constraints between components. Every variant has solver
/// equations — placeholder variants without equations are deliberately not
/// modeled (a host implements them as combinations of these, or waits until
/// they get real equations).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstraintType {
    /// Fix component in place (6 DOF removed)
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
    /// Concentric axes: two body axes collinear (4 DOF removed). The mate
    /// pins an axis DIRECTION and an axis POINT per body (e.g. a shaft's
    /// bore axis); residual rows are the direction cross (rank 2) plus the
    /// perpendicular offset of the two axis points (rank 2).
    Concentric,
    /// A point lies on a plane (1 DOF removed). The plane is carried by the
    /// FIRST body ([`Constraint::point1`] + unit [`Constraint::dir1`]); the
    /// point is [`Constraint::point2`] on the SECOND body. Residual row is
    /// the signed perpendicular distance `(w2 - w1) . n1w`.
    PointOnPlane,
    /// Lock the relative pose: the SECOND body's pose equals the FIRST
    /// body's pose composed with a fixed relative transform (6 DOF removed).
    /// The target relative translation/rotation are
    /// [`Constraint::lock_translation`] / [`Constraint::lock_rotation`].
    Lock,
}

/// A constraint between two components.
#[derive(Debug, Clone)]
pub struct Constraint {
    /// Unique constraint identifier
    pub id: usize,
    /// Type of constraint
    pub constraint_type: ConstraintType,
    /// First component
    pub component1: ComponentId,
    /// Reference point on first component
    pub point1: Point3,
    /// Second component (None for fixed constraint)
    pub component2: Option<ComponentId>,
    /// Reference point on second component
    pub point2: Option<Point3>,
    /// Distance value (for distance constraint)
    pub distance: Option<f64>,
    /// Angle value in radians (for angle constraint)
    pub angle: Option<f64>,
    /// World orientation target for [`ConstraintType::Fixed`] (`None` = translation-only fix).
    pub fixed_orientation: Option<DQuat>,
    /// Local axis direction on the first component
    /// ([`ConstraintType::Concentric`] / [`ConstraintType::PointOnPlane`];
    /// unit by authoring contract).
    pub dir1: Option<DVec3>,
    /// Local axis direction on the second component
    /// ([`ConstraintType::Concentric`]; unit by authoring contract).
    pub dir2: Option<DVec3>,
    /// Relative translation of the second body in the first body's frame
    /// ([`ConstraintType::Lock`]).
    pub lock_translation: Option<DVec3>,
    /// Relative rotation of the second body in the first body's frame
    /// ([`ConstraintType::Lock`]).
    pub lock_rotation: Option<DQuat>,
}

impl Constraint {
    /// Create a fixed constraint.
    pub fn fixed(id: usize, component: ComponentId, point: Point3) -> Self {
        Self {
            id,
            constraint_type: ConstraintType::Fixed,
            component1: component,
            point1: point,
            component2: None,
            point2: None,
            distance: None,
            angle: None,
            fixed_orientation: None,
            dir1: None,
            dir2: None,
            lock_translation: None,
            lock_rotation: None,
        }
    }

    /// Ground a component at a pose (position + world orientation quaternion).
    pub fn fixed_with_orientation(
        id: usize,
        component: ComponentId,
        position: Point3,
        orientation: DQuat,
    ) -> Self {
        Self {
            id,
            constraint_type: ConstraintType::Fixed,
            component1: component,
            point1: position,
            component2: None,
            point2: None,
            distance: None,
            angle: None,
            fixed_orientation: Some(orientation.normalize()),
            dir1: None,
            dir2: None,
            lock_translation: None,
            lock_rotation: None,
        }
    }

    /// Create a coincident constraint.
    pub fn coincident(
        id: usize,
        comp1: ComponentId,
        pt1: Point3,
        comp2: ComponentId,
        pt2: Point3,
    ) -> Self {
        Self {
            id,
            constraint_type: ConstraintType::Coincident,
            component1: comp1,
            point1: pt1,
            component2: Some(comp2),
            point2: Some(pt2),
            distance: None,
            angle: None,
            fixed_orientation: None,
            dir1: None,
            dir2: None,
            lock_translation: None,
            lock_rotation: None,
        }
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
        Self {
            id,
            constraint_type: ConstraintType::Distance,
            component1: comp1,
            point1: pt1,
            component2: Some(comp2),
            point2: Some(pt2),
            distance: Some(dist),
            angle: None,
            fixed_orientation: None,
            dir1: None,
            dir2: None,
            lock_translation: None,
            lock_rotation: None,
        }
    }

    /// Create an angle constraint.
    pub fn angle(
        id: usize,
        comp1: ComponentId,
        pt1: Point3,
        comp2: ComponentId,
        pt2: Point3,
        ang: f64,
    ) -> Self {
        Self {
            id,
            constraint_type: ConstraintType::Angle,
            component1: comp1,
            point1: pt1,
            component2: Some(comp2),
            point2: Some(pt2),
            distance: None,
            angle: Some(ang),
            fixed_orientation: None,
            dir1: None,
            dir2: None,
            lock_translation: None,
            lock_rotation: None,
        }
    }

    /// Create a parallel constraint.
    pub fn parallel(
        id: usize,
        comp1: ComponentId,
        pt1: Point3,
        comp2: ComponentId,
        pt2: Point3,
    ) -> Self {
        Self {
            id,
            constraint_type: ConstraintType::Parallel,
            component1: comp1,
            point1: pt1,
            component2: Some(comp2),
            point2: Some(pt2),
            distance: None,
            angle: None,
            fixed_orientation: None,
            dir1: None,
            dir2: None,
            lock_translation: None,
            lock_rotation: None,
        }
    }

    /// Create a perpendicular constraint (direction vectors in component space).
    pub fn perpendicular(
        id: usize,
        comp1: ComponentId,
        pt1: Point3,
        comp2: ComponentId,
        pt2: Point3,
    ) -> Self {
        Self {
            id,
            constraint_type: ConstraintType::Perpendicular,
            component1: comp1,
            point1: pt1,
            component2: Some(comp2),
            point2: Some(pt2),
            distance: None,
            angle: None,
            fixed_orientation: None,
            dir1: None,
            dir2: None,
            lock_translation: None,
            lock_rotation: None,
        }
    }

    /// Create a concentric-axes constraint: the two body axes (a point on
    /// the axis plus the unit axis direction, both in component space) must
    /// be collinear — the shaft/bore mate (4 DOF removed).
    pub fn concentric(
        id: usize,
        comp1: ComponentId,
        axis_point1: Point3,
        axis_dir1: Point3,
        comp2: ComponentId,
        axis_point2: Point3,
        axis_dir2: Point3,
    ) -> Self {
        Self {
            id,
            constraint_type: ConstraintType::Concentric,
            component1: comp1,
            point1: axis_point1,
            component2: Some(comp2),
            point2: Some(axis_point2),
            distance: None,
            angle: None,
            fixed_orientation: None,
            dir1: Some(axis_dir1),
            dir2: Some(axis_dir2),
            lock_translation: None,
            lock_rotation: None,
        }
    }

    /// Create a point-on-plane constraint: the plane is carried by the
    /// FIRST body (`plane_point` + unit `plane_normal`), the point is
    /// `point` on the SECOND body. One DOF removed.
    pub fn point_on_plane(
        id: usize,
        comp1: ComponentId,
        plane_point: Point3,
        plane_normal: Point3,
        comp2: ComponentId,
        point: Point3,
    ) -> Self {
        Self {
            id,
            constraint_type: ConstraintType::PointOnPlane,
            component1: comp1,
            point1: plane_point,
            component2: Some(comp2),
            point2: Some(point),
            distance: None,
            angle: None,
            fixed_orientation: None,
            dir1: Some(plane_normal),
            dir2: None,
            lock_translation: None,
            lock_rotation: None,
        }
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
        Self {
            id,
            constraint_type: ConstraintType::Lock,
            component1: comp1,
            point1: DVec3::ZERO,
            component2: Some(comp2),
            point2: None,
            distance: None,
            angle: None,
            fixed_orientation: None,
            dir1: None,
            dir2: None,
            lock_translation: Some(lock_translation),
            lock_rotation: Some(lock_rotation.normalize()),
        }
    }

    /// Get the number of DOF removed by this constraint.
    pub fn dof_removed(&self) -> usize {
        match self.constraint_type {
            ConstraintType::Fixed => {
                if self.fixed_orientation.is_some() {
                    6
                } else {
                    3
                }
            }
            ConstraintType::Coincident => 3,
            ConstraintType::Distance => 1,
            ConstraintType::Angle => 1,
            ConstraintType::Parallel => 2,
            ConstraintType::Perpendicular => 1,
            ConstraintType::Concentric => 4,
            ConstraintType::PointOnPlane => 1,
            ConstraintType::Lock => 6,
        }
    }
}

/// A system of constraints.
#[derive(Debug, Clone, Default)]
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
            .filter(|c| c.component1 == component || c.component2 == Some(component))
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
        assert_eq!(c.constraint_type, ConstraintType::Fixed);
        assert_eq!(c.dof_removed(), 3);
    }

    #[test]
    fn test_fixed_with_orientation_dof() {
        let c = Constraint::fixed_with_orientation(0, 0, DVec3::ZERO, DQuat::IDENTITY);
        assert_eq!(c.dof_removed(), 6);
        assert!(c.fixed_orientation.is_some());
    }

    #[test]
    fn test_coincident_constraint() {
        let c = Constraint::coincident(0, 1, DVec3::ZERO, 2, DVec3::new(1.0, 0.0, 0.0));
        assert_eq!(c.constraint_type, ConstraintType::Coincident);
        assert_eq!(c.dof_removed(), 3);
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
    }
}
