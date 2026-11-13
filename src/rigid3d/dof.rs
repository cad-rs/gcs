//! Degree of Freedom analysis for assemblies.

use crate::rigid3d::ConstraintSystem;

/// Degree of freedom type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegreeOfFreedom {
    /// Translation along X
    TranslateX,
    /// Translation along Y
    TranslateY,
    /// Translation along Z
    TranslateZ,
    /// Rotation about X
    RotateX,
    /// Rotation about Y
    RotateY,
    /// Rotation about Z
    RotateZ,
}

impl DegreeOfFreedom {
    /// Get all 6 DOFs.
    pub fn all() -> [Self; 6] {
        [
            Self::TranslateX,
            Self::TranslateY,
            Self::TranslateZ,
            Self::RotateX,
            Self::RotateY,
            Self::RotateZ,
        ]
    }
}

/// Result of DOF analysis.
#[derive(Debug, Clone)]
pub struct DOFResult {
    /// Total DOF in the system
    pub total_dof: usize,
    /// DOF per component
    pub component_dof: Vec<usize>,
    /// Is the system fully constrained
    pub is_fully_constrained: bool,
    /// Is the system over-constrained
    pub is_over_constrained: bool,
    /// Free DOFs by component
    pub free_dofs: Vec<Vec<DegreeOfFreedom>>,
}

/// DOF analyzer (count-based; see [`crate::rigid3d::diagnose_assembly`] for
/// the rank-based account).
#[derive(Debug, Clone)]
pub struct DOFAnalyzer {
    // Keeps the type constructible only through `new` / `default`.
    _private: (),
}

impl DOFAnalyzer {
    /// Create a new analyzer.
    pub fn new() -> Self {
        Self { _private: () }
    }

    /// Analyze DOF of a constraint system.
    pub fn analyze(&self, system: &ConstraintSystem) -> DOFResult {
        let num_components = system.num_components;
        let max_dof = num_components * 6;
        let total_removed = system.total_dof_removed();

        // Simple DOF calculation (ignores coupling)
        let remaining = max_dof.saturating_sub(total_removed);

        // Calculate DOF per component
        let mut component_dof = vec![6; num_components];
        let mut free_dofs: Vec<Vec<DegreeOfFreedom>> = (0..num_components)
            .map(|_| DegreeOfFreedom::all().to_vec())
            .collect();

        for constraint in &system.constraints {
            let dof_removed = constraint.dof_removed();

            // Remove DOFs from first component
            let comp1 = constraint.component1();
            for _ in 0..dof_removed.min(component_dof[comp1]) {
                if !free_dofs[comp1].is_empty() {
                    free_dofs[comp1].pop();
                }
            }
            component_dof[comp1] = free_dofs[comp1].len();

            // For two-component constraints, also affect second component
            if let Some(comp2) = constraint.component2() {
                // The second component loses DOFs as it's constrained relative to first
                let remaining_dof = constraint.constraint_type().dof_for_second();
                component_dof[comp2] = component_dof[comp2].min(remaining_dof);
                free_dofs[comp2].truncate(component_dof[comp2]);
            }
        }

        let is_fully_constrained = remaining == 0;
        let is_over_constrained = total_removed > max_dof;

        DOFResult {
            total_dof: remaining,
            component_dof,
            is_fully_constrained,
            is_over_constrained,
            free_dofs,
        }
    }
}

impl Default for DOFAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

/// Extension trait for ConstraintType to get DOF info.
trait ConstraintTypeExt {
    fn dof_for_second(&self) -> usize;
}

impl ConstraintTypeExt for crate::rigid3d::ConstraintType {
    fn dof_for_second(&self) -> usize {
        match self {
            Self::Fixed => 0,
            Self::Coincident => 3, // Second component still has 3 rotational DOFs
            Self::Distance => 5,   // Can rotate around the distance axis
            Self::Angle => 5,
            Self::Parallel => 4,
            Self::Perpendicular => 5,
            Self::Concentric => 2, // Second component: slide along + spin about the axis
            Self::PointOnPlane => 5, // One translational DOF removed
            Self::Lock => 0,       // Relative pose fully fixed
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rigid3d::Point3;
    use crate::rigid3d::{Constraint, ConstraintSystem};

    #[test]
    fn test_dof_all() {
        let all = DegreeOfFreedom::all();
        assert_eq!(all.len(), 6);
    }

    #[test]
    fn test_dof_unconstrained() {
        let system = ConstraintSystem::with_components(1);
        let analyzer = DOFAnalyzer::new();
        let result = analyzer.analyze(&system);

        assert_eq!(result.total_dof, 6);
        assert!(!result.is_fully_constrained);
    }

    #[test]
    fn test_dof_fixed_translation_only() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));

        let analyzer = DOFAnalyzer::new();
        let result = analyzer.analyze(&system);

        assert_eq!(result.total_dof, 3);
        assert!(!result.is_fully_constrained);
    }

    #[test]
    fn test_dof_fixed_full_pose() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed_with_orientation(
            0,
            0,
            Point3::ZERO,
            glam::DQuat::IDENTITY,
        ));

        let analyzer = DOFAnalyzer::new();
        let result = analyzer.analyze(&system);

        assert_eq!(result.total_dof, 0);
        assert!(result.is_fully_constrained);
    }

    #[test]
    fn test_dof_two_components() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        system.add_constraint(Constraint::coincident(
            1,
            0,
            Point3::ZERO,
            1,
            Point3::new(1.0, 0.0, 0.0),
        ));

        let analyzer = DOFAnalyzer::new();
        let result = analyzer.analyze(&system);

        // First component is fixed, second is constrained
        assert!(result.component_dof[0] < result.component_dof[1] || result.component_dof[1] == 0);
    }
}
