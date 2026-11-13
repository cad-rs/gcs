//! Degree of freedom analysis for assemblies — the former count-based
//! analyzer, kept as a wrapper over the rank-based account.
//!
//! New code reads the degrees of freedom from
//! [`diagnose_assembly`][crate::rigid3d::diagnose_assembly] and the
//! coordinates that can still move from
//! [`assembly_mobility`][crate::rigid3d::assembly_mobility].

use crate::rigid3d::diagnosis::{diagnose, mobility};
use crate::rigid3d::{ConstraintSolver, ConstraintSystem, PoseSolveOptions, PoseSolveResult};

/// Degree of freedom type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegreeOfFreedom {
    /// Translation along world X
    TranslateX,
    /// Translation along world Y
    TranslateY,
    /// Translation along world Z
    TranslateZ,
    /// Rotation about the body's own x axis
    RotateX,
    /// Rotation about the body's own y axis
    RotateY,
    /// Rotation about the body's own z axis
    RotateZ,
}

impl DegreeOfFreedom {
    /// Get all 6 DOFs, in the order of a body's six coordinates.
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
    /// Degrees of freedom of the assembly: six per component minus the
    /// rank of the constraints (the `dof` of the rank diagnosis).
    pub total_dof: usize,
    /// Per component, how many of its six coordinates can still move (the
    /// length of its entry in `free_dofs`). This counts coordinates, not
    /// independent motions: a body left free to slide along a diagonal has
    /// two movable translations and one freedom.
    pub component_dof: Vec<usize>,
    /// No freedom is left.
    pub is_fully_constrained: bool,
    /// Some constraint is redundant or conflicts with the others.
    pub is_over_constrained: bool,
    /// Per component, the coordinates that can still move.
    pub free_dofs: Vec<Vec<DegreeOfFreedom>>,
}

/// DOF analyzer.
///
/// This used to subtract a fixed number of freedoms per constraint, which
/// is wrong as soon as two constraints overlap. It now solves the assembly
/// and reads the answer off the rank diagnosis, so it agrees with
/// [`diagnose_assembly`][crate::rigid3d::diagnose_assembly]; it is kept so
/// that existing callers keep working.
#[deprecated(
    since = "0.2.0",
    note = "use `diagnose_assembly` for the degrees of freedom and `assembly_mobility` for the coordinates that can move"
)]
#[derive(Debug, Clone)]
pub struct DOFAnalyzer {
    // Keeps the type constructible only through `new` / `default`.
    _private: (),
}

/// Solve an assembly tightly from its seed pose, for analysis. Rank is
/// read at a solution when there is one: redundancy is only meaningful
/// there. A cluster that cannot be solved is analyzed at the seed pose,
/// which is where a failed solve leaves it.
pub(crate) fn solve_for_analysis(system: &ConstraintSystem) -> PoseSolveResult {
    ConstraintSolver::with_options(PoseSolveOptions {
        tolerance: 1.0e-10,
        ..PoseSolveOptions::default()
    })
    .solve(system)
}

#[allow(deprecated)]
impl DOFAnalyzer {
    /// Create a new analyzer.
    pub fn new() -> Self {
        Self { _private: () }
    }

    /// Analyze DOF of a constraint system.
    ///
    /// The assembly is solved from its seed pose (the one
    /// [`ConstraintSolver::solve`] starts from) and analyzed where that
    /// leaves it.
    pub fn analyze(&self, system: &ConstraintSystem) -> DOFResult {
        let pose = solve_for_analysis(system);
        let diagnosis = diagnose(system, &pose.positions, &pose.rotations);
        let free_dofs: Vec<Vec<DegreeOfFreedom>> =
            mobility(system, &pose.positions, &pose.rotations)
                .iter()
                .map(|movable| {
                    DegreeOfFreedom::all()
                        .into_iter()
                        .zip(movable)
                        .filter(|(_, &can_move)| can_move)
                        .map(|(dof, _)| dof)
                        .collect()
                })
                .collect();

        DOFResult {
            total_dof: diagnosis.dof.max(0) as usize,
            component_dof: free_dofs.iter().map(Vec::len).collect(),
            is_fully_constrained: diagnosis.dof == 0,
            is_over_constrained: !diagnosis.redundant.is_empty()
                || !diagnosis.conflicting.is_empty(),
            free_dofs,
        }
    }
}

#[allow(deprecated)]
impl Default for DOFAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
#[allow(deprecated)]
mod tests {
    use super::*;
    use crate::rigid3d::Point3;
    use crate::rigid3d::{Constraint, ConstraintSystem};
    use glam::DQuat;

    #[test]
    fn test_dof_all() {
        let all = DegreeOfFreedom::all();
        assert_eq!(all.len(), 6);
    }

    #[test]
    fn test_dof_unconstrained() {
        let system = ConstraintSystem::with_components(1);
        let result = DOFAnalyzer::new().analyze(&system);

        assert_eq!(result.total_dof, 6);
        assert!(!result.is_fully_constrained);
        assert!(!result.is_over_constrained);
        assert_eq!(result.component_dof, vec![6]);
        assert_eq!(result.free_dofs[0], DegreeOfFreedom::all().to_vec());
    }

    #[test]
    fn test_dof_fixed_translation_only() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        let result = DOFAnalyzer::new().analyze(&system);

        assert_eq!(result.total_dof, 3);
        assert!(!result.is_fully_constrained);
        assert_eq!(
            result.free_dofs[0],
            vec![
                DegreeOfFreedom::RotateX,
                DegreeOfFreedom::RotateY,
                DegreeOfFreedom::RotateZ
            ]
        );
    }

    #[test]
    fn test_dof_fixed_full_pose() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed_with_orientation(
            0,
            0,
            Point3::ZERO,
            DQuat::IDENTITY,
        ));
        let result = DOFAnalyzer::new().analyze(&system);

        assert_eq!(result.total_dof, 0);
        assert!(result.is_fully_constrained);
        assert_eq!(result.component_dof, vec![0]);
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
        let result = DOFAnalyzer::new().analyze(&system);

        // The first body can only turn about its pinned origin. The second
        // hangs on it by one point: it can turn about that point, and that
        // moves its origin.
        assert_eq!(result.total_dof, 12 - 3 - 3);
        assert_eq!(result.component_dof[0], 3);
        assert!(result.component_dof[1] > result.component_dof[0]);
    }

    /// The count the old analyzer got wrong: overlapping constraints. Three
    /// coincident pairs pin a body with nine equations; six freedoms are
    /// gone, not nine, and the surplus shows as over-constraint.
    #[test]
    fn overlapping_constraints_are_not_counted_twice() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed_with_orientation(
            0,
            0,
            Point3::ZERO,
            DQuat::IDENTITY,
        ));
        for (id, local) in [
            Point3::ZERO,
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 3.0, 1.0),
        ]
        .into_iter()
        .enumerate()
        {
            let world = local + Point3::new(1.0, 1.0, 1.0);
            system.add_constraint(Constraint::coincident(id + 1, 0, world, 1, local));
        }
        assert_eq!(system.total_dof_removed(), 6 + 9, "what counting says");

        let result = DOFAnalyzer::new().analyze(&system);
        assert_eq!(result.total_dof, 0);
        assert!(result.is_fully_constrained);
        assert!(result.is_over_constrained);
        assert_eq!(result.component_dof, vec![0, 0]);
    }
}
