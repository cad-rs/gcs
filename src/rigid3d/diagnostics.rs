//! Assembly diagnostics as a list of issues — the former pattern-matching
//! analyzer, kept as a wrapper over the rank diagnosis.
//!
//! New code reads the structured result of
//! [`diagnose_assembly`][crate::rigid3d::diagnose_assembly] and
//! [`assembly_mobility`][crate::rigid3d::assembly_mobility] directly; the
//! messages here are a convenience on top of them.

use crate::rigid3d::diagnosis::{diagnose, mobility};
use crate::rigid3d::dof::solve_for_analysis;
use crate::rigid3d::{ConstraintKind, ConstraintSystem};

/// Type of diagnostic issue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiagnosticType {
    /// Constraints that depend on one another: some add nothing.
    Redundant,
    /// Constraints that cannot hold together.
    Conflict,
    /// Degrees of freedom remain.
    UnderConstrained,
    /// No longer reported: over-constraint shows as [`Self::Redundant`] or
    /// [`Self::Conflict`], which name the constraints involved.
    OverConstrained,
    /// No longer reported.
    Invalid,
    /// No longer reported: a closed loop of mates is an ordinary assembly,
    /// not a defect.
    CyclicDependency,
    /// Components that no constraint refers to.
    IncompleteChain,
}

/// Severity of a diagnostic issue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiagnosticSeverity {
    /// Error - prevents solving
    Error,
    /// Warning - may affect results
    Warning,
    /// Info - informational only
    Info,
}

impl DiagnosticSeverity {
    /// Get icon/label for this severity
    pub fn label(&self) -> &'static str {
        match self {
            Self::Error => "ERROR",
            Self::Warning => "WARNING",
            Self::Info => "INFO",
        }
    }
}

/// A recommendation for fixing a diagnostic issue.
#[derive(Debug, Clone)]
pub struct DiagnosticRecommendation {
    /// Description of the fix
    pub description: String,
    /// Constraint IDs to modify/remove (if applicable)
    pub affected_constraints: Vec<usize>,
    /// Priority (lower is more important)
    pub priority: u32,
}

/// A diagnostic issue found in the constraint system.
#[derive(Debug, Clone)]
pub struct DiagnosticIssue {
    /// Type of issue
    pub issue_type: DiagnosticType,
    /// Severity of the issue
    pub severity: DiagnosticSeverity,
    /// Constraint IDs involved
    pub constraint_ids: Vec<usize>,
    /// Component IDs involved
    pub component_ids: Vec<usize>,
    /// Human-readable message
    pub message: String,
    /// Recommendations for fixing the issue
    pub recommendations: Vec<DiagnosticRecommendation>,
}

/// Solver diagnostics.
///
/// This used to match patterns in the constraint list (two `Fixed` on one
/// body, more than three coincident pairs, any loop of mates) and count
/// freedoms by subtraction. It now solves the assembly and reports what the
/// rank diagnosis finds there, so it agrees with
/// [`diagnose_assembly`][crate::rigid3d::diagnose_assembly]; it is kept so
/// that existing callers keep working.
#[deprecated(
    since = "0.2.0",
    note = "use `diagnose_assembly` and `assembly_mobility`, which return the same findings as data"
)]
#[derive(Debug, Clone, Default)]
pub struct SolverDiagnostics {
    /// Issues found
    pub issues: Vec<DiagnosticIssue>,
    /// Summary statistics
    pub stats: DiagnosticStats,
}

/// Diagnostic statistics for the assembly.
#[derive(Debug, Clone, Default)]
pub struct DiagnosticStats {
    /// Number of components (removed ones are not counted)
    pub total_components: usize,
    /// Number of enforced constraints
    pub total_constraints: usize,
    /// Degrees of freedom the components have between them (six each)
    pub total_dof: usize,
    /// Degrees of freedom the constraints remove: the rank of their
    /// Jacobian, so overlapping constraints are not counted twice
    pub dof_removed: usize,
    /// Degrees of freedom remaining
    pub dof_remaining: usize,
    /// Number of components that carry a `Fixed` constraint
    pub fixed_components: usize,
    /// Number of components pinned in full (position and orientation)
    pub grounded_components: usize,
}

#[allow(deprecated)]
impl SolverDiagnostics {
    /// Create empty diagnostics.
    pub fn new() -> Self {
        Self::default()
    }

    /// Analyze a constraint system for issues.
    ///
    /// The assembly is solved from its seed pose (the one
    /// [`ConstraintSolver::solve`][crate::rigid3d::ConstraintSolver::solve]
    /// starts from) and diagnosed where that leaves it: at the solution,
    /// where redundancy and remaining freedoms are read, or — for a cluster
    /// that cannot be solved — back at the seed pose, where the conflict
    /// shows.
    pub fn analyze(&mut self, system: &ConstraintSystem) {
        self.issues.clear();

        let pose = solve_for_analysis(system);
        let diagnosis = diagnose(system, &pose.positions, &pose.rotations);
        let movable = mobility(system, &pose.positions, &pose.rotations);

        let mut fixed = vec![false; system.num_components];
        let mut grounded = vec![false; system.num_components];
        let mut referenced = vec![false; system.num_components];
        for (_, constraint) in system.active_constraints() {
            let (first, second) = constraint.kind.bodies();
            for body in std::iter::once(first).chain(second) {
                referenced[body] = true;
            }
            if let ConstraintKind::Fixed {
                body, orientation, ..
            } = constraint.kind
            {
                fixed[body] = true;
                grounded[body] |= orientation.is_some();
            }
        }
        let count = |flags: &[bool]| flags.iter().filter(|flag| **flag).count();

        let enforced = system.active_constraints().count();
        self.stats = DiagnosticStats {
            total_components: system.live_components(),
            total_constraints: enforced,
            total_dof: system.max_dof(),
            dof_removed: diagnosis.rank,
            dof_remaining: diagnosis.dof.max(0) as usize,
            fixed_components: count(&fixed),
            grounded_components: count(&grounded),
        };

        // The components the constraints with these ids act on.
        let components_of = |ids: &[usize]| -> Vec<usize> {
            let mut out: Vec<usize> = system
                .active_constraints()
                .map(|(_, c)| c)
                .filter(|c| ids.contains(&c.id))
                .flat_map(|c| {
                    let (first, second) = c.kind.bodies();
                    std::iter::once(first).chain(second)
                })
                .collect();
            out.sort_unstable();
            out.dedup();
            out
        };

        // Conflicts.
        if !diagnosis.conflicting.is_empty() {
            let ids = diagnosis.conflicting.clone();
            self.issues.push(DiagnosticIssue {
                issue_type: DiagnosticType::Conflict,
                severity: DiagnosticSeverity::Error,
                component_ids: components_of(&ids),
                message: format!("Constraints {ids:?} cannot hold together"),
                recommendations: vec![DiagnosticRecommendation {
                    description: "Remove or change one of the conflicting constraints".to_string(),
                    affected_constraints: ids.clone(),
                    priority: 1,
                }],
                constraint_ids: ids,
            });
        } else if !pose.success {
            self.issues.push(DiagnosticIssue {
                issue_type: DiagnosticType::Conflict,
                severity: DiagnosticSeverity::Error,
                constraint_ids: vec![],
                component_ids: vec![],
                message: format!(
                    "The assembly could not be solved (error {:.3e}); no single set of \
                     conflicting constraints was identified",
                    pose.error
                ),
                recommendations: vec![DiagnosticRecommendation {
                    description: "Check the most recently added constraints".to_string(),
                    affected_constraints: vec![],
                    priority: 1,
                }],
            });
        }

        // Redundancy: one issue per group of mutually dependent constraints,
        // or a single one when no group could be extracted.
        let groups: Vec<Vec<usize>> = if diagnosis.dependency_groups.is_empty() {
            let mut all = diagnosis.redundant.clone();
            all.extend(&diagnosis.partially_redundant);
            if all.is_empty() {
                Vec::new()
            } else {
                vec![all]
            }
        } else {
            diagnosis.dependency_groups.clone()
        };
        for group in groups {
            // The members the diagnosis singled out as removable.
            let mut removable: Vec<usize> = group
                .iter()
                .copied()
                .filter(|id| diagnosis.redundant.contains(id))
                .collect();
            if removable.is_empty() {
                removable = group.clone();
            }
            self.issues.push(DiagnosticIssue {
                issue_type: DiagnosticType::Redundant,
                severity: DiagnosticSeverity::Warning,
                component_ids: components_of(&group),
                message: format!("Constraints {group:?} depend on one another"),
                recommendations: vec![DiagnosticRecommendation {
                    description: "Remove one of the constraints that add nothing".to_string(),
                    affected_constraints: removable,
                    priority: 1,
                }],
                constraint_ids: group,
            });
        }

        // Remaining freedoms. An assembly without any constraint is not
        // reported: there is nothing to complete yet.
        if diagnosis.dof > 0 && enforced > 0 {
            let moving: Vec<usize> = (0..system.num_components)
                .filter(|&body| movable[body].iter().any(|can_move| *can_move))
                .collect();
            let mut recommendations = vec![DiagnosticRecommendation {
                description: "Add constraints to the components that can still move".to_string(),
                affected_constraints: vec![],
                priority: 1,
            }];
            if count(&grounded) == 0 {
                recommendations.push(DiagnosticRecommendation {
                    description: "Ground one component with a fixed constraint as reference"
                        .to_string(),
                    affected_constraints: vec![],
                    priority: 2,
                });
            }
            self.issues.push(DiagnosticIssue {
                issue_type: DiagnosticType::UnderConstrained,
                severity: DiagnosticSeverity::Warning,
                constraint_ids: vec![],
                component_ids: moving,
                message: format!(
                    "System is under-constrained: {} DOF remaining",
                    diagnosis.dof
                ),
                recommendations,
            });
        }

        // Components nothing refers to.
        if enforced > 0 {
            let isolated: Vec<usize> = (0..system.num_components)
                .filter(|&body| !referenced[body] && !system.is_component_removed(body))
                .collect();
            if !isolated.is_empty() {
                self.issues.push(DiagnosticIssue {
                    issue_type: DiagnosticType::IncompleteChain,
                    severity: DiagnosticSeverity::Warning,
                    constraint_ids: vec![],
                    message: format!(
                        "{} component(s) are not connected to the assembly: {:?}",
                        isolated.len(),
                        isolated
                    ),
                    component_ids: isolated,
                    recommendations: vec![DiagnosticRecommendation {
                        description: "Add mates to connect isolated components to the assembly"
                            .to_string(),
                        affected_constraints: vec![],
                        priority: 1,
                    }],
                });
            }
        }
    }

    /// Check if system has any errors.
    pub fn has_errors(&self) -> bool {
        self.issues
            .iter()
            .any(|i| i.severity == DiagnosticSeverity::Error)
    }

    /// Check if system has any warnings.
    pub fn has_warnings(&self) -> bool {
        self.issues
            .iter()
            .any(|i| i.severity == DiagnosticSeverity::Warning)
    }

    /// Check if system has any info messages.
    pub fn has_info(&self) -> bool {
        self.issues
            .iter()
            .any(|i| i.severity == DiagnosticSeverity::Info)
    }

    /// Get all error issues.
    pub fn errors(&self) -> Vec<&DiagnosticIssue> {
        self.issues
            .iter()
            .filter(|i| i.severity == DiagnosticSeverity::Error)
            .collect()
    }

    /// Get all warning issues.
    pub fn warnings(&self) -> Vec<&DiagnosticIssue> {
        self.issues
            .iter()
            .filter(|i| i.severity == DiagnosticSeverity::Warning)
            .collect()
    }

    /// Get issues by type.
    pub fn issues_by_type(&self, issue_type: DiagnosticType) -> Vec<&DiagnosticIssue> {
        self.issues
            .iter()
            .filter(|i| i.issue_type == issue_type)
            .collect()
    }

    /// Get a summary report as a string.
    pub fn summary(&self) -> String {
        let mut report = String::new();

        report.push_str(&format!(
            "Diagnostic Summary:\n  Components: {}\n  Constraints: {}\n  DOF: {} remaining, {} removed\n",
            self.stats.total_components,
            self.stats.total_constraints,
            self.stats.dof_remaining,
            self.stats.dof_removed
        ));

        if self.issues.is_empty() {
            report.push_str("  No issues found.\n");
        } else {
            report.push_str(&format!("  Issues: {}\n", self.issues.len()));
            for issue in &self.issues {
                report.push_str(&format!(
                    "  [{}] {}: {}\n",
                    issue.severity.label(),
                    issue.issue_type_name(),
                    issue.message
                ));
            }
        }

        report
    }
}

impl DiagnosticIssue {
    /// Get a human-readable name for this issue type.
    pub fn issue_type_name(&self) -> &'static str {
        match self.issue_type {
            DiagnosticType::Redundant => "Redundant",
            DiagnosticType::Conflict => "Conflict",
            DiagnosticType::UnderConstrained => "Under-Constrained",
            DiagnosticType::OverConstrained => "Over-Constrained",
            DiagnosticType::Invalid => "Invalid",
            DiagnosticType::CyclicDependency => "Cyclic Dependency",
            DiagnosticType::IncompleteChain => "Incomplete Chain",
        }
    }
}

#[cfg(test)]
#[allow(deprecated)]
mod tests {
    use super::*;
    use crate::rigid3d::Point3;
    use crate::rigid3d::{Constraint, ConstraintSystem};
    use glam::DQuat;

    fn analyze(system: &ConstraintSystem) -> SolverDiagnostics {
        let mut diag = SolverDiagnostics::new();
        diag.analyze(system);
        diag
    }

    fn ground(system: &mut ConstraintSystem, id: usize, body: usize) {
        system.add_constraint(Constraint::fixed_with_orientation(
            id,
            body,
            Point3::ZERO,
            DQuat::IDENTITY,
        ));
    }

    #[test]
    fn two_pins_that_disagree_are_a_conflict() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        system.add_constraint(Constraint::fixed(1, 0, Point3::new(1.0, 0.0, 0.0)));
        let diag = analyze(&system);

        assert!(diag.has_errors());
        let conflicts = diag.issues_by_type(DiagnosticType::Conflict);
        assert_eq!(conflicts.len(), 1);
        assert_eq!(conflicts[0].severity, DiagnosticSeverity::Error);
        assert!(!conflicts[0].constraint_ids.is_empty());
        assert_eq!(conflicts[0].component_ids, vec![0]);
        assert_eq!(
            conflicts[0].recommendations[0].affected_constraints,
            conflicts[0].constraint_ids
        );
        assert_eq!(diag.errors().len(), 1);
    }

    #[test]
    fn two_pins_that_agree_are_redundant_not_conflicting() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, Point3::new(1.0, 2.0, 3.0)));
        system.add_constraint(Constraint::fixed(1, 0, Point3::new(1.0, 2.0, 3.0)));
        let diag = analyze(&system);

        assert!(!diag.has_errors());
        let redundant = diag.issues_by_type(DiagnosticType::Redundant);
        assert_eq!(redundant.len(), 1, "{:?}", diag.issues);
        assert_eq!(redundant[0].constraint_ids, vec![0, 1]);
        assert_eq!(redundant[0].severity, DiagnosticSeverity::Warning);
        // Only three freedoms are gone, however many times they are pinned.
        assert_eq!(diag.stats.dof_removed, 3);
        assert_eq!(diag.stats.dof_remaining, 3);
    }

    #[test]
    fn an_assembly_without_constraints_reports_nothing() {
        let diag = analyze(&ConstraintSystem::with_components(2));
        assert!(diag.issues.is_empty());
        assert_eq!(diag.stats.dof_remaining, 12);
        assert!(diag.summary().contains("No issues found"));
    }

    #[test]
    fn remaining_freedoms_are_a_warning_naming_what_can_move() {
        let mut system = ConstraintSystem::with_components(2);
        ground(&mut system, 0, 0);
        system.add_constraint(Constraint::coincident(
            1,
            0,
            Point3::new(1.0, 0.0, 0.0),
            1,
            Point3::ZERO,
        ));
        let diag = analyze(&system);

        assert!(!diag.has_errors());
        let under = diag.issues_by_type(DiagnosticType::UnderConstrained);
        assert_eq!(under.len(), 1);
        assert_eq!(under[0].component_ids, vec![1], "the ground cannot move");
        assert!(under[0].message.contains("3 DOF"));
        // A ground exists, so grounding one is not suggested.
        assert_eq!(under[0].recommendations.len(), 1);
        assert!(diag.has_warnings());
        assert!(!diag.has_info());
    }

    #[test]
    fn a_fully_determined_assembly_has_no_issues() {
        let mut system = ConstraintSystem::with_components(2);
        ground(&mut system, 0, 0);
        system.add_constraint(Constraint::lock(
            1,
            0,
            1,
            Point3::new(1.0, 2.0, 3.0),
            DQuat::from_rotation_z(0.5),
        ));
        let diag = analyze(&system);
        assert!(diag.issues.is_empty(), "{:?}", diag.issues);
        assert_eq!(diag.stats.dof_remaining, 0);
        assert_eq!(diag.stats.dof_removed, 12);
    }

    /// A closed loop of mates is an ordinary assembly: the old analyzer
    /// called it a cyclic dependency.
    #[test]
    fn a_closed_loop_is_not_an_issue() {
        let mut system = ConstraintSystem::with_components(3);
        ground(&mut system, 0, 0);
        let offset = Point3::new(1.0, 0.0, 0.0);
        system.add_constraint(Constraint::lock(1, 0, 1, offset, DQuat::IDENTITY));
        system.add_constraint(Constraint::lock(2, 1, 2, offset, DQuat::IDENTITY));
        // Closing the loop with a distance that already holds.
        system.add_constraint(Constraint::distance(
            3,
            2,
            Point3::ZERO,
            0,
            Point3::ZERO,
            2.0,
        ));
        let diag = analyze(&system);
        assert!(diag
            .issues_by_type(DiagnosticType::CyclicDependency)
            .is_empty());
        assert!(!diag.has_errors());
        // The closing distance adds nothing, and that is what is reported.
        let redundant = diag.issues_by_type(DiagnosticType::Redundant);
        assert_eq!(redundant.len(), 1, "{:?}", diag.issues);
        assert!(redundant[0].constraint_ids.contains(&3));
        assert_eq!(
            redundant[0].recommendations[0].affected_constraints,
            vec![3]
        );
    }

    #[test]
    fn test_diagnostic_stats() {
        let mut system = ConstraintSystem::with_components(3);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        system.add_constraint(Constraint::coincident(
            1,
            0,
            Point3::ZERO,
            1,
            Point3::new(1.0, 0.0, 0.0),
        ));
        let diag = analyze(&system);

        assert_eq!(diag.stats.total_components, 3);
        assert_eq!(diag.stats.total_constraints, 2);
        assert_eq!(diag.stats.total_dof, 18);
        assert_eq!(diag.stats.dof_removed, 6);
        assert_eq!(diag.stats.dof_remaining, 12);
        assert_eq!(diag.stats.fixed_components, 1);
        assert_eq!(diag.stats.grounded_components, 0);
    }

    #[test]
    fn components_nothing_refers_to_are_reported() {
        let mut system = ConstraintSystem::with_components(3);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        let diag = analyze(&system);

        let isolated = diag.issues_by_type(DiagnosticType::IncompleteChain);
        assert_eq!(isolated.len(), 1);
        assert_eq!(isolated[0].component_ids, vec![1, 2]);
    }

    #[test]
    fn test_diagnostic_summary() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        let summary = analyze(&system).summary();
        assert!(summary.contains("Components: 2"));
        assert!(summary.contains("Constraints: 1"));
        assert!(summary.contains("9 remaining, 3 removed"));
        assert!(summary.contains("[WARNING] Under-Constrained"));
    }
}
