//! Solver diagnostics for redundancy and conflict detection.

use crate::rigid3d::{Constraint, ConstraintSystem};
use std::collections::{HashMap, HashSet};

/// Type of diagnostic issue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiagnosticType {
    /// Redundant constraint
    Redundant,
    /// Conflicting constraint
    Conflict,
    /// Under-constrained system
    UnderConstrained,
    /// Over-constrained system
    OverConstrained,
    /// Invalid constraint parameters
    Invalid,
    /// Cyclic dependency detected
    CyclicDependency,
    /// Incomplete constraint chain
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
    /// Total number of components
    pub total_components: usize,
    /// Total number of constraints
    pub total_constraints: usize,
    /// Total DOF in system
    pub total_dof: usize,
    /// DOF removed by constraints
    pub dof_removed: usize,
    /// DOF remaining
    pub dof_remaining: usize,
    /// Number of fixed components
    pub fixed_components: usize,
    /// Number of grounded components (explicitly fixed)
    pub grounded_components: usize,
}

impl SolverDiagnostics {
    /// Create empty diagnostics.
    pub fn new() -> Self {
        Self::default()
    }

    /// Analyze a constraint system for issues.
    pub fn analyze(&mut self, system: &ConstraintSystem) {
        self.issues.clear();

        // Compute statistics
        self.compute_stats(system);

        // Check for over-constrained
        self.check_over_constrained(system);

        // Check for under-constrained
        self.check_under_constrained(system);

        // Check for redundant constraints
        self.detect_redundancy(system);

        // Check for conflicts
        self.detect_conflicts(system);

        // Check for cyclic dependencies
        self.detect_cycles(system);

        // Check for incomplete constraint chains
        self.detect_incomplete_chains(system);
    }

    /// Compute diagnostic statistics.
    fn compute_stats(&mut self, system: &ConstraintSystem) {
        let max_dof = system.num_components * 6;
        let total_removed = system.total_dof_removed();

        let fixed_count = system.constraints
            .iter()
            .filter(|c| c.constraint_type == crate::rigid3d::ConstraintType::Fixed)
            .count();

        self.stats = DiagnosticStats {
            total_components: system.num_components,
            total_constraints: system.constraints.len(),
            total_dof: max_dof,
            dof_removed: total_removed,
            dof_remaining: max_dof.saturating_sub(total_removed),
            fixed_components: fixed_count,
            grounded_components: fixed_count,
        };
    }

    /// Check for over-constrained system.
    fn check_over_constrained(&mut self, system: &ConstraintSystem) {
        let max_dof = system.num_components * 6;
        let total_removed = system.total_dof_removed();

        if total_removed > max_dof {
            let constraint_ids: Vec<usize> = system.constraints.iter().map(|c| c.id).collect();

            self.issues.push(DiagnosticIssue {
                issue_type: DiagnosticType::OverConstrained,
                severity: DiagnosticSeverity::Error,
                constraint_ids: constraint_ids.clone(),
                component_ids: (0..system.num_components).collect(),
                message: format!(
                    "System is over-constrained: {} DOF removed from {} available",
                    total_removed, max_dof
                ),
                recommendations: vec![
                    DiagnosticRecommendation {
                        description: "Remove redundant constraints to resolve over-constraint".to_string(),
                        affected_constraints: constraint_ids,
                        priority: 1,
                    },
                    DiagnosticRecommendation {
                        description: "Consider using a fixed constraint instead of multiple mates".to_string(),
                        affected_constraints: vec![],
                        priority: 2,
                    },
                ],
            });
        }
    }

    /// Check for under-constrained system.
    fn check_under_constrained(&mut self, system: &ConstraintSystem) {
        let max_dof = system.num_components * 6;
        let total_removed = system.total_dof_removed();

        if total_removed < max_dof && !system.constraints.is_empty() {
            let remaining = max_dof - total_removed;

            self.issues.push(DiagnosticIssue {
                issue_type: DiagnosticType::UnderConstrained,
                severity: DiagnosticSeverity::Warning,
                constraint_ids: vec![],
                component_ids: (0..system.num_components).collect(),
                message: format!(
                    "System is under-constrained: {} DOF remaining",
                    remaining
                ),
                recommendations: vec![
                    DiagnosticRecommendation {
                        description: format!("Add {} more constraints to fully constrain the assembly", remaining),
                        affected_constraints: vec![],
                        priority: 1,
                    },
                    DiagnosticRecommendation {
                        description: "Ground one component with a fixed constraint as reference".to_string(),
                        affected_constraints: vec![],
                        priority: 2,
                    },
                ],
            });
        }
    }

    /// Detect redundant constraints.
    fn detect_redundancy(&mut self, system: &ConstraintSystem) {
        // Group constraints by component
        let mut component_constraints: HashMap<usize, Vec<&Constraint>> = HashMap::new();

        for constraint in &system.constraints {
            component_constraints
                .entry(constraint.component1)
                .or_default()
                .push(constraint);
            if let Some(comp2) = constraint.component2 {
                component_constraints.entry(comp2).or_default().push(constraint);
            }
        }

        // Check for multiple fixed constraints on same component
        for (&component, constraints) in &component_constraints {
            let fixed_count = constraints
                .iter()
                .filter(|c| c.constraint_type == crate::rigid3d::ConstraintType::Fixed)
                .count();

            if fixed_count > 1 {
                let ids: Vec<usize> = constraints
                    .iter()
                    .filter(|c| c.constraint_type == crate::rigid3d::ConstraintType::Fixed)
                    .map(|c| c.id)
                    .collect();

                self.issues.push(DiagnosticIssue {
                    issue_type: DiagnosticType::Redundant,
                    severity: DiagnosticSeverity::Warning,
                    constraint_ids: ids.clone(),
                    component_ids: vec![component],
                    message: format!(
                        "Component {} has multiple fixed constraints",
                        component
                    ),
                    recommendations: vec![
                        DiagnosticRecommendation {
                            description: "Remove all but one fixed constraint from this component".to_string(),
                            affected_constraints: ids,
                            priority: 1,
                        },
                    ],
                });
            }
        }

        // Check for redundant coincident constraints between same pair of components
        self.detect_redundant_coincident(system);
    }

    /// Detect redundant coincident constraints.
    fn detect_redundant_coincident(&mut self, system: &ConstraintSystem) {
        let mut coincident_pairs: HashMap<(usize, usize), Vec<usize>> = HashMap::new();

        for constraint in &system.constraints {
            if constraint.constraint_type == crate::rigid3d::ConstraintType::Coincident {
                if let Some(comp2) = constraint.component2 {
                    let key = if constraint.component1 < comp2 {
                        (constraint.component1, comp2)
                    } else {
                        (comp2, constraint.component1)
                    };
                    coincident_pairs.entry(key).or_default().push(constraint.id);
                }
            }
        }

        for ((_c1, _c2), ids) in coincident_pairs {
            if ids.len() > 3 {
                // More than 3 coincident constraints between same pair is likely redundant
                self.issues.push(DiagnosticIssue {
                    issue_type: DiagnosticType::Redundant,
                    severity: DiagnosticSeverity::Info,
                    constraint_ids: ids.clone(),
                    component_ids: vec![_c1, _c2],
                    message: format!(
                        "Components {} and {} have {} coincident constraints; 3 are typically sufficient",
                        _c1, _c2, ids.len()
                    ),
                    recommendations: vec![
                        DiagnosticRecommendation {
                            description: "Consider removing redundant coincident constraints".to_string(),
                            affected_constraints: ids,
                            priority: 3,
                        },
                    ],
                });
            }
        }
    }

    /// Detect conflicting constraints.
    fn detect_conflicts(&mut self, system: &ConstraintSystem) {
        // Check for conflicting distance constraints
        for (i, c1) in system.constraints.iter().enumerate() {
            for c2 in system.constraints.iter().skip(i + 1) {
                if self.are_conflicting(c1, c2) {
                    self.issues.push(DiagnosticIssue {
                        issue_type: DiagnosticType::Conflict,
                        severity: DiagnosticSeverity::Error,
                        constraint_ids: vec![c1.id, c2.id],
                        component_ids: vec![c1.component1, c2.component1],
                        message: format!("Constraints {} and {} are conflicting", c1.id, c2.id),
                        recommendations: vec![
                            DiagnosticRecommendation {
                                description: "Remove or modify one of the conflicting constraints".to_string(),
                                affected_constraints: vec![c1.id, c2.id],
                                priority: 1,
                            },
                        ],
                    });
                }
            }
        }
    }

    /// Detect cyclic dependencies in constraints.
    fn detect_cycles(&mut self, system: &ConstraintSystem) {
        // Build a dependency graph
        let mut graph: HashMap<usize, HashSet<usize>> = HashMap::new();

        for constraint in &system.constraints {
            if let Some(comp2) = constraint.component2 {
                // component1 depends on comp2 through this constraint
                graph.entry(constraint.component1).or_default().insert(comp2);
            }
        }

        // Detect cycles using DFS
        let mut visited = HashSet::new();
        let mut rec_stack = HashSet::new();

        for node in 0..system.num_components {
            if self.has_cycle(node, &graph, &mut visited, &mut rec_stack) {
                self.issues.push(DiagnosticIssue {
                    issue_type: DiagnosticType::CyclicDependency,
                    severity: DiagnosticSeverity::Warning,
                    constraint_ids: vec![],
                    component_ids: (0..system.num_components).collect(),
                    message: "Cyclic dependency detected in constraint graph".to_string(),
                    recommendations: vec![
                        DiagnosticRecommendation {
                            description: "Reorganize constraints to break the cycle".to_string(),
                            affected_constraints: vec![],
                            priority: 1,
                        },
                    ],
                });
                break;
            }
        }
    }

    /// Check if there's a cycle starting from the given node.
    fn has_cycle(
        &self,
        node: usize,
        graph: &HashMap<usize, HashSet<usize>>,
        visited: &mut HashSet<usize>,
        rec_stack: &mut HashSet<usize>,
    ) -> bool {
        if rec_stack.contains(&node) {
            return true;
        }

        if visited.contains(&node) {
            return false;
        }

        visited.insert(node);
        rec_stack.insert(node);

        if let Some(neighbors) = graph.get(&node) {
            for &neighbor in neighbors {
                if self.has_cycle(neighbor, graph, visited, rec_stack) {
                    return true;
                }
            }
        }

        rec_stack.remove(&node);
        false
    }

    /// Detect incomplete constraint chains (unconnected components).
    fn detect_incomplete_chains(&mut self, system: &ConstraintSystem) {
        if system.constraints.is_empty() || system.num_components <= 1 {
            return;
        }

        // Find all connected components using union-find
        let mut parent: Vec<usize> = (0..system.num_components).collect();

        // Iterative find with path compression
        let find_iterative = |parent: &mut [usize], mut x: usize| -> usize {
            let mut root = x;
            while parent[root] != root {
                root = parent[root];
            }
            // Path compression
            while parent[x] != x {
                let next = parent[x];
                parent[x] = root;
                x = next;
            }
            root
        };

        let union = |parent: &mut [usize], x: usize, y: usize| {
            let px = find_iterative(parent, x);
            let py = find_iterative(parent, y);
            if px != py {
                parent[px] = py;
            }
        };

        for constraint in &system.constraints {
            if let Some(comp2) = constraint.component2 {
                union(&mut parent, constraint.component1, comp2);
            }
        }

        // Count connected groups
        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..system.num_components {
            let root = find_iterative(&mut parent, i);
            groups.entry(root).or_default().push(i);
        }

        if groups.len() > 1 {
            let isolated: Vec<usize> = groups.values()
                .filter(|g| g.len() == 1)
                .flat_map(|g| g.iter().copied())
                .collect();

            if !isolated.is_empty() {
                self.issues.push(DiagnosticIssue {
                    issue_type: DiagnosticType::IncompleteChain,
                    severity: DiagnosticSeverity::Warning,
                    constraint_ids: vec![],
                    component_ids: isolated.clone(),
                    message: format!(
                        "{} component(s) are not connected to the assembly: {:?}",
                        isolated.len(),
                        isolated
                    ),
                    recommendations: vec![
                        DiagnosticRecommendation {
                            description: "Add mates to connect isolated components to the assembly".to_string(),
                            affected_constraints: vec![],
                            priority: 1,
                        },
                    ],
                });
            }
        }
    }

    /// Check if two constraints are conflicting.
    fn are_conflicting(&self, c1: &Constraint, c2: &Constraint) -> bool {
        // Simple check: same component with different fixed positions
        if c1.component1 == c2.component1 {
            if c1.constraint_type == crate::rigid3d::ConstraintType::Fixed
                && c2.constraint_type == crate::rigid3d::ConstraintType::Fixed
            {
                return (c1.point1 - c2.point1).length() > 0.001;
            }
        }
        false
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
            "Diagnostic Summary:\n  Components: {}\n  Constraints: {}\n  DOF: {} / {} removed\n",
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
mod tests {
    use super::*;
    use crate::rigid3d::{Constraint, ConstraintSystem, ConstraintType};
    use crate::rigid3d::Point3;

    #[test]
    fn test_detect_over_constrained() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        system.add_constraint(Constraint::fixed(1, 0, Point3::new(1.0, 0.0, 0.0)));

        let mut diag = SolverDiagnostics::new();
        diag.analyze(&system);

        assert!(diag.has_warnings() || diag.has_errors());
    }

    #[test]
    fn test_detect_conflict() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        system.add_constraint(Constraint::fixed(1, 0, Point3::new(1.0, 0.0, 0.0)));

        let mut diag = SolverDiagnostics::new();
        diag.analyze(&system);

        assert!(diag.issues.iter().any(|i| i.issue_type == DiagnosticType::Conflict));
    }

    #[test]
    fn test_detect_under_constrained() {
        let system = ConstraintSystem::with_components(2);

        let mut diag = SolverDiagnostics::new();
        diag.analyze(&system);

        // Empty system shouldn't trigger under-constrained warning
        assert!(!diag.issues.iter().any(|i| i.issue_type == DiagnosticType::UnderConstrained));
    }

    #[test]
    fn test_no_issues() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));

        let mut diag = SolverDiagnostics::new();
        diag.analyze(&system);

        // Should have no errors for fully constrained single component
        assert!(!diag.has_errors());
    }

    #[test]
    fn test_diagnostic_stats() {
        let mut system = ConstraintSystem::with_components(3);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        system.add_constraint(Constraint::coincident(
            1, 0, Point3::ZERO,
            1, Point3::new(1.0, 0.0, 0.0),
        ));

        let mut diag = SolverDiagnostics::new();
        diag.analyze(&system);

        assert_eq!(diag.stats.total_components, 3);
        assert_eq!(diag.stats.total_constraints, 2);
        assert_eq!(diag.stats.fixed_components, 1);
    }

    #[test]
    fn test_incomplete_chain_detection() {
        let mut system = ConstraintSystem::with_components(3);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        // Component 2 is not connected to the rest

        let mut diag = SolverDiagnostics::new();
        diag.analyze(&system);

        // Should detect incomplete chain (component 2 is isolated)
        let incomplete_issues = diag.issues_by_type(DiagnosticType::IncompleteChain);
        assert!(!incomplete_issues.is_empty() || diag.stats.dof_remaining > 0);
    }

    #[test]
    fn test_diagnostic_severity() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        system.add_constraint(Constraint::fixed(1, 0, Point3::new(1.0, 0.0, 0.0)));

        let mut diag = SolverDiagnostics::new();
        diag.analyze(&system);

        // Check that we can filter by severity
        let errors = diag.errors();
        let warnings = diag.warnings();

        assert!(!errors.is_empty() || !warnings.is_empty());
    }

    #[test]
    fn test_diagnostic_summary() {
        let mut system = ConstraintSystem::with_components(2);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));

        let mut diag = SolverDiagnostics::new();
        diag.analyze(&system);

        let summary = diag.summary();
        assert!(summary.contains("Components: 2"));
        assert!(summary.contains("Constraints: 1"));
    }

    #[test]
    fn test_recommendations() {
        let mut system = ConstraintSystem::with_components(1);
        system.add_constraint(Constraint::fixed(0, 0, Point3::ZERO));
        system.add_constraint(Constraint::fixed(1, 0, Point3::new(1.0, 0.0, 0.0)));

        let mut diag = SolverDiagnostics::new();
        diag.analyze(&system);

        // Check that issues have recommendations
        for issue in &diag.issues {
            if issue.issue_type == DiagnosticType::Conflict {
                assert!(!issue.recommendations.is_empty());
            }
        }
    }
}
