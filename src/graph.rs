//! Constraint-graph decomposition (solver roadmap M3).
//!
//! Entities and constraints form a bipartite graph: a constraint is
//! adjacent to every parameter it touches. The connected components of
//! that graph partition the constraint set into subsystems that share no
//! parameters — they can be solved, diagnosed, and skipped independently.
//! This is the structural precondition for DCM-style localized solving:
//! full solve = solving each component, and incremental re-solve = solving
//! only the components a constraint edit landed in.
//!
//! The parameter sets come from [`Constraint::jacobian`], whose rows list
//! every touched parameter (derivatives may be zero there, but the entry
//! is present) — so the decomposition stays in lockstep with the equation
//! set without a third per-variant method.

use crate::constraint::Constraint;
use crate::entity::Entity;
use crate::solver::SketchEquations;
use crate::system;

/// Connected components of the constraint graph.
///
/// Returns groups of constraint indices (positions in `constraints`).
/// Constraints sharing any parameter land in the same group; groups are
/// ordered by first appearance, so the partition is deterministic.
pub fn constraint_components(
    params: &[f64],
    entities: &[Entity],
    constraints: &[Constraint],
) -> Vec<Vec<usize>> {
    system::components(&SketchEquations::new(entities, constraints), params)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::{Entity, EntityKind, PointRef};
    use crate::sketch::Sketch;

    #[test]
    fn disjoint_groups_split_two_components() {
        let mut sk = Sketch::new();
        // Group A: a fixed point + a line pinned to it.
        let pa = sk.add_point(0.0, 0.0);
        let la = sk.add_line(0.0, 0.0, 5.0, 0.0);
        sk.add_constraint(Constraint::fix_point(PointRef::Point(pa), 0.0, 0.0));
        sk.add_constraint(Constraint::coincident(PointRef::LineStart(la), pa));
        sk.add_constraint(Constraint::LineLength {
            line: la,
            length: 5.0,
        });
        // Group B: a self-contained circle.
        let cb = sk.add_circle(9.0, 9.0, 2.0);
        sk.add_constraint(Constraint::Radius {
            circle: cb,
            radius: 2.0,
        });

        let comps = constraint_components(&sk.params, &sk.entities, &sk.constraints);
        assert_eq!(comps.len(), 2, "two disjoint groups: {comps:?}");
        // First-appearance order: group A's three constraints first.
        assert_eq!(comps[0], vec![0, 1, 2]);
        assert_eq!(comps[1], vec![3]);
    }

    #[test]
    fn shared_parameter_merges_components() {
        let entities = vec![
            Entity::new(EntityKind::Point, 0),
            Entity::new(EntityKind::Point, 2),
        ];
        let params = vec![0.0, 0.0, 5.0, 5.0];
        let constraints = vec![
            Constraint::fix_point(PointRef::Point(0), 0.0, 0.0),
            Constraint::point_distance(PointRef::Point(0), PointRef::Point(1), 5.0),
        ];
        let comps = constraint_components(&params, &entities, &constraints);
        assert_eq!(comps.len(), 1, "shared point merges: {comps:?}");
    }

    /// Regression (M6.1, found by the assembly port): a constraint that merges
    /// two established groups must pull EARLIER constraints of both groups
    /// into one component — grouping by the root at insertion time split
    /// connected systems when the merge happened after the first group formed.
    #[test]
    fn late_merge_pulls_earlier_constraints_together() {
        let mut sk = Sketch::new();
        let p0 = sk.add_point(0.0, 0.0);
        let p1 = sk.add_point(1.0, 1.0);
        let p2 = sk.add_point(2.0, 2.0);
        let p3 = sk.add_point(3.0, 3.0);
        sk.add_constraint(Constraint::fix_point(PointRef::Point(p0), 0.0, 0.0)); // c0: p0
        sk.add_constraint(Constraint::coincident(
            PointRef::Point(p1),
            PointRef::Point(p2),
        )); // c1: p1,p2
        sk.add_constraint(Constraint::coincident(
            PointRef::Point(p0),
            PointRef::Point(p1),
        )); // c2: merges
        sk.add_constraint(Constraint::coincident(
            PointRef::Point(p2),
            PointRef::Point(p3),
        )); // c3: p3 joins
        let comps = crate::graph::constraint_components(&sk.params, &sk.entities, &sk.constraints);
        assert_eq!(
            comps.len(),
            1,
            "all four constraints share one component: {comps:?}"
        );
        assert_eq!(comps[0], vec![0, 1, 2, 3]);
    }
}
