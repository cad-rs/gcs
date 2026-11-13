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
use crate::entity::PointRef;
use crate::sketch::Sketch;
use crate::entity::Entity;
use std::collections::HashMap;

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
    let n_params = params.len();
    // Union-find over parameters.
    let mut parent: Vec<usize> = (0..n_params).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        let mut root = x;
        while parent[root] != root {
            root = parent[root];
        }
        let mut cur = x;
        while parent[cur] != root {
            let next = parent[cur];
            parent[cur] = root;
            cur = next;
        }
        root
    }
    fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[ra] = rb;
        }
    }

    // Pass 1: every constraint unions the parameters it touches (from its
    // Jacobian rows — the complete touched set regardless of derivative
    // values).
    let mut scratch: Vec<Vec<(usize, f64)>> = Vec::new();
    for c in constraints {
        let cnt = c.equation_count();
        scratch.clear();
        scratch.resize_with(cnt, Vec::new);
        c.jacobian(params, entities, &mut scratch);
        let mut first: Option<usize> = None;
        for row in &scratch {
            for &(p, _) in row {
                if let Some(f) = first {
                    union(&mut parent, f, p);
                } else {
                    first = Some(p);
                }
            }
        }
    }

    // Pass 2: group constraints by the FINAL root of their first touched
    // parameter. Grouping during pass 1 would key earlier constraints by
    // stale roots when a later constraint merges two established groups
    // (regression-tested).
    let mut order: Vec<usize> = Vec::new();
    let mut by_root: HashMap<usize, Vec<usize>> = HashMap::new();
    scratch.clear();
    for (idx, c) in constraints.iter().enumerate() {
        let cnt = c.equation_count();
        scratch.clear();
        scratch.resize_with(cnt, Vec::new);
        c.jacobian(params, entities, &mut scratch);
        let root = scratch
            .first()
            .and_then(|row| row.first())
            .map(|&(p, _)| find(&mut parent, p));
        match root {
            Some(root) => {
                let group = by_root.entry(root).or_insert_with(|| {
                    order.push(root);
                    Vec::new()
                });
                group.push(idx);
            }
            // A constraint touching no parameters has no equations to
            // satisfy through the solver; it forms its own (trivial)
            // group so indices are never dropped.
            None => {
                order.push(usize::MAX - idx);
                by_root.insert(usize::MAX - idx, vec![idx]);
            }
        }
    }
    order.into_iter()
        .filter_map(|root| by_root.remove(&root))
        .collect()
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
        sk.add_constraint(Constraint::LineLength { line: la, length: 5.0 });
        // Group B: a self-contained circle.
        let cb = sk.add_circle(9.0, 9.0, 2.0);
        sk.add_constraint(Constraint::Radius { circle: cb, radius: 2.0 });

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
    sk.add_constraint(Constraint::coincident(PointRef::Point(p1), PointRef::Point(p2))); // c1: p1,p2
    sk.add_constraint(Constraint::coincident(PointRef::Point(p0), PointRef::Point(p1))); // c2: merges
    sk.add_constraint(Constraint::coincident(PointRef::Point(p2), PointRef::Point(p3))); // c3: p3 joins
    let comps = crate::graph::constraint_components(&sk.params, &sk.entities, &sk.constraints);
    assert_eq!(comps.len(), 1, "all four constraints share one component: {comps:?}");
    assert_eq!(comps[0], vec![0, 1, 2, 3]);
}
