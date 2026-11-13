//! Constraint-graph decomposition for 3D space sketches.
//!
//! Mirrors the 2D [`crate::graph`] module: a constraint is adjacent to every
//! parameter it touches, and the connected components of that bipartite
//! graph partition the constraint set into subsystems that share no
//! parameters. Components can be solved, diagnosed, and skipped
//! independently — the structural precondition for localized (incremental)
//! re-solve.
//!
//! The touch sets come from [`SpaceConstraint::jacobian`], whose rows list
//! every touched parameter (derivatives may be zero there, but the entry is
//! present), so the decomposition stays in lockstep with the equation set.

use std::collections::HashMap;

use crate::space3d::constraint::SpaceConstraint;
use crate::space3d::entity::SpaceEntity;

/// Connected components of the 3D constraint graph.
///
/// Returns groups of constraint indices (positions in `constraints`).
/// Constraints sharing any parameter land in the same group; groups are
/// ordered by first appearance, so the partition is deterministic.
pub fn constraint_components(
    params: &[f64],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
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
    // stale roots when a later constraint merges two established groups.
    let mut order: Vec<usize> = Vec::new();
    let mut by_root: HashMap<usize, Vec<usize>> = HashMap::new();
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
            // satisfy through the solver; it forms its own (trivial) group
            // so indices are never dropped.
            None => {
                order.push(usize::MAX - idx);
                by_root.insert(usize::MAX - idx, vec![idx]);
            }
        }
    }
    order
        .into_iter()
        .filter_map(|root| by_root.remove(&root))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space3d::constraint::SpaceConstraint;
    use crate::space3d::entity::{SpaceEntityKind, SpacePointRef};
    use crate::space3d::sketch::SpaceSketch;

    #[test]
    fn disjoint_groups_split_two_components() {
        let mut sk = SpaceSketch::new();
        // Group A: a fixed point + a line pinned to it.
        let pa = sk.add_point(0.0, 0.0, 0.0);
        let la = sk.add_line(0.0, 0.0, 0.0, 5.0, 0.0, 0.0);
        sk.add_constraint(SpaceConstraint::fix_point(pa, 0.0, 0.0, 0.0));
        sk.add_constraint(SpaceConstraint::coincident(SpacePointRef::LineStart(la), pa));
        sk.add_constraint(SpaceConstraint::LineLength { line: la, length: 5.0 });
        // Group B: a self-contained sphere radius.
        let sb = sk.add_sphere(9.0, 9.0, 9.0, 2.0);
        sk.add_constraint(SpaceConstraint::SphereRadius { sphere: sb, radius: 2.0 });

        let comps = constraint_components(&sk.params, &sk.entities, &sk.constraints);
        assert_eq!(comps.len(), 2, "two disjoint groups: {comps:?}");
        assert_eq!(comps[0], vec![0, 1, 2]);
        assert_eq!(comps[1], vec![3]);
    }

    /// A curve constraint joins the line and the circle into one component;
    /// an unrelated curve stays in its own component. The touch sets come
    /// from the new circle rows, so this also guards their completeness.
    #[test]
    fn curve_constraints_join_components() {
        use crate::space3d::entity::SpaceFrame3;
        let mut sk = SpaceSketch::new();
        let line = sk.add_line(0.0, 0.0, 0.0, 3.0, 0.0, 0.0);
        let circle =
            sk.add_circle(0.0, 1.0, 0.0, 1.0, SpaceFrame3::sanitized(glam::DVec3::Z, glam::DVec3::X, 0.0, 0.0));
        let other =
            sk.add_circle(9.0, 9.0, 9.0, 2.0, SpaceFrame3::sanitized(glam::DVec3::Y, glam::DVec3::X, 0.0, 0.0));
        sk.add_constraint(SpaceConstraint::LineTangentCircle { line, circle });
        sk.add_constraint(SpaceConstraint::CircleRadius { circle: other, radius: 2.0 });
        let comps = constraint_components(&sk.params, &sk.entities, &sk.constraints);
        assert_eq!(comps.len(), 2, "line+circle and the unrelated circle: {comps:?}");
        assert_eq!(comps[0], vec![0]);
        assert_eq!(comps[1], vec![1]);
    }

    #[test]
    fn shared_parameter_merges_components() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 3),
        ];
        let params = vec![0.0, 0.0, 0.0, 5.0, 5.0, 5.0];
        let constraints = vec![
            SpaceConstraint::fix_point(SpacePointRef::Point(0), 0.0, 0.0, 0.0),
            SpaceConstraint::point_distance(SpacePointRef::Point(0), SpacePointRef::Point(1), 5.0),
        ];
        let comps = constraint_components(&params, &entities, &constraints);
        assert_eq!(comps.len(), 1, "shared point merges: {comps:?}");
    }

    /// A constraint that merges two established groups must pull EARLIER
    /// constraints of both groups into one component.
    #[test]
    fn late_merge_pulls_earlier_constraints_together() {
        let mut sk = SpaceSketch::new();
        let p0 = sk.add_point(0.0, 0.0, 0.0);
        let p1 = sk.add_point(1.0, 1.0, 1.0);
        let p2 = sk.add_point(2.0, 2.0, 2.0);
        let p3 = sk.add_point(3.0, 3.0, 3.0);
        sk.add_constraint(SpaceConstraint::fix_point(p0, 0.0, 0.0, 0.0));
        sk.add_constraint(SpaceConstraint::coincident(p1, p2));
        sk.add_constraint(SpaceConstraint::coincident(p0, p1));
        sk.add_constraint(SpaceConstraint::coincident(p2, p3));
        let comps = constraint_components(&sk.params, &sk.entities, &sk.constraints);
        assert_eq!(comps.len(), 1, "all four constraints share one component: {comps:?}");
        assert_eq!(comps[0], vec![0, 1, 2, 3]);
    }
}
