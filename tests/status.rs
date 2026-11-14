//! One state per entity and per constraint (P4.1). The corpus gate
//! (`entity_statuses_match_the_labels` in `tests/corpus.rs`) checks every
//! entity of the 2D corpus against labels; these cover the states the
//! corpus does not reach.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::constraint::Constraint::*;
use rcad_gcs::entity::PointRef;
use rcad_gcs::{ConstraintStatus, EntityStatus, Sketch};

/// A point fixed at the origin and a line hanging off it with a length.
fn anchored_line() -> (Sketch, usize, usize, usize) {
    let mut sk = Sketch::new();
    let anchor = sk.add_point(0.0, 0.0);
    let line = sk.add_line(0.2, 0.1, 3.0, 1.5);
    sk.add_constraint(Constraint::fix_point(anchor, 0.0, 0.0));
    sk.add_constraint(Coincident(
        PointRef::LineStart(line),
        PointRef::Point(anchor),
    ));
    let length = sk.add_constraint(LineLength { line, length: 4.0 });
    (sk, anchor, line, length)
}

#[test]
fn geometry_not_solved_yet_is_unsolved_and_solved_it_is_not() {
    let (mut sk, anchor, line, length) = anchored_line();
    let before = sk.status();
    // The anchor stands where it is fixed; the line's start does not meet
    // it, nor does its length hold.
    assert_eq!(before.entities[line], Some(EntityStatus::Unsolved));
    assert_eq!(before.entities[anchor], Some(EntityStatus::Unsolved));
    assert_eq!(before.constraints[0], ConstraintStatus::Met);
    assert_eq!(before.constraints[1], ConstraintStatus::Unmet);
    assert_eq!(before.constraints[length], ConstraintStatus::Unmet);

    assert!(sk.solve().converged);
    let after = sk.status();
    // The line can still turn about the anchor.
    assert_eq!(after.entities[anchor], Some(EntityStatus::WellDefined));
    assert_eq!(after.entities[line], Some(EntityStatus::UnderDefined));
    assert!(after
        .constraints
        .iter()
        .all(|&c| c == ConstraintStatus::Met));
    assert_eq!(after.dof, 1);
}

#[test]
fn a_removed_entity_has_no_state_and_a_suppressed_constraint_is_inactive() {
    let (mut sk, _, line, length) = anchored_line();
    let spare = sk.add_point(5.0, 5.0);
    sk.set_constraint_enabled(length, false);
    let removed = sk.remove_entity(spare);
    assert!(removed.is_empty());
    assert!(sk.solve().converged);
    let status = sk.status();
    assert_eq!(status.entities[spare], None);
    assert_eq!(status.constraints[length], ConstraintStatus::Inactive);
    // Without its length the line's end is free in both directions.
    assert_eq!(status.entities[line], Some(EntityStatus::UnderDefined));
}

#[test]
fn every_member_of_a_conflict_is_conflicting_and_the_rest_keeps_its_state() {
    // Two points pinned 1 apart and asked to be 3 apart: the distance and
    // both fixings conflict together. A third point elsewhere is free.
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(1.0, 0.0);
    let c = sk.add_point(4.0, 4.0);
    sk.add_constraint(Constraint::fix_point(a, 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(b, 1.0, 0.0));
    sk.add_constraint(Constraint::point_distance(a, b, 3.0));
    let status = sk.status();
    assert_eq!(status.entities_in(EntityStatus::Conflicting), vec![a, b]);
    assert_eq!(status.entities[c], Some(EntityStatus::UnderDefined));
    assert_eq!(
        status.constraints_in(ConstraintStatus::Conflicting),
        vec![0, 1, 2]
    );
}

#[test]
fn a_redundant_dimension_overdefines_only_its_own_geometry() {
    let (mut sk, anchor, line, length) = anchored_line();
    sk.add_constraint(Horizontal(line));
    // The same length again, as a horizontal distance: redundant with the
    // length once the line is horizontal.
    let again = sk.add_constraint(Constraint::distance_x(
        PointRef::LineEnd(line),
        PointRef::LineStart(line),
        4.0,
    ));
    let other = sk.add_point(9.0, 9.0);
    sk.add_constraint(Constraint::fix_point(other, 9.0, 9.0));
    assert!(sk.solve().converged);
    let status = sk.status();
    assert_eq!(status.entities[line], Some(EntityStatus::OverDefined));
    assert_eq!(status.entities[other], Some(EntityStatus::WellDefined));
    assert_eq!(status.constraints[length], ConstraintStatus::Redundant);
    assert_eq!(status.constraints[again], ConstraintStatus::Redundant);
    // The dependency is between the two dimensions of the line alone.
    assert_eq!(status.entities[anchor], Some(EntityStatus::WellDefined));
    assert_eq!(status.dependency_groups, vec![vec![length, again]]);
}

#[test]
fn the_states_are_ordered_by_seriousness() {
    use EntityStatus::*;
    assert!(WellDefined < UnderDefined);
    assert!(UnderDefined < OverDefined);
    assert!(OverDefined < Unsolved);
    assert!(Unsolved < Conflicting);
    assert!(ConstraintStatus::Met < ConstraintStatus::Redundant);
    assert!(ConstraintStatus::Unmet < ConstraintStatus::Conflicting);
}

// ── 3D ──────────────────────────────────────────────────────────────────────

mod space {
    use rcad_gcs::space3d::constraint::SpaceConstraint;
    use rcad_gcs::space3d::entity::SpacePointRef;
    use rcad_gcs::space3d::SpaceSketch;
    use rcad_gcs::EntityStatus;

    #[test]
    fn a_pinned_plane_is_well_defined_whatever_its_normal_length() {
        // Four stored numbers, three freedoms: the normal's length is how
        // the plane is stored, not a way it can move.
        let mut sk = SpaceSketch::new();
        let plane = sk.add_plane(0.1, -0.2, 1.3, 0.4);
        let p = sk.add_point(0.0, 0.0, 2.0);
        sk.fix_entity(p);
        sk.add_constraint(SpaceConstraint::PlaneNormal {
            plane,
            nx: 0.0,
            ny: 0.0,
            nz: 1.0,
        });
        sk.add_constraint(SpaceConstraint::PointOnPlane {
            point: SpacePointRef::Point(p),
            plane,
        });
        assert!(sk.solve().converged);
        let status = sk.status();
        assert_eq!(
            status.entities[plane],
            Some(EntityStatus::WellDefined),
            "{status:?}"
        );
        assert_eq!(status.dof, 0);

        let mut free = SpaceSketch::new();
        let plane = free.add_plane(0.3, -0.4, 0.5, 1.2);
        assert_eq!(
            free.status().entities[plane],
            Some(EntityStatus::UnderDefined)
        );
    }

    #[test]
    fn a_cylinder_through_five_points_is_well_defined_though_its_axis_point_slides() {
        // Five held points on the surface of a cylinder of radius 1 about
        // the z axis fix its five freedoms; sliding the stored axis point
        // along the axis, or rescaling the axis, moves nothing.
        let mut sk = SpaceSketch::new();
        let c = sk.add_cylinder(0.0, 0.0, 0.5, 0.0, 0.0, 2.0, 1.0);
        for (x, y, z) in [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 1.0),
            (-1.0, 0.0, 2.0),
            (0.0, -1.0, 3.0),
            (0.6, 0.8, 4.0),
        ] {
            let p = sk.add_point(x, y, z);
            sk.fix_entity(p);
            sk.add_constraint(SpaceConstraint::PointOnCylinder {
                point: SpacePointRef::Point(p),
                cylinder: c,
            });
        }
        assert!(sk.solve().converged);
        let status = sk.status();
        assert_eq!(status.dof, 0, "{status:?}");
        assert_eq!(
            status.entities[c],
            Some(EntityStatus::WellDefined),
            "{status:?}"
        );

        // With one point fewer it has a freedom left.
        let mut fewer = SpaceSketch::new();
        let c = fewer.add_cylinder(0.0, 0.0, 0.5, 0.0, 0.0, 2.0, 1.0);
        for (x, y, z) in [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 1.0),
            (-1.0, 0.0, 2.0),
            (0.0, -1.0, 3.0),
        ] {
            let p = fewer.add_point(x, y, z);
            fewer.fix_entity(p);
            fewer.add_constraint(SpaceConstraint::PointOnCylinder {
                point: SpacePointRef::Point(p),
                cylinder: c,
            });
        }
        assert_eq!(fewer.status().entities[c], Some(EntityStatus::UnderDefined));
    }
}
