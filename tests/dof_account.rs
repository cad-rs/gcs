//! One DOF account: `dof()` is the rank-based count that `diagnose()`
//! reports, and the cheap count survives as `structural_dof_estimate()`.

use rcad_gcs::constraint::Constraint;
use rcad_gcs::space3d::constraint::SpaceConstraint;
use rcad_gcs::space3d::SpaceSketch;
use rcad_gcs::Sketch;

#[test]
fn redundant_rows_do_not_drive_dof_negative() {
    // One point pinned twice at the same place: four equations, rank two.
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 4.0);
    sk.add_constraint(Constraint::fix_point(p, 3.0, 4.0));
    sk.add_constraint(Constraint::fix_point(p, 3.0, 4.0));

    assert_eq!(sk.dof(), 0, "the point has no freedom left");
    assert_eq!(sk.dof(), sk.diagnose().dof, "one account");
    assert_eq!(
        sk.structural_dof_estimate(),
        -2,
        "2 parameters − 4 equations"
    );

    let report = sk.dof_analysis();
    assert_eq!(report.free_params, 2);
    assert_eq!(report.equations, 4);
    assert_eq!(report.rank, 2);
    assert_eq!(report.dof, 0);
    assert_eq!(report.structural_dof, -2);
}

#[test]
fn a_surplus_in_one_place_does_not_hide_a_freedom_elsewhere() {
    // The doubly pinned point has two surplus equations; the free point has
    // two freedoms. The count says 0, the rank says 2.
    let mut sk = Sketch::new();
    let pinned = sk.add_point(0.0, 0.0);
    let _free = sk.add_point(5.0, 5.0);
    sk.add_constraint(Constraint::fix_point(pinned, 0.0, 0.0));
    sk.add_constraint(Constraint::fix_point(pinned, 0.0, 0.0));

    assert_eq!(sk.structural_dof_estimate(), 0);
    assert_eq!(sk.dof(), 2);
}

#[test]
fn independent_rows_agree_with_the_count() {
    let mut sk = Sketch::new();
    let line = sk.add_line(0.0, 0.0, 3.0, 4.0);
    sk.add_constraint(Constraint::LineLength { line, length: 5.0 });

    assert_eq!(sk.dof(), 3);
    assert_eq!(sk.structural_dof_estimate(), 3);
}

// ── 3D: storage parameters that are not physical freedoms ───────────────────

/// A cylinder with a deliberately non-unit, non-axis-aligned direction.
fn tilted_cylinder(sk: &mut SpaceSketch) -> usize {
    sk.add_cylinder(1.0, -2.0, 0.5, 0.3, 0.2, 1.0, 2.0)
}

#[test]
fn free_cylinder_has_five_degrees_of_freedom() {
    // Seven storage parameters; sliding the axis point along the axis and
    // rescaling the direction do not change the cylinder.
    let mut sk = SpaceSketch::new();
    tilted_cylinder(&mut sk);

    assert_eq!(sk.structural_dof_estimate(), 7);
    assert_eq!(sk.dof(), 5);
}

#[test]
fn cylinder_radius_takes_one_freedom() {
    let mut sk = SpaceSketch::new();
    let cylinder = tilted_cylinder(&mut sk);
    sk.add_constraint(SpaceConstraint::CylinderRadius {
        cylinder,
        radius: 2.0,
    });

    assert_eq!(sk.dof(), 4, "axis line (4) remains: {:?}", sk.diagnose());
}

#[test]
fn point_on_cylinder_takes_one_freedom() {
    let mut sk = SpaceSketch::new();
    let cylinder = tilted_cylinder(&mut sk);
    let point = sk.add_point(3.5, -1.6, 1.5);
    sk.fix_entity(point);
    sk.add_constraint(SpaceConstraint::PointOnCylinder {
        point: point.into(),
        cylinder,
    });

    let d = sk.diagnose();
    assert_eq!(d.rank, 1, "{d:?}");
    assert_eq!(d.dof, 4, "a cylinder through a fixed point: {d:?}");
}

#[test]
fn concentric_pins_the_axis_point_including_its_slide() {
    // Concentric3D reads the cylinder's axis point as a centre, so the slide
    // along the axis is no longer free; only the direction rescale is gauge.
    let mut sk = SpaceSketch::new();
    let cylinder = tilted_cylinder(&mut sk);
    let sphere = sk.add_sphere(1.0, -2.0, 0.5, 1.0);
    sk.fix_entity(sphere);
    sk.add_constraint(SpaceConstraint::Concentric3D(cylinder, sphere));

    let d = sk.diagnose();
    assert_eq!(d.rank, 3, "{d:?}");
    assert_eq!(d.dof, 3, "axis direction (2) + radius (1): {d:?}");
    assert!(d.conflicting.is_empty() && d.redundant.is_empty(), "{d:?}");
}

#[test]
fn fixed_cylinder_has_no_freedom_and_no_gauge() {
    let mut sk = SpaceSketch::new();
    let cylinder = tilted_cylinder(&mut sk);
    sk.fix_entity(cylinder);

    assert_eq!(sk.dof(), 0);
}

#[test]
fn free_plane_still_has_three() {
    let mut sk = SpaceSketch::new();
    sk.add_plane(0.0, 0.0, 2.0, 1.0);

    assert_eq!(sk.structural_dof_estimate(), 4);
    assert_eq!(sk.dof(), 3);
}
