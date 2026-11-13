//! Serialized constraint-model round-trip (serde feature only).
#![cfg(feature = "serde")]

use rcad_gcs::prelude::*;

/// The serialized constraint model round-trips through JSON (serde feature).
#[test]
fn constraint_model_json_roundtrip() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 10.0, 0.0);
    sk.add_constraint(Constraint::Horizontal(l));
    sk.add_constraint(Constraint::LineLength {
        line: l,
        length: 10.0,
    });

    let json = serde_json::to_string(&sk).expect("serialize");
    let back: Sketch = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(back.entities.len(), sk.entities.len());
    assert_eq!(back.constraints().len(), sk.constraints().len());

    // The round-tripped sketch still solves to the same geometry.
    let mut sk2 = back;
    let result = sk2.solve();
    assert!(result.converged, "residual {}", result.residual);
}

/// An assembly — one constraint of every kind — round-trips through JSON
/// unchanged, and the copy solves to the same poses.
///
/// Every number here is exactly representable in binary. `serde_json`
/// reads long decimals with a parser that may be one unit in the last place
/// off (unless its `float_roundtrip` feature is on), and this test is about
/// the model, not about that.
#[test]
fn assembly_json_roundtrip() {
    use glam::{DQuat, DVec3};
    use rcad_gcs::rigid3d::{Constraint, ConstraintSolver, ConstraintSystem};

    // Half a turn about Z.
    let turn = DQuat::from_xyzw(0.0, 0.0, 1.0, 0.0);
    let mut system = ConstraintSystem::with_components(4);
    for constraint in [
        Constraint::fixed_with_orientation(10, 0, DVec3::ZERO, DQuat::IDENTITY),
        Constraint::fixed(11, 1, DVec3::new(2.0, 0.0, 0.0)),
        Constraint::coincident(12, 0, DVec3::X, 2, DVec3::ZERO),
        Constraint::distance(13, 0, DVec3::ZERO, 3, DVec3::ZERO, 3.0),
        Constraint::angle(14, 0, DVec3::X, 1, DVec3::X, 0.75),
        Constraint::parallel(15, 0, DVec3::Z, 2, DVec3::Z),
        Constraint::perpendicular(16, 0, DVec3::X, 3, DVec3::Y),
        Constraint::concentric(17, 0, DVec3::ZERO, DVec3::Z, 1, DVec3::ZERO, DVec3::Z),
        Constraint::point_on_plane(18, 0, DVec3::ZERO, DVec3::Z, 3, DVec3::ZERO),
        Constraint::lock(19, 0, 2, DVec3::X, turn),
    ] {
        system.add_constraint(constraint);
    }

    let json = serde_json::to_string(&system).expect("serialize");
    let back: ConstraintSystem = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(back, system);

    let solver = ConstraintSolver::new();
    let (a, b) = (solver.solve(&system), solver.solve(&back));
    assert_eq!(a.positions, b.positions);
    assert_eq!(a.rotations, b.rotations);
}
