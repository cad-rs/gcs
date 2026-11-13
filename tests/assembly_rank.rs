//! A mate contributes as many equations as it removes freedoms.
//!
//! `Parallel` removes two and `Concentric` four. Written with more rows than
//! that, a mate carries a dependency among its own rows, and the diagnosis
//! reports a perfectly sound mate as redundant.

use glam::{DQuat, DVec3};
use rcad_gcs::rigid3d::{
    diagnose_assembly, AssemblyDiagnosis, Constraint, ConstraintSolver, ConstraintSystem,
    SolveOptions,
};

/// Solve tightly, then diagnose at the solved pose.
fn solve_and_diagnose(system: &ConstraintSystem) -> AssemblyDiagnosis {
    let solver = ConstraintSolver::with_options(SolveOptions {
        tolerance: 1.0e-12,
        ..SolveOptions::default()
    });
    let solved = solver.solve(system);
    assert!(solved.success, "error {:.3e}", solved.error);
    diagnose_assembly(system, &solved.positions, &solved.rotations)
}

fn grounded_pair() -> ConstraintSystem {
    let mut system = ConstraintSystem::with_components(2);
    system.add_constraint(Constraint::fixed_with_orientation(
        0,
        0,
        DVec3::ZERO,
        DQuat::IDENTITY,
    ));
    system
}

#[test]
fn a_lone_parallel_mate_is_not_redundant() {
    for local in [
        DVec3::new(0.3, 0.1, 1.0),
        DVec3::new(1.0, -0.4, 0.2),
        DVec3::new(-0.2, 1.0, 0.5),
    ] {
        let mut system = grounded_pair();
        system.add_constraint(Constraint::parallel(1, 0, DVec3::Z, 1, local));
        let d = solve_and_diagnose(&system);
        assert_eq!(d.redundant, Vec::<usize>::new(), "{local}");
        assert_eq!(d.partially_redundant, Vec::<usize>::new(), "{local}");
        assert_eq!(d.conflicting, Vec::<usize>::new(), "{local}");
        assert!(d.dependency_groups.is_empty(), "{local}");
        // Two bodies, the first pinned in full, two freedoms taken.
        assert_eq!(d.rank, 6 + 2);
        assert_eq!(d.dof, 4);
    }
}

#[test]
fn a_lone_concentric_mate_is_not_redundant() {
    for (point, axis) in [
        (DVec3::new(0.5, -0.2, 0.1), DVec3::new(0.2, 0.3, 1.0)),
        (DVec3::new(-1.0, 0.4, 2.0), DVec3::new(1.0, 0.1, -0.3)),
    ] {
        let mut system = grounded_pair();
        system.add_constraint(Constraint::concentric(
            1,
            0,
            DVec3::new(1.0, 2.0, 0.0),
            DVec3::Z,
            1,
            point,
            axis.normalize(),
        ));
        let d = solve_and_diagnose(&system);
        assert_eq!(d.redundant, Vec::<usize>::new(), "{axis}");
        assert_eq!(d.partially_redundant, Vec::<usize>::new(), "{axis}");
        assert_eq!(d.conflicting, Vec::<usize>::new(), "{axis}");
        assert!(d.dependency_groups.is_empty(), "{axis}");
        // Sliding along the axis and spinning about it remain.
        assert_eq!(d.rank, 6 + 4);
        assert_eq!(d.dof, 2);
    }
}

/// Two parallel mates on the same pair of directions are redundant with
/// each other, and that is still reported.
#[test]
fn a_repeated_parallel_mate_is_redundant() {
    let mut system = grounded_pair();
    let local = DVec3::new(0.3, 0.1, 1.0);
    system.add_constraint(Constraint::parallel(1, 0, DVec3::Z, 1, local));
    system.add_constraint(Constraint::parallel(2, 0, DVec3::Z, 1, local));
    let d = solve_and_diagnose(&system);
    assert_eq!(d.conflicting, Vec::<usize>::new());
    assert_eq!(d.rank, 6 + 2);
    assert_eq!(d.dof, 4);
    assert!(
        d.redundant == vec![1] || d.redundant == vec![2],
        "one of the pair: {:?}",
        d.redundant
    );
    assert_eq!(d.dependency_groups, vec![vec![1, 2]]);
}

/// One of every kind of mate between the ground (body 0) and a free body.
fn one_of_each() -> Vec<(&'static str, Constraint)> {
    let turn = DQuat::from_euler(glam::EulerRot::XYZ, 0.4, -0.3, 0.7);
    vec![
        ("fixed", Constraint::fixed(1, 1, DVec3::new(1.0, 2.0, 3.0))),
        (
            "fixed with orientation",
            Constraint::fixed_with_orientation(1, 1, DVec3::new(1.0, 2.0, 3.0), turn),
        ),
        (
            "coincident",
            Constraint::coincident(
                1,
                0,
                DVec3::new(1.0, 0.5, 0.0),
                1,
                DVec3::new(0.2, 0.0, 0.4),
            ),
        ),
        (
            "distance",
            Constraint::distance(
                1,
                0,
                DVec3::new(1.0, 0.5, 0.0),
                1,
                DVec3::new(0.2, 0.0, 0.4),
                2.5,
            ),
        ),
        (
            "angle",
            Constraint::angle(1, 0, DVec3::X, 1, DVec3::new(1.0, 1.0, 0.2), 1.0),
        ),
        (
            "parallel",
            Constraint::parallel(
                1,
                0,
                DVec3::new(0.2, -0.3, 1.0),
                1,
                DVec3::new(0.5, 0.1, 1.0),
            ),
        ),
        (
            "perpendicular",
            Constraint::perpendicular(1, 0, DVec3::X, 1, DVec3::new(1.0, 1.0, 0.2)),
        ),
        (
            "concentric",
            Constraint::concentric(
                1,
                0,
                DVec3::new(1.0, 2.0, 0.0),
                DVec3::new(0.2, -0.3, 1.0),
                1,
                DVec3::new(0.5, -0.2, 0.1),
                DVec3::new(0.5, 0.1, 1.0),
            ),
        ),
        (
            "point on plane",
            Constraint::point_on_plane(
                1,
                0,
                DVec3::new(0.0, 0.0, 1.5),
                DVec3::new(0.1, 0.2, 1.0),
                1,
                DVec3::new(0.3, 0.2, -0.4),
            ),
        ),
        (
            "lock",
            Constraint::lock(1, 0, 1, DVec3::new(0.5, -1.0, 2.0), turn),
        ),
    ]
}

/// Alone, every kind of mate has exactly as many independent equations as
/// it removes freedoms: nothing in it is redundant, and the freedoms left
/// over are six minus what it takes.
#[test]
fn every_mate_alone_takes_exactly_its_freedoms() {
    let mates = one_of_each();
    assert_eq!(mates.len(), 10, "one per kind, `Fixed` in both forms");
    for (name, mate) in mates {
        let mut system = grounded_pair();
        system.add_constraint(mate);
        let d = solve_and_diagnose(&system);
        let removed = mate.dof_removed();
        assert_eq!(d.rank, 6 + removed, "{name}: rank");
        assert_eq!(d.dof, 6 - removed as i64, "{name}: freedoms left");
        assert_eq!(d.redundant, Vec::<usize>::new(), "{name}");
        assert_eq!(d.partially_redundant, Vec::<usize>::new(), "{name}");
        assert_eq!(d.conflicting, Vec::<usize>::new(), "{name}");
        assert!(d.dependency_groups.is_empty(), "{name}: {d:?}");
    }
}
