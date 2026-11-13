//! Orientation residuals of the rigid-pose family: `Fixed` with an
//! orientation and `Lock`.
//!
//! A unit quaternion and its negative are the same rotation. Which of the
//! two a host hands over is an accident of how it converted its matrices,
//! and must not matter.

use glam::{DQuat, DVec3};
use rcad_gcs::rigid3d::{
    diagnose_assembly, Constraint, ConstraintSolver, ConstraintSystem, PoseSolveOptions,
};

fn negated(q: DQuat) -> DQuat {
    DQuat::from_xyzw(-q.x, -q.y, -q.z, -q.w)
}

fn solver() -> ConstraintSolver {
    ConstraintSolver::with_options(PoseSolveOptions {
        tolerance: 1.0e-12,
        ..PoseSolveOptions::default()
    })
}

/// A ground and a part locked to it at a relative pose.
fn locked(translation: DVec3, rotation: DQuat) -> ConstraintSystem {
    let mut system = ConstraintSystem::with_components(2);
    system.add_constraint(Constraint::fixed_with_orientation(
        0,
        0,
        DVec3::ZERO,
        DQuat::IDENTITY,
    ));
    system.add_constraint(Constraint::lock(1, 0, 1, translation, rotation));
    system
}

/// The part already sits exactly where the lock wants it, with its
/// quaternion written in the other sign. Nothing is wrong and nothing
/// should move.
#[test]
fn a_lock_that_holds_is_satisfied_in_either_quaternion_sign() {
    let translation = DVec3::new(1.0, -2.0, 0.5);
    for rotation in [
        DQuat::from_rotation_z(0.4),
        DQuat::from_axis_angle(DVec3::new(1.0, 2.0, -1.0).normalize(), 2.9),
        DQuat::IDENTITY,
    ] {
        let system = locked(translation, rotation);
        for start in [rotation, negated(rotation)] {
            let result = solver().solve_from(
                &system,
                vec![DVec3::ZERO, translation],
                vec![DQuat::IDENTITY, start],
            );
            assert!(result.success, "error {:.3e} from {start:?}", result.error);
            assert!(result.error < 1e-12, "error {:.3e}", result.error);
            assert_eq!(result.iterations, 0, "there is nothing to do");
            assert_eq!(result.positions[1], translation);
            assert_eq!(result.rotations[1], start, "the part is left alone");
        }
    }
}

/// The same for a body held at a pose by a second, full `Fixed` (the first
/// one grounds it; the second is an equation that has to hold).
#[test]
fn a_full_fix_that_holds_is_satisfied_in_either_quaternion_sign() {
    let position = DVec3::new(0.5, 1.5, -1.0);
    let orientation = DQuat::from_axis_angle(DVec3::new(0.0, 1.0, 1.0).normalize(), 1.3);
    let mut system = ConstraintSystem::with_components(1);
    system.add_constraint(Constraint::fixed_with_orientation(
        0,
        0,
        position,
        orientation,
    ));
    system.add_constraint(Constraint::fixed_with_orientation(
        1,
        0,
        position,
        negated(orientation),
    ));
    let result = solver().solve(&system);
    assert!(result.success, "error {:.3e}", result.error);
    assert!(result.error < 1e-12, "error {:.3e}", result.error);
}

/// A lock is reached from a start most of half a turn away, whichever sign
/// the start is written in.
#[test]
fn a_lock_is_reached_from_far_away_in_either_sign() {
    let translation = DVec3::new(1.0, -2.0, 0.5);
    let rotation = DQuat::from_rotation_z(0.4);
    let system = locked(translation, rotation);
    // 2.5 rad away from the target about Z.
    let away = DQuat::from_rotation_z(0.4 + 2.5);
    for start in [away, negated(away)] {
        let result = solver().solve_from(
            &system,
            vec![DVec3::ZERO, DVec3::new(0.0, 1.0, 0.0)],
            vec![DQuat::IDENTITY, start],
        );
        assert!(result.success, "error {:.3e} from {start:?}", result.error);
        assert!(
            result.rotations[1].angle_between(rotation) < 1e-9,
            "ends at {:?}",
            result.rotations[1]
        );
        assert!((result.positions[1] - translation).length() < 1e-9);
        let d = diagnose_assembly(&system, &result.positions, &result.rotations);
        assert_eq!(d.dof, 0);
        assert!(d.conflicting.is_empty());
    }
}

/// The error of a lock that is off by an angle is that angle, whichever
/// sign the part's quaternion is written in — not the rest of the turn.
#[test]
fn the_error_of_a_lock_is_the_short_rotation_in_either_sign() {
    use rcad_gcs::FailurePolicy;
    let translation = DVec3::new(1.0, -2.0, 0.5);
    let rotation = DQuat::from_rotation_z(0.4);
    let system = locked(translation, rotation);
    // Report the start without iterating.
    let look = ConstraintSolver::with_options(PoseSolveOptions {
        max_iterations: 0,
        on_failure: FailurePolicy::KeepLast,
        ..PoseSolveOptions::default()
    });
    for off in [1.0e-3, 0.5, 2.5, 3.1] {
        let turned = DQuat::from_rotation_z(0.4 + off);
        for start in [turned, negated(turned)] {
            let at_start = look.solve_from(
                &system,
                vec![DVec3::ZERO, translation],
                vec![DQuat::IDENTITY, start],
            );
            assert!(
                (at_start.error - off).abs() < 1e-9,
                "{off} off with w = {:+.3} reads as {}",
                start.w,
                at_start.error
            );
        }
    }
}

/// Two constraints that ask a part for the same orientation, one of them
/// with the quaternion written in the other sign, agree with each other.
/// At a pose a hair off both, that is redundancy, not a conflict.
#[test]
fn one_orientation_written_two_ways_does_not_conflict() {
    let translation = DVec3::new(1.0, -2.0, 0.5);
    let rotation = DQuat::from_rotation_z(0.4);
    let mut system = locked(translation, rotation);
    system.add_constraint(Constraint::fixed_with_orientation(
        2,
        1,
        translation,
        negated(rotation),
    ));

    let nearly = DQuat::from_rotation_z(0.4 + 1.0e-6);
    for start in [nearly, negated(nearly)] {
        let d = diagnose_assembly(
            &system,
            &[DVec3::ZERO, translation],
            &[DQuat::IDENTITY, start],
        );
        assert_eq!(d.dof, 0);
        assert_eq!(d.conflicting, Vec::<usize>::new(), "w = {:+.3}", start.w);
        // Eighteen equations over twelve freedoms: six of them repeat what
        // the others say.
        assert_eq!(d.rank, 12);
        assert!(!d.redundant.is_empty());
    }
}
