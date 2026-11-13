//! Grounded bodies in the rigid-pose solve.
//!
//! A body that a `Fixed` constraint pins in full (position and orientation)
//! is placed at that pose and is not an unknown. Parts mated only to the
//! ground therefore do not interact.

use glam::{DQuat, DVec3};
use rcad_gcs::rigid3d::{Constraint, ConstraintSolver, ConstraintSystem};

/// The ground's pose in these tests: deliberately not the identity.
fn ground_pose() -> (DVec3, DQuat) {
    (
        DVec3::new(1.0, 2.0, 3.0),
        DQuat::from_rotation_z(0.3).normalize(),
    )
}

/// A system with body 0 grounded at [`ground_pose`].
fn grounded_system() -> ConstraintSystem {
    let (t, q) = ground_pose();
    let mut system = ConstraintSystem::with_components(1);
    system.add_constraint(Constraint::fixed_with_orientation(0, 0, t, q));
    system
}

/// Add a body pinned to the ground at `pose` by three coincident point
/// pairs (nine consistent equations for its six freedoms). Returns its
/// component index.
fn pin_to_ground(system: &mut ConstraintSystem, pose: (DVec3, DQuat), first_id: usize) -> usize {
    let (ground_t, ground_q) = ground_pose();
    let body = system.num_components;
    system.num_components += 1;
    let locals = [
        DVec3::ZERO,
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(0.0, 3.0, 1.0),
    ];
    for (k, local) in locals.into_iter().enumerate() {
        // Where the body point sits in the world, seen from the ground frame.
        let world = pose.0 + pose.1.mul_vec3(local);
        let on_ground = ground_q.inverse().mul_vec3(world - ground_t);
        system.add_constraint(Constraint::coincident(
            first_id + k,
            0,
            on_ground,
            body,
            local,
        ));
    }
    body
}

/// Every body at the origin with identity orientation — the ground too.
fn all_at_origin(system: &ConstraintSystem) -> (Vec<DVec3>, Vec<DQuat>) {
    (
        vec![DVec3::ZERO; system.num_components],
        vec![DQuat::IDENTITY; system.num_components],
    )
}

fn bits(v: DVec3) -> [u64; 3] {
    [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()]
}

#[test]
fn a_grounded_body_is_placed_exactly_from_any_start() {
    let mut system = grounded_system();
    let part_pose = (
        DVec3::new(4.0, -1.0, 0.5),
        DQuat::from_rotation_x(0.4).normalize(),
    );
    let part = pin_to_ground(&mut system, part_pose, 10);

    // The warm start has the ground somewhere else entirely.
    let (positions, rots) = all_at_origin(&system);
    let result = ConstraintSolver::new().solve_from(&system, positions, rots);

    assert!(result.success, "error {}", result.error);
    let (ground_t, ground_q) = ground_pose();
    assert_eq!(
        bits(result.positions[0]),
        bits(ground_t),
        "the ground is placed, not approached"
    );
    assert_eq!(result.rotations[0], ground_q);

    assert!((result.positions[part] - part_pose.0).length() < 1e-9);
    assert!(result.rotations[part].angle_between(part_pose.1) < 1e-9);
}

#[test]
fn parts_pinned_to_the_ground_do_not_interact() {
    let pose_a = (
        DVec3::new(4.0, -1.0, 0.5),
        DQuat::from_rotation_x(0.4).normalize(),
    );

    // Part A alone.
    let mut alone = grounded_system();
    let a_alone = pin_to_ground(&mut alone, pose_a, 10);
    let (positions, rots) = all_at_origin(&alone);
    let solved_alone = ConstraintSolver::new().solve_from(&alone, positions, rots);
    assert!(solved_alone.success);

    // Part A next to a part B that cannot be assembled: two of its points
    // are 2 apart, the ground points they must meet are 5 apart.
    let mut both = grounded_system();
    let a = pin_to_ground(&mut both, pose_a, 10);
    let b = both.num_components;
    both.num_components += 1;
    both.add_constraint(Constraint::coincident(20, 0, DVec3::ZERO, b, DVec3::ZERO));
    both.add_constraint(Constraint::coincident(
        21,
        0,
        DVec3::new(5.0, 0.0, 0.0),
        b,
        DVec3::new(2.0, 0.0, 0.0),
    ));
    let (positions, rots) = all_at_origin(&both);
    let solved_both = ConstraintSolver::new().solve_from(&both, positions, rots);

    assert!(!solved_both.success, "part B has no solution");
    assert!(solved_both.error > 1.0, "error {}", solved_both.error);
    // Part A does not notice: its pose is bit for bit what it is alone.
    assert_eq!(
        bits(solved_both.positions[a]),
        bits(solved_alone.positions[a_alone])
    );
    assert_eq!(solved_both.rotations[a], solved_alone.rotations[a_alone]);
}

#[test]
fn a_second_full_fix_at_another_pose_fails_the_solve() {
    let mut system = grounded_system();
    system.add_constraint(Constraint::fixed_with_orientation(
        1,
        0,
        DVec3::new(9.0, 9.0, 9.0),
        DQuat::IDENTITY,
    ));

    let result = ConstraintSolver::new().solve(&system);

    assert!(!result.success);
    assert!(result.error > 1.0, "error {}", result.error);
    // The first full fix is the one that defines the body.
    assert_eq!(bits(result.positions[0]), bits(ground_pose().0));
}

#[test]
fn a_second_full_fix_at_the_same_pose_is_harmless() {
    let (t, q) = ground_pose();
    let mut system = grounded_system();
    system.add_constraint(Constraint::fixed_with_orientation(1, 0, t, q));

    let result = ConstraintSolver::new().solve(&system);

    assert!(result.success, "error {}", result.error);
    assert_eq!(result.error, 0.0);
}

#[test]
fn a_translation_only_fix_is_still_an_equation() {
    // No orientation given: the body keeps three rotational freedoms and is
    // solved for, not placed.
    let target = DVec3::new(1.0, 2.0, 3.0);
    let mut system = ConstraintSystem::with_components(1);
    system.add_constraint(Constraint::fixed(0, 0, target));

    let start_rotation = DQuat::from_rotation_y(0.7).normalize();
    let result = ConstraintSolver::new().solve_from(
        &system,
        vec![DVec3::new(0.5, 0.5, 0.5)],
        vec![start_rotation],
    );

    assert!(result.success, "error {}", result.error);
    assert!((result.positions[0] - target).length() < 1e-9);
    assert!(
        result.rotations[0].angle_between(start_rotation) < 1e-12,
        "nothing asks the body to turn"
    );
}
