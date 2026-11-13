//! Editing an assembly in place: removing and suppressing constraints,
//! removing components, changing dimension values, measuring, checking.
//!
//! The rigid-pose counterpart of `tests/editing.rs`. Constraint indices and
//! component ids are stable throughout: nothing an edit does shifts the
//! index of something that is still there.

use glam::{DQuat, DVec3};
use rcad_gcs::rigid3d::{
    assembly_mobility, diagnose_assembly, AssemblyError, Constraint, ConstraintSolver,
    ConstraintSystem, PoseSolveOptions, PoseSolveResult,
};
use rcad_gcs::ConstraintState;

const GROUND: usize = 0;
const DOOR: usize = 1;
const SLIDER: usize = 2;

/// A ground, a door hinged on it about Z and held open by an angle
/// dimension, and a slider whose origin rides the ground's X axis at a
/// distance dimension from the origin (its orientation is left free).
/// Returns the system and the indices of the angle and the distance.
fn assembly(angle: f64, distance: f64) -> (ConstraintSystem, usize, usize) {
    let mut s = ConstraintSystem::with_components(3);
    s.add_constraint(Constraint::fixed_with_orientation(
        0,
        GROUND,
        DVec3::ZERO,
        DQuat::IDENTITY,
    ));
    // The hinge: the door's own z axis on the ground's, and no sliding
    // along it.
    s.add_constraint(Constraint::concentric(
        1,
        GROUND,
        DVec3::ZERO,
        DVec3::Z,
        DOOR,
        DVec3::ZERO,
        DVec3::Z,
    ));
    s.add_constraint(Constraint::point_on_plane(
        2,
        GROUND,
        DVec3::ZERO,
        DVec3::Z,
        DOOR,
        DVec3::ZERO,
    ));
    let angle_dim = s.add_constraint(Constraint::angle(
        3,
        GROUND,
        DVec3::X,
        DOOR,
        DVec3::X,
        angle,
    ));
    // The slider's origin on the X axis.
    s.add_constraint(Constraint::point_on_plane(
        4,
        GROUND,
        DVec3::ZERO,
        DVec3::Y,
        SLIDER,
        DVec3::ZERO,
    ));
    s.add_constraint(Constraint::point_on_plane(
        5,
        GROUND,
        DVec3::ZERO,
        DVec3::Z,
        SLIDER,
        DVec3::ZERO,
    ));
    let distance_dim = s.add_constraint(Constraint::distance(
        6,
        GROUND,
        DVec3::ZERO,
        SLIDER,
        DVec3::ZERO,
        distance,
    ));
    (s, angle_dim, distance_dim)
}

/// Where the parts are before the first solve: the door ajar and the
/// slider somewhere off the axis. (An angle between directions that start
/// out parallel has no gradient to follow.)
fn start() -> (Vec<DVec3>, Vec<DQuat>) {
    (
        vec![DVec3::ZERO, DVec3::ZERO, DVec3::new(1.0, 0.2, -0.1)],
        vec![
            DQuat::IDENTITY,
            DQuat::from_rotation_z(0.3),
            DQuat::from_rotation_y(0.2),
        ],
    )
}

fn solver() -> ConstraintSolver {
    ConstraintSolver::with_options(PoseSolveOptions {
        tolerance: 1.0e-12,
        ..PoseSolveOptions::default()
    })
}

fn solve_from_start(system: &ConstraintSystem) -> PoseSolveResult {
    let (positions, rotations) = start();
    let result = solver().solve_from(system, positions, rotations);
    assert!(result.success, "error {:.3e}", result.error);
    result
}

fn resolve(system: &ConstraintSystem, from: &PoseSolveResult) -> PoseSolveResult {
    solver().solve_from(system, from.positions.clone(), from.rotations.clone())
}

/// How far the door stands open.
fn door_angle(result: &PoseSolveResult) -> f64 {
    (result.rotations[DOOR] * DVec3::X)
        .dot(DVec3::X)
        .clamp(-1.0, 1.0)
        .acos()
}

fn dof(system: &ConstraintSystem, at: &PoseSolveResult) -> i64 {
    diagnose_assembly(system, &at.positions, &at.rotations).dof
}

#[test]
fn add_constraint_returns_stable_indices() {
    let (mut s, angle_dim, distance_dim) = assembly(0.5, 2.0);
    assert_eq!((angle_dim, distance_dim), (3, 6));

    assert!(s.remove_constraint(angle_dim));
    // The slots after the removed one keep their indices, and their ids.
    assert_eq!(s.constraints.len(), 7);
    assert_eq!(s.constraint_state(angle_dim), ConstraintState::Removed);
    assert_eq!(s.constraint_state(distance_dim), ConstraintState::Active);
    assert_eq!(s.constraints[distance_dim].id, 6);
    assert_eq!(s.active_constraints().count(), 6);
    // A new constraint gets a new index, not the freed one.
    let again = s.add_constraint(Constraint::angle(7, GROUND, DVec3::X, DOOR, DVec3::X, 0.5));
    assert_eq!(again, 7);
    assert!(!s.remove_constraint(angle_dim), "already removed");
    assert!(s.is_active(again));
}

#[test]
fn a_removed_constraint_is_no_longer_enforced_or_diagnosed() {
    let (mut s, angle_dim, _) = assembly(0.5, 2.0);
    let solved = solve_from_start(&s);
    assert!((door_angle(&solved) - 0.5).abs() < 1e-9);
    // The slider's orientation is all that is left free.
    assert_eq!(dof(&s, &solved), 3);

    assert!(s.remove_constraint(angle_dim));
    assert_eq!(dof(&s, &solved), 4, "the door swings again");
    assert!(
        ConstraintSolver::components(&s)
            .iter()
            .flatten()
            .all(|&c| c != angle_dim),
        "a removed constraint belongs to no cluster"
    );

    // Swing the door by hand: nothing pulls it back to 0.5.
    let mut moved = solved.clone();
    moved.rotations[DOOR] = DQuat::from_rotation_z(1.25);
    let result = resolve(&s, &moved);
    assert!(result.success);
    assert!((door_angle(&result) - 1.25).abs() < 1e-9);
    let d = diagnose_assembly(&s, &result.positions, &result.rotations);
    assert!(d.redundant.is_empty() && d.conflicting.is_empty());
}

#[test]
fn a_suppressed_dimension_is_a_reference_dimension() {
    let (mut s, angle_dim, _) = assembly(0.5, 2.0);
    let solved = solve_from_start(&s);

    // Stop driving the angle and swing the door by hand.
    assert!(s.set_constraint_enabled(angle_dim, false));
    assert_eq!(s.constraint_state(angle_dim), ConstraintState::Suppressed);
    let mut moved = solved.clone();
    moved.rotations[DOOR] = DQuat::from_rotation_z(1.25);
    let result = resolve(&s, &moved);
    assert!(result.success);

    // The dimension still says 0.5, but it measures 1.25 and holds nothing.
    assert_eq!(s.constraint_value(angle_dim), Some(0.5));
    let measured = s
        .measure(angle_dim, &result.positions, &result.rotations)
        .expect("an angle can be measured");
    assert!((measured - 1.25).abs() < 1e-9, "measured {measured}");

    // Driving it again pulls the door back.
    assert!(s.set_constraint_enabled(angle_dim, true));
    let result = resolve(&s, &result);
    assert!(result.success);
    assert!((door_angle(&result) - 0.5).abs() < 1e-9);

    // A removed constraint cannot be re-enabled.
    s.remove_constraint(angle_dim);
    assert!(!s.set_constraint_enabled(angle_dim, true));
    assert_eq!(s.constraint_state(angle_dim), ConstraintState::Removed);
}

#[test]
fn changing_a_dimension_re_solves_from_where_the_assembly_stands() {
    let (mut s, angle_dim, distance_dim) = assembly(0.5, 2.0);
    let original = solve_from_start(&s);
    let revision = s.topology_revision();

    // Open the door and pull the slider out, then put both back: a round
    // trip through warm-started solves, with no system rebuilt in between.
    assert!(s.set_constraint_value(angle_dim, 1.2));
    assert!(s.set_constraint_value(distance_dim, 3.5));
    let moved = resolve(&s, &original);
    assert!(moved.success);
    assert!((door_angle(&moved) - 1.2).abs() < 1e-9);
    assert!((moved.positions[SLIDER] - DVec3::new(3.5, 0.0, 0.0)).length() < 1e-9);
    let measured = s
        .measure(distance_dim, &moved.positions, &moved.rotations)
        .unwrap();
    assert!((measured - 3.5).abs() < 1e-9);

    assert!(s.set_constraint_value(angle_dim, 0.5));
    assert!(s.set_constraint_value(distance_dim, 2.0));
    let back = resolve(&s, &moved);
    assert!(back.success);

    let worst = (0..3)
        .map(|body| {
            let turned = back.rotations[body].angle_between(original.rotations[body]);
            (back.positions[body] - original.positions[body])
                .length()
                .max(turned)
        })
        .fold(0.0_f64, f64::max);
    assert!(worst <= 1e-9, "round-trip drift {worst:.3e}");
    assert_eq!(
        s.topology_revision(),
        revision,
        "changing values is not a structural change"
    );

    // A relation has no value to set.
    assert!(!s.set_constraint_value(1, 1.0), "Concentric has no value");
    assert_eq!(s.constraint_value(1), None);
    assert_eq!(s.measure(1, &back.positions, &back.rotations), None);
}

#[test]
fn removing_a_component_takes_its_constraints_with_it() {
    let (mut s, angle_dim, distance_dim) = assembly(0.5, 2.0);
    let solved = solve_from_start(&s);
    assert_eq!(dof(&s, &solved), 3);

    let gone = s.remove_component(SLIDER);
    assert_eq!(gone, vec![4, 5, distance_dim]);
    assert!(s.is_component_removed(SLIDER));
    for &index in &gone {
        assert_eq!(s.constraint_state(index), ConstraintState::Removed);
    }
    assert_eq!(s.constraint_state(angle_dim), ConstraintState::Active);
    assert_eq!(s.live_components(), 2);
    assert_eq!(s.max_dof(), 12);

    // The removed part counts no freedom and cannot move; its slot, and its
    // place in the pose arrays, stay.
    assert_eq!(dof(&s, &solved), 0);
    let movable = assembly_mobility(&s, &solved.positions, &solved.rotations);
    assert_eq!(movable[SLIDER], [false; 6]);
    assert_eq!(s.validate(), Ok(()));
    let result = resolve(&s, &solved);
    assert!(result.success);
    assert_eq!(result.positions.len(), 3);
    assert_eq!(result.positions[SLIDER], solved.positions[SLIDER]);
    assert!(s.remove_component(SLIDER).is_empty(), "already removed");

    // Components added afterwards get new ids and work as usual.
    let part = s.add_component();
    assert_eq!(part, 3);
    s.add_constraint(Constraint::lock(
        7,
        GROUND,
        part,
        DVec3::new(1.0, 2.0, 3.0),
        DQuat::IDENTITY,
    ));
    let mut positions = result.positions.clone();
    let mut rotations = result.rotations.clone();
    positions.push(DVec3::ZERO);
    rotations.push(DQuat::IDENTITY);
    let result = solver().solve_from(&s, positions, rotations);
    assert!(result.success);
    assert!((result.positions[part] - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-9);
    assert_eq!(dof(&s, &result), 0);
}

#[test]
fn structural_edits_bump_the_topology_revision() {
    let mut s = ConstraintSystem::new();
    let mut last = s.topology_revision();
    let mut bumped = |s: &ConstraintSystem, what: &str| {
        let now = s.topology_revision();
        assert!(now > last, "{what} must change the revision");
        last = now;
    };

    let a = s.add_component();
    bumped(&s, "adding a component");
    let b = s.add_component();
    let c = s.add_constraint(Constraint::distance(0, a, DVec3::ZERO, b, DVec3::ZERO, 1.0));
    bumped(&s, "adding a constraint");
    s.set_constraint_enabled(c, false);
    bumped(&s, "suppressing");
    s.set_constraint_enabled(c, true);
    bumped(&s, "enabling");
    s.replace_constraint(c, Constraint::coincident(0, a, DVec3::ZERO, b, DVec3::ZERO));
    bumped(&s, "replacing");
    s.remove_constraint(c);
    bumped(&s, "removing a constraint");
    s.remove_component(b);
    bumped(&s, "removing a component");

    // Enabling what is already enabled is not a change.
    let c2 = s.add_constraint(Constraint::fixed(1, a, DVec3::ZERO));
    let before = s.topology_revision();
    s.set_constraint_enabled(c2, true);
    assert_eq!(s.topology_revision(), before);
}

#[test]
fn validate_reports_what_is_wrong_with_the_model() {
    let (mut s, _, _) = assembly(0.5, 2.0);
    assert_eq!(s.validate(), Ok(()));

    // A component that was never added.
    let missing = s.add_constraint(Constraint::distance(
        10,
        GROUND,
        DVec3::ZERO,
        7,
        DVec3::ZERO,
        1.0,
    ));
    // A value that is not a number.
    let nan = s.add_constraint(Constraint::distance(
        11,
        GROUND,
        DVec3::ZERO,
        DOOR,
        DVec3::ZERO,
        f64::NAN,
    ));
    // A direction that points nowhere.
    let no_direction = s.add_constraint(Constraint::parallel(
        12,
        GROUND,
        DVec3::Z,
        DOOR,
        DVec3::ZERO,
    ));
    // The same part on both sides.
    let both_sides = s.add_constraint(Constraint::coincident(13, DOOR, DVec3::X, DOOR, DVec3::Y));

    let errors = s.validate().expect_err("four defects");
    assert_eq!(
        errors,
        vec![
            AssemblyError::ComponentOutOfRange {
                constraint: missing,
                component: 7,
            },
            AssemblyError::NonFiniteValue { constraint: nan },
            AssemblyError::ZeroLength {
                constraint: no_direction,
            },
            AssemblyError::SameComponent {
                constraint: both_sides,
                component: DOOR,
            },
        ]
    );
    assert!(errors[0].to_string().contains("component 7"));

    // Removing the offenders clears the report.
    for index in [missing, nan, no_direction, both_sides] {
        s.remove_constraint(index);
    }
    assert_eq!(s.validate(), Ok(()));

    // A constraint added against a component that is already gone.
    s.remove_component(SLIDER);
    let stale = s.add_constraint(Constraint::coincident(
        14,
        GROUND,
        DVec3::ZERO,
        SLIDER,
        DVec3::ZERO,
    ));
    assert_eq!(
        s.validate(),
        Err(vec![AssemblyError::ComponentRemoved {
            constraint: stale,
            component: SLIDER,
        }])
    );
}

/// The deprecated analyzers see the edited assembly too.
#[test]
#[allow(deprecated)]
fn the_old_analyzers_follow_the_edits() {
    use rcad_gcs::rigid3d::{DOFAnalyzer, SolverDiagnostics};

    let (mut s, angle_dim, _) = assembly(0.5, 2.0);
    s.remove_component(SLIDER);
    s.set_constraint_enabled(angle_dim, false);

    let mut diag = SolverDiagnostics::new();
    diag.analyze(&s);
    assert_eq!(diag.stats.total_components, 2);
    assert_eq!(diag.stats.total_constraints, 3);
    assert_eq!(diag.stats.total_dof, 12);
    assert_eq!(diag.stats.dof_remaining, 1, "the door swings");

    let dof = DOFAnalyzer::new().analyze(&s);
    assert_eq!(dof.total_dof, 1);
    assert_eq!(dof.component_dof, vec![0, 1, 0]);
}

#[cfg(feature = "serde")]
#[test]
fn the_edit_state_round_trips_through_json() {
    let (mut s, angle_dim, distance_dim) = assembly(0.5, 2.0);
    s.set_constraint_enabled(angle_dim, false);
    s.remove_constraint(distance_dim);

    let json = serde_json::to_string(&s).expect("serialize");
    let back: ConstraintSystem = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(back, s);
    assert_eq!(
        back.constraint_state(angle_dim),
        ConstraintState::Suppressed
    );
    assert_eq!(
        back.constraint_state(distance_dim),
        ConstraintState::Removed
    );
    assert_eq!(back.topology_revision(), s.topology_revision());

    // A system written without any edit state reads as all active.
    let plain = r#"{
        "constraints": [
            {"id": 5, "kind": {"Fixed": {"body": 0, "position": [1.0, 2.0, 3.0], "orientation": null}}}
        ],
        "num_components": 1
    }"#;
    let system: ConstraintSystem = serde_json::from_str(plain).expect("deserialize");
    assert!(system.is_active(0));
    assert!(!system.is_component_removed(0));
    let result = ConstraintSolver::new().solve(&system);
    assert!(result.success);
    assert_eq!(result.positions[0], DVec3::new(1.0, 2.0, 3.0));
}
