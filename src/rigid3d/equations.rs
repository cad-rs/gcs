//! The rigid-pose family as the shared pipeline's [`Equations`].
//!
//! A body is stored as seven numbers — translation `(x, y, z)` and
//! orientation quaternion `(x, y, z, w)` — and has six columns: translation
//! and a body-frame rotation vector. A step in those columns is applied by
//! [`PoseSpace`]: the translation adds, the orientation is rotated by the
//! step and renormalised. Counting in columns is what makes the degrees of
//! freedom come out right without any gauge bookkeeping.

use std::cell::RefCell;

use glam::{DQuat, DVec3};

use super::constraint::{Constraint, ConstraintKind, ConstraintSystem};
use super::solver::{ConstraintSolver, ResidualJacobian};
use crate::numerics::StepSpace;
use crate::system::{Equations, SparseRow};

/// Stored numbers per body: translation (3) + quaternion `(x, y, z, w)` (4).
pub(crate) const STORED_PER_BODY: usize = 7;
/// Columns per body: translation (3) + body-frame rotation vector (3).
pub(crate) const COLUMNS_PER_BODY: usize = 6;

/// Flatten poses into the pipeline's parameter vector.
pub(crate) fn pack(positions: &[DVec3], rots: &[DQuat]) -> Vec<f64> {
    let mut params = Vec::with_capacity(positions.len() * STORED_PER_BODY);
    for (t, q) in positions.iter().zip(rots) {
        params.extend_from_slice(&[t.x, t.y, t.z, q.x, q.y, q.z, q.w]);
    }
    params
}

/// The poses held in a parameter vector written by [`pack`].
pub(crate) fn unpack(params: &[f64]) -> (Vec<DVec3>, Vec<DQuat>) {
    let bodies = params.len() / STORED_PER_BODY;
    (
        (0..bodies).map(|b| position(params, b)).collect(),
        (0..bodies).map(|b| rotation(params, b)).collect(),
    )
}

fn position(params: &[f64], body: usize) -> DVec3 {
    let s = &params[body * STORED_PER_BODY..];
    DVec3::new(s[0], s[1], s[2])
}

fn rotation(params: &[f64], body: usize) -> DQuat {
    let s = &params[body * STORED_PER_BODY..];
    DQuat::from_xyzw(s[3], s[4], s[5], s[6])
}

/// An assembly's constraints as the pipeline's [`Equations`].
pub(crate) struct RigidEquations<'a> {
    system: &'a ConstraintSystem,
    /// Per body: grounded, so not an unknown (see [`RigidEquations::grounded`]).
    held: Vec<bool>,
    /// Per constraint: the `Fixed` that defines a grounded body.
    definition: Vec<bool>,
    /// Per constraint: the constraint as a system of its own, over just its
    /// own bodies.
    alone: Vec<Alone>,
    /// Per constraint, its last evaluation. The Newton loop asks for the
    /// residuals and then the Jacobian at the same poses, and a constraint
    /// evaluation yields both, so the second request is served from here.
    /// The buffers are reused.
    last: RefCell<Vec<Evaluation>>,
}

/// One constraint as a system over its own one or two bodies.
struct Alone {
    /// The constraint, its component indices renumbered to `0` / `1`.
    system: ConstraintSystem,
    /// The constraint's bodies in the full system, ascending.
    bodies: [usize; 2],
    /// How many bodies there are (one or two).
    count: usize,
}

/// One constraint evaluated at particular poses of its own bodies.
struct Evaluation {
    /// Whether anything has been evaluated yet.
    filled: bool,
    /// The stored parameters of the constraint's bodies it was evaluated at.
    inputs: [f64; 2 * STORED_PER_BODY],
    /// Residuals, and the Jacobian over the columns of the constraint's
    /// bodies.
    rj: ResidualJacobian,
}

impl<'a> RigidEquations<'a> {
    /// Every body is an unknown and every constraint an equation — the
    /// view the diagnosis takes.
    pub(crate) fn new(system: &'a ConstraintSystem) -> Self {
        let alone = system
            .constraints
            .iter()
            .map(|constraint| {
                let (first, second) = constraint.kind.bodies();
                let (bodies, count) = match second {
                    Some(second) if second != first => ([first.min(second), first.max(second)], 2),
                    _ => ([first, first], 1),
                };
                let local_index = |body: usize| usize::from(count == 2 && body == bodies[1]);
                let local = Constraint::new(constraint.id, constraint.kind.map_bodies(local_index));
                let mut alone = ConstraintSystem::with_components(count);
                alone.add_constraint(local);
                Alone {
                    system: alone,
                    bodies,
                    count,
                }
            })
            .collect();
        let last = (0..system.constraints.len())
            .map(|_| Evaluation {
                filled: false,
                inputs: [0.0; 2 * STORED_PER_BODY],
                rj: ResidualJacobian {
                    r: Vec::new(),
                    j: Vec::new(),
                },
            })
            .collect();
        Self {
            system,
            held: vec![false; system.num_components],
            definition: vec![false; system.constraints.len()],
            alone,
            last: RefCell::new(last),
        }
    }

    /// The view the solve takes: a body that a `Fixed` constraint pins in
    /// full (position and orientation) is grounded. Its pose is known, so
    /// it is not an unknown: its columns are left out, and the constraint
    /// that defines it contributes no equation (it holds by construction
    /// once [`RigidEquations::snap_grounded`] has placed the body). Mates
    /// to a grounded body then couple nothing — the clusters hanging off
    /// the ground separate.
    ///
    /// The first full `Fixed` of a body is its definition. A second one, or
    /// any other constraint between grounded bodies only, stays an
    /// equation with nothing left to move: it holds or the solve fails.
    pub(crate) fn grounded(system: &'a ConstraintSystem) -> Self {
        let mut equations = Self::new(system);
        for (index, constraint) in system.active_constraints() {
            if let ConstraintKind::Fixed {
                body,
                orientation: Some(_),
                ..
            } = constraint.kind
            {
                if !equations.held[body] {
                    equations.held[body] = true;
                    equations.definition[index] = true;
                }
            }
        }
        equations
    }

    /// Place the grounded bodies that the constraints `group` define at
    /// their fixed pose.
    pub(crate) fn snap_grounded(&self, group: &[usize], params: &mut [f64]) {
        for &index in group {
            if !self.definition[index] {
                continue;
            }
            let ConstraintKind::Fixed {
                body,
                position: t,
                orientation: Some(orientation),
            } = self.system.constraints[index].kind
            else {
                continue;
            };
            let q = orientation.normalize();
            let base = body * STORED_PER_BODY;
            params[base..base + STORED_PER_BODY]
                .copy_from_slice(&[t.x, t.y, t.z, q.x, q.y, q.z, q.w]);
        }
    }

    /// L2 norm of the residuals of every constraint at `params`, the
    /// grounding definitions included: the error of the system as authored,
    /// whichever view this is.
    pub(crate) fn residual_norm_of_every_constraint(&self, params: &[f64]) -> f64 {
        let mut sum = 0.0_f64;
        for (index, _) in self.system.active_constraints() {
            self.with_evaluation(index, params, |rj, _| {
                for r in &rj.r {
                    sum += r * r;
                }
            });
        }
        sum.sqrt()
    }

    /// Run `read` on the evaluation of constraint `index` at `params`:
    /// residuals and Jacobian over the columns of its own bodies only.
    ///
    /// The numbers come from the same builder the whole-system assembly
    /// uses, run on the constraint alone: a constraint's rows depend on no
    /// other body.
    fn with_evaluation<R>(
        &self,
        index: usize,
        params: &[f64],
        read: impl FnOnce(&ResidualJacobian, &Alone) -> R,
    ) -> R {
        let alone = &self.alone[index];
        let mut inputs = [0.0_f64; 2 * STORED_PER_BODY];
        for (k, &body) in alone.bodies[..alone.count].iter().enumerate() {
            inputs[k * STORED_PER_BODY..(k + 1) * STORED_PER_BODY]
                .copy_from_slice(&params[body * STORED_PER_BODY..(body + 1) * STORED_PER_BODY]);
        }

        let mut last = self.last.borrow_mut();
        let seen = &mut last[index];
        // Bitwise comparison: the same poses give the same numbers.
        if !seen.filled || seen.inputs.map(f64::to_bits) != inputs.map(f64::to_bits) {
            let positions = [
                position(params, alone.bodies[0]),
                position(params, alone.bodies[1]),
            ];
            let rots = [
                rotation(params, alone.bodies[0]),
                rotation(params, alone.bodies[1]),
            ];
            ConstraintSolver::fill_residual_jacobian(
                &alone.system,
                &positions[..alone.count],
                &rots[..alone.count],
                &mut seen.rj,
            );
            seen.inputs = inputs;
            seen.filled = true;
        }
        read(&seen.rj, alone)
    }
}

impl Equations for RigidEquations<'_> {
    fn constraint_count(&self) -> usize {
        self.system.constraints.len()
    }

    /// A suppressed or removed constraint has no equation.
    fn equation_count(&self, index: usize) -> usize {
        if self.system.is_active(index) {
            ConstraintSolver::equation_count_of(&self.system.constraints[index])
        } else {
            0
        }
    }

    fn residuals(&self, index: usize, params: &[f64], out: &mut [f64]) {
        if self.definition[index] || !self.system.is_active(index) {
            out.fill(0.0);
            return;
        }
        self.with_evaluation(index, params, |rj, _| out.copy_from_slice(&rj.r));
    }

    fn jacobian(&self, index: usize, params: &[f64], out: &mut [SparseRow]) {
        if self.definition[index] || !self.system.is_active(index) {
            return;
        }
        self.with_evaluation(index, params, |rj, alone| {
            for (row, entries) in rj.j.iter().zip(out.iter_mut()) {
                for (k, &body) in alone.bodies[..alone.count].iter().enumerate() {
                    if self.held[body] {
                        continue;
                    }
                    for d in 0..COLUMNS_PER_BODY {
                        entries.push((body * COLUMNS_PER_BODY + d, row[k * COLUMNS_PER_BODY + d]));
                    }
                }
            }
        });
    }

    fn column_count(&self, _params: &[f64]) -> usize {
        self.system.num_components * COLUMNS_PER_BODY
    }

    fn step_space(&self, free: &[usize]) -> Box<dyn StepSpace + '_> {
        Box::new(PoseSpace::new(free))
    }

    /// The family has its own analytic-vs-finite-difference gate test, and
    /// the per-solve check would cost two whole-cluster residual
    /// evaluations per column here.
    fn debug_check_jacobian(&self) -> bool {
        false
    }
}

/// The step space of rigid poses: translations add, orientations rotate.
pub(crate) struct PoseSpace {
    /// Bodies with at least one free column, ascending.
    bodies: Vec<usize>,
    /// For each of those bodies, where each of its six columns sits in the
    /// step vector (`None`: the column is held).
    slots: Vec<[Option<usize>; COLUMNS_PER_BODY]>,
    /// The stored parameters of those bodies.
    storage: Vec<usize>,
    dim: usize,
}

impl PoseSpace {
    /// `free` lists the free columns in ascending order.
    pub(crate) fn new(free: &[usize]) -> Self {
        let mut bodies: Vec<usize> = Vec::new();
        let mut slots: Vec<[Option<usize>; COLUMNS_PER_BODY]> = Vec::new();
        for (k, &col) in free.iter().enumerate() {
            let body = col / COLUMNS_PER_BODY;
            if bodies.last() != Some(&body) {
                bodies.push(body);
                slots.push([None; COLUMNS_PER_BODY]);
            }
            if let Some(last) = slots.last_mut() {
                last[col % COLUMNS_PER_BODY] = Some(k);
            }
        }
        let storage = bodies
            .iter()
            .flat_map(|&b| b * STORED_PER_BODY..(b + 1) * STORED_PER_BODY)
            .collect();
        Self {
            bodies,
            slots,
            storage,
            dim: free.len(),
        }
    }
}

impl StepSpace for PoseSpace {
    fn dim(&self) -> usize {
        self.dim
    }

    fn storage(&self) -> &[usize] {
        &self.storage
    }

    fn apply(&self, params: &mut [f64], saved: &[f64], step: &[f64]) {
        for (n, (&body, slots)) in self.bodies.iter().zip(&self.slots).enumerate() {
            let from = &saved[n * STORED_PER_BODY..(n + 1) * STORED_PER_BODY];
            let base = body * STORED_PER_BODY;

            // Translation: add. A held column keeps its value exactly.
            for d in 0..3 {
                params[base + d] = match slots[d] {
                    Some(k) => from[d] + step[k],
                    None => from[d],
                };
            }

            // Orientation: rotate by the body-frame rotation vector and
            // renormalise. No rotation leaves the quaternion bit for bit.
            let at = |d: usize| slots[d].map_or(0.0, |k| step[k]);
            let turn = DVec3::new(at(3), at(4), at(5));
            let start = DQuat::from_xyzw(from[3], from[4], from[5], from[6]);
            let q = if turn == DVec3::ZERO {
                start
            } else {
                (start * DQuat::from_scaled_axis(turn)).normalize()
            };
            params[base + 3] = q.x;
            params[base + 4] = q.y;
            params[base + 5] = q.z;
            params[base + 6] = q.w;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pack_and_unpack_round_trip_bitwise() {
        let positions = vec![DVec3::new(1.0, -2.5, 0.125), DVec3::new(0.0, 7.0, -3.0)];
        let rots = vec![
            DQuat::from_axis_angle(DVec3::new(1.0, 2.0, 3.0).normalize(), 0.7),
            DQuat::IDENTITY,
        ];
        let (p, r) = unpack(&pack(&positions, &rots));
        assert_eq!(p, positions);
        assert_eq!(r, rots);
    }

    #[test]
    fn a_step_translates_and_rotates_one_body() {
        let positions = vec![DVec3::ZERO, DVec3::new(1.0, 2.0, 3.0)];
        let rots = vec![DQuat::IDENTITY, DQuat::IDENTITY];
        let mut params = pack(&positions, &rots);

        // Body 1 only: all six of its columns free.
        let free: Vec<usize> = (6..12).collect();
        let space = PoseSpace::new(&free);
        assert_eq!(space.dim(), 6);
        assert_eq!(space.storage(), (7..14).collect::<Vec<_>>().as_slice());

        let saved: Vec<f64> = space.storage().iter().map(|&i| params[i]).collect();
        let quarter_turn = std::f64::consts::FRAC_PI_2;
        space.apply(
            &mut params,
            &saved,
            &[0.5, 0.0, -1.0, 0.0, 0.0, quarter_turn],
        );

        let (p, r) = unpack(&params);
        assert_eq!(p[0], DVec3::ZERO, "body 0 is not part of the space");
        assert_eq!(r[0], DQuat::IDENTITY);
        assert_eq!(p[1], DVec3::new(1.5, 2.0, 2.0));
        let x_axis = r[1].mul_vec3(DVec3::X);
        assert!(
            (x_axis - DVec3::Y).length() < 1e-12,
            "a quarter turn about z"
        );
        assert!((r[1].length() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn a_zero_step_leaves_the_pose_bit_for_bit() {
        // A quaternion that is unit only to rounding: renormalising it
        // would move the last bits.
        let q = DQuat::from_xyzw(0.1, 0.2, 0.3, 0.927_361_849_549_570_4);
        let mut params = pack(&[DVec3::new(-0.0, 1.0, 2.0)], &[q]);
        let before = params.clone();
        let space = PoseSpace::new(&[0, 1, 2, 3, 4, 5]);
        let saved: Vec<f64> = space.storage().iter().map(|&i| params[i]).collect();

        // Held columns and a zero rotation: only the y translation moves.
        let space_y = PoseSpace::new(&[1]);
        let saved_y: Vec<f64> = space_y.storage().iter().map(|&i| params[i]).collect();
        space_y.apply(&mut params, &saved_y, &[0.0]);
        assert_eq!(
            params.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            before.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );

        // And the full space reproduces the saved state when asked to.
        space.apply(&mut params, &saved, &[1.0, 1.0, 1.0, 0.3, 0.0, 0.0]);
        for (&i, &v) in space.storage().iter().zip(&saved) {
            params[i] = v;
        }
        assert_eq!(params, before);
    }
}
