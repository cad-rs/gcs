//! Which branch a configuration is on: the signs a continuous motion
//! keeps.
//!
//! The solution set of a constraint system is rarely one connected piece.
//! A tangency can be met on either side of its line, an angle can open
//! either way, a rectangle held by horizontal and vertical sides has a
//! mirrored twin on the far side of zero width. Which of those the
//! geometry is on cannot be read from the residuals — they are zero on
//! all of them — so it is read from the geometry itself: each constraint
//! that has a side to be on contributes a signed, scale-free quantity
//! (`Equations::branch_witnesses`), and the signs of those quantities are
//! the signature of the configuration.
//!
//! A drag takes the signature of the configuration it starts from and
//! holds the geometry to it, shortening its step rather than crossing
//! over; a solve compares the signature before and after and says, in
//! [`SolveStatus::BranchFlip`](crate::SolveStatus::BranchFlip), when the
//! answer came out on another branch.
//!
//! What a signature does **not** cover: two solutions that are not
//! separated by a sign of this kind — a four-bar linkage's elbow, which
//! turns over as the linkage straightens, is a mirror image with no
//! constraint of its own to witness it. Measured on this crate's scripted
//! drags, such a crossing is a perfectly continuous motion (the geometry
//! moves no further in that frame than in its neighbours, 0.38 against a
//! median of 2.0 times the cursor's own travel), so no continuity test
//! separates it either. What keeps a drag off it is not crossing into the
//! configuration where the linkage is straight in the first place — the
//! business of walking only as far as the feasible set reaches.

use crate::numerics::SolveStatus;
use crate::system::Equations;

/// How small a witness counts as zero: the quantities are scale-free, so
/// this is an absolute number. A configuration exactly on the boundary —
/// a tangency through its line, an angle at zero — has no side, and the
/// signature neither holds it to one nor reads a crossing from it: the
/// sign on the far side is compared with the last one that was not zero.
const FLAT: f64 = 1.0e-9;

/// The signs of the branch witnesses of the constraints `indices`, in
/// order.
pub(crate) fn signature<E: Equations + ?Sized>(
    eqs: &E,
    indices: &[usize],
    params: &[f64],
) -> Vec<i8> {
    let mut signs = Vec::new();
    let mut witnesses: Vec<f64> = Vec::new();
    for &index in indices {
        witnesses.clear();
        eqs.branch_witnesses(index, params, &mut witnesses);
        signs.extend(witnesses.iter().map(|&w| {
            if w > FLAT {
                1
            } else if w < -FLAT {
                -1
            } else {
                0
            }
        }));
    }
    signs
}

/// The signature of the whole system: what a solve takes before and after
/// itself.
pub(crate) fn of_system<E: Equations + ?Sized>(eqs: &E, params: &[f64]) -> Vec<i8> {
    signature(eqs, &crate::system::all_indices(eqs), params)
}

/// Raise `status` to [`SolveStatus::BranchFlip`] when the system has
/// turned over since `before` was taken: what every solve does with the
/// signature it kept (the solution it found is the one it keeps; putting
/// it back is a drag's business).
pub(crate) fn flag_turn_over<E: Equations + ?Sized>(
    eqs: &E,
    params: &[f64],
    before: &[i8],
    status: &mut SolveStatus,
) {
    if turned_over(before, &of_system(eqs, params)) {
        *status = (*status).max(SolveStatus::BranchFlip);
    }
}

/// Whether `now` has turned any of `was` over: a sign that was one thing
/// and is now the other. A witness that is flat at either end has no side
/// to have changed.
pub(crate) fn turned_over(was: &[i8], now: &[i8]) -> bool {
    was.iter()
        .zip(now)
        .any(|(was, is)| *was != 0 && *is != 0 && was != is)
}

/// Carry forward the signs of `now` that are not flat, leaving the last
/// known side in place for the ones that are: a configuration passing
/// through a boundary keeps the side it came from until it has a side
/// again.
pub(crate) fn carry_forward(signature: &mut [i8], now: &[i8]) {
    for (slot, &sign) in signature.iter_mut().zip(now) {
        if sign != 0 {
            *slot = sign;
        }
    }
}
