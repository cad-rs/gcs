//! Dimension jumps: reaching a value a solve cannot go straight to.
//!
//! A host types a new number into a dimension and asks for a solve.
//! Newton starts from the geometry as it stands, and that is usually
//! enough: a width of 40 instead of 42 is two Newton steps away. A jump
//! of a different size is not. Driving the same width from 8 to 80 can
//! leave Newton in a place where it stops gaining (`Stalled`), and
//! driving it through zero takes the geometry to the mirror image of
//! itself — a solution of the constraints, but not the configuration the
//! host was looking at ([`crate::branch`]).
//!
//! What a hand would do in a sketcher is drag the dimension there, a bit
//! at a time, and that is what this module does when a solve cannot do it
//! in one go: it walks the value from where the last solve met it to
//! where the host put it, solving at each stop, halving a step that fails
//! and giving up below [`MIN_STEP`] of the whole jump. A step that lands
//! on a configuration with no side to it — a rectangle of exactly zero
//! width, a linkage exactly straight — is a legitimate place to stand:
//! the walk goes on from there, and a value the host asked for on the
//! other side of it is reached continuously rather than jumped to.
//!
//! The walk only ever happens after a solve that did not arrive, so a
//! sketch that solves is solved exactly as it was before this module
//! existed, down to the last bit and the iteration count.

use crate::numerics::{SolveOptions, SolveResult, SolveStatus};

/// The smallest fraction of a jump a walk will try before it gives up.
/// A step this short that still fails says the way is blocked, not that
/// the step is too long: every jump in the 2D corpus that can be walked
/// at all is walked in halves (measured, see `docs/CATCHUP_PLAN.md`
/// §13.7).
pub(crate) const MIN_STEP: f64 = 1.0 / 32.0;

/// What a family offers so that its jumps can be walked. The one
/// implementor per family is its sketch type.
pub(crate) trait Jumping {
    /// The parameters a solve moves.
    fn jump_params(&mut self) -> &mut Vec<f64>;

    /// Put `value` in dimension slot `slot`. This is the walk setting its
    /// own waypoints, not a host changing a dimension: it does not count
    /// as a jump of its own.
    fn jump_set_value(&mut self, slot: usize, value: f64);

    /// One solve of the whole sketch from the geometry as it stands, with
    /// no walking of its own.
    fn jump_solve(&mut self, options: &SolveOptions) -> SolveResult;

    /// The residual of the geometry as it stands against the values as
    /// they stand.
    fn jump_residual(&self) -> f64;
}

/// Whether a solve arrived where it was asked to go: the constraints are
/// met, in the configuration they were met in before.
pub(crate) fn arrived(result: &SolveResult) -> bool {
    result.converged && result.status != SolveStatus::BranchFlip
}

/// Walk `jump` — one `(slot, the value the last solve met, the value the
/// host asked for)` per dimension the host changed — from where it was to
/// where it is asked to be.
///
/// `before` is the geometry the failed solve started from, and `direct`
/// its result. The walk starts at half the jump, since the whole of it is
/// what `direct` just tried.
///
/// Arriving returns the last step's result, with the iterations of the
/// whole walk.
///
/// Giving up leaves the values the host asked for — a library does not
/// put a dimension back — and then depends on what the direct solve had
/// in hand. A solve that did not converge had nothing: the geometry stays
/// at the furthest configuration the walk reached and the status is the
/// direct one, so what the host hears is that their value was not
/// reached and what they see is how far it got. A solve that came out on
/// another branch had a solution of the constraints, just not the one
/// that was asked for; that is worth more than a geometry which is
/// halfway and satisfies nothing, so it is given back as it was, flip and
/// all. (The internal tangency of two spheres is the case in hand: the
/// radii can only pass each other through the configuration where the
/// centres coincide, where Newton has nowhere to step, so the swap can
/// be solved for but not walked to.)
pub(crate) fn walk<T: Jumping + ?Sized>(
    sketch: &mut T,
    jump: &[(usize, f64, f64)],
    before: &[f64],
    direct: SolveResult,
    options: &SolveOptions,
) -> SolveResult {
    // What the direct solve left behind, in case the walk does worse.
    let found: Vec<f64> = sketch.jump_params().clone();
    sketch.jump_params().copy_from_slice(before);
    let mut best: Vec<f64> = before.to_vec();
    let mut iterations = direct.iterations;
    let mut reached = 0.0_f64;
    let mut step = 0.5_f64;
    let mut last: Option<SolveResult> = None;

    while reached < 1.0 {
        let stop = (reached + step).min(1.0);
        for &(slot, was, asked) in jump {
            sketch.jump_set_value(slot, was + (asked - was) * stop);
        }
        let result = sketch.jump_solve(options);
        iterations += result.iterations;
        if arrived(&result) {
            reached = stop;
            best.copy_from_slice(sketch.jump_params());
            last = Some(result);
        } else {
            // This step never happened. Try half of it.
            sketch.jump_params().copy_from_slice(&best);
            step *= 0.5;
            if step < MIN_STEP {
                break;
            }
        }
    }

    // The values are the host's either way.
    for &(slot, _, asked) in jump {
        sketch.jump_set_value(slot, asked);
    }

    match last.filter(|_| reached >= 1.0) {
        Some(mut result) => {
            result.iterations = iterations;
            result
        }
        None if direct.converged => {
            // A branch flip is still a solution: back to the one the
            // direct solve found, said of in the same words as before.
            sketch.jump_params().copy_from_slice(&found);
            SolveResult {
                iterations,
                ..direct
            }
        }
        None => {
            sketch.jump_params().copy_from_slice(&best);
            SolveResult {
                converged: false,
                residual: sketch.jump_residual(),
                iterations,
                status: direct.status,
                satisfied_level: 0,
            }
        }
    }
}
