//! Drag sessions: the cursor as a soft target, the constraints kept hard.
//!
//! A drag is not a constraint. Pinning the dragged point where the cursor
//! is adds two equations, which a sketch with fewer free motions than that
//! cannot take: a one-degree-of-freedom four-bar linkage, whose coupler
//! point may only travel along one curve, is over-constrained by a cursor
//! that is anywhere else, and the solve fails on nearly every frame
//! (measured: 95 of 97). Writing the cursor into the point's parameters
//! instead and re-solving keeps the constraints satisfiable, but nothing
//! then holds the point near the cursor — the solve is free to put the
//! whole sketch back where it was, and a closed drag loop does not come
//! back to the configuration it started from.
//!
//! So the two are kept apart, in layers. Every frame:
//!
//! 1. **Predictor.** Over the motions the constraints leave free at the
//!    current configuration — the null space of the Jacobian — take the
//!    shortest motion that moves the dragged points as near their cursors
//!    as those motions allow (`system::tangent_step`). Nothing
//!    else moves any further than that asks.
//! 2. **Corrector.** Re-solve the components the drag touches, warm, from
//!    there: the linearisation of the predictor leaves a second-order
//!    residual, and the solve takes the hard constraints back to zero.
//!
//! A frame walks that pair until the cursor is reached, the walk stops
//! gaining, or the corrector declines to converge — the last configuration
//! that did converge is what the sketch keeps. A target outside what the
//! constraints can reach therefore leaves the geometry at the nearest
//! configuration the walk found, with [`DragFrame::reached`] false, rather
//! than deforming the sketch toward it.

use crate::entity::PointRef;
use crate::numerics::{CoreOptions, SolveStatus};
use crate::system::{self, Equations, Plans, Pull, SparseRow};

/// How many predictor–corrector pairs one frame may walk. A cursor that
/// moves by a frame's worth of pixels is one step; the cap is for a frame
/// that jumps (a tablet, a dropped frame, a host that drags by keyboard).
const MAX_SUBSTEPS: usize = 8;

/// How far one substep may ask a dragged point to travel, as a fraction of
/// the length scale of what it is attached to. The predictor is a
/// linearisation: asked for more than this it would aim along a tangent
/// that has already turned, and the corrector would spend the gain putting
/// the constraints back. The walk covers a longer jump in several steps,
/// re-linearising at each.
const ADVANCE: f64 = 0.25;

/// When the dragged points are this near their cursors, relative to the
/// length scale, the frame has arrived: the predictor's own accuracy is
/// about the conditioning of the regularised projector, so asking for
/// more would spend substeps without moving anything.
const TRACK_TOL: f64 = 1.0e-9;

/// A substep has to cut the tracking error by this fraction of what was
/// left, or the walk stops: either the cursor is out of reach, or the
/// linearisation has nothing more to give at this configuration.
const GAIN_MIN: f64 = 1.0e-3;

/// The token of an open drag session: what [`crate::Sketch::begin_drag`]
/// returns and the other two take.
///
/// It names the points being dragged and which session it belongs to; the
/// state of the drag (the configuration to fall back to, and what the
/// session has to keep between frames) stays in the sketch, so holding one
/// of these borrows nothing and a drag may be diagnosed while it runs.
///
/// Beginning a drag ends any drag already running, which leaves the older
/// token stale: a frame asked for with it moves nothing and says
/// [`SolveStatus::Stalled`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Drag {
    /// Which session: compared against the sketch's own.
    pub(crate) id: u64,
    /// The dragged points, in the order their cursors are given.
    pub(crate) handles: Vec<PointRef>,
}

impl Drag {
    /// The points this drag moves, in the order their cursors are given.
    pub fn handles(&self) -> &[PointRef] {
        &self.handles
    }
}

/// What one frame of a drag did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DragFrame {
    /// How the corrector of the last substep ended, or
    /// [`SolveStatus::Converged`] for a frame whose cursors were already
    /// reached and which therefore solved nothing.
    pub status: SolveStatus,
    /// The RMS residual of the constraints of the components the drag
    /// touches, as the frame leaves them. The rest of the sketch is not
    /// looked at: a drag does not move it.
    pub residual: f64,
    /// Newton iterations the frame's correctors took, all substeps
    /// together.
    pub iterations: usize,
    /// How far the dragged points ended up from their cursors: the
    /// largest of the distances, zero when every one arrived.
    pub tracking_error: f64,
    /// Whether the cursors were reached. False means the constraints do
    /// not allow it and the geometry is at the nearest configuration the
    /// frame found — not that the solve failed; `status` says that.
    pub reached: bool,
}

/// What a sketch keeps while a drag runs.
#[derive(Debug, Clone)]
pub(crate) struct Session {
    /// Which session this is; an older token does not match.
    pub(crate) id: u64,
    /// The parameters the cursors pull, as many per handle as a point has
    /// coordinates, in the order of the handles.
    pub(crate) pulled: Vec<usize>,
    /// A gradient row per pulled parameter (the parameter itself), kept so
    /// that every substep hands the same rows to the predictor.
    rows: Vec<SparseRow>,
    /// The components the pulled parameters belong to, as ordinals into
    /// the plans' components, ascending.
    components: Vec<usize>,
    /// The topology revision `components` was worked out at.
    revision: u64,
}

impl Session {
    /// A session dragging `pulled`, over the components of `plans`.
    pub(crate) fn new<E: Equations + ?Sized>(
        id: u64,
        pulled: Vec<usize>,
        revision: u64,
        eqs: &E,
        params: &[f64],
        plans: &Plans,
    ) -> Self {
        let rows = pulled.iter().map(|&p| vec![(p, 1.0)]).collect();
        let components = components_of(eqs, params, plans, &pulled);
        Self {
            id,
            pulled,
            rows,
            components,
            revision,
        }
    }

    /// Take the components afresh if the structure has changed under the
    /// drag (a host may add or remove a constraint between frames).
    fn refresh<E: Equations + ?Sized>(
        &mut self,
        revision: u64,
        eqs: &E,
        params: &[f64],
        plans: &Plans,
    ) {
        if self.revision != revision {
            self.components = components_of(eqs, params, plans, &self.pulled);
            self.revision = revision;
        }
    }
}

/// The components of `plans` that any of `pulled` belongs to, as ordinals,
/// ascending.
fn components_of<E: Equations + ?Sized>(
    eqs: &E,
    params: &[f64],
    plans: &Plans,
    pulled: &[usize],
) -> Vec<usize> {
    plans
        .components()
        .iter()
        .enumerate()
        .filter(|(_, group)| {
            let touched = system::touched_params(eqs, group, params);
            pulled.iter().any(|p| touched.contains(p))
        })
        .map(|(ordinal, _)| ordinal)
        .collect()
}

/// One frame of a drag: walk the predictor–corrector pair until the
/// `cursors` (one per pulled parameter, in the order of
/// [`Session::pulled`]) are reached or the walk stops gaining.
///
/// `distance` turns the per-parameter errors into the one tracking error
/// the frame reports — for a family whose points have two coordinates, the
/// largest of the planar distances. `weight_of` says how reluctant each
/// parameter is to move (see [`system::tangent_step`]).
#[allow(clippy::too_many_arguments)]
pub(crate) fn frame<E: Equations + ?Sized>(
    eqs: &E,
    params: &mut [f64],
    fixed: &[bool],
    plans: &Plans,
    session: &mut Session,
    revision: u64,
    cursors: &[f64],
    distance: &dyn Fn(&[f64]) -> f64,
    weight_of: &dyn Fn(usize) -> f64,
    options: &CoreOptions,
) -> DragFrame {
    session.refresh(revision, eqs, params, plans);
    // A pulled parameter no constraint touches is free to go where the
    // cursor is, and no solve is needed for it.
    let mut group: Vec<usize> = Vec::new();
    for &ordinal in &session.components {
        group.extend(plans.components()[ordinal].iter().copied());
    }
    let touched = system::touched_params(eqs, &group, params);
    for (&p, &cursor) in session.pulled.iter().zip(cursors) {
        if !fixed[p] && !touched.contains(&p) {
            params[p] = cursor;
        }
    }

    let errors_at = |params: &[f64]| -> Vec<f64> {
        session
            .pulled
            .iter()
            .zip(cursors)
            .map(|(&p, &cursor)| cursor - params[p])
            .collect()
    };
    let scale = eqs.length_scale(&group, params).unwrap_or(1.0).max(1.0);
    let (cap, track_tol) = (ADVANCE * scale, TRACK_TOL * scale);

    let mut status = SolveStatus::Converged;
    let mut iterations = 0usize;
    let mut remaining = distance(&errors_at(params));
    for _ in 0..MAX_SUBSTEPS {
        if group.is_empty() || remaining <= track_tol {
            break;
        }
        // Aim no further ahead than the linearisation carries.
        let errors = errors_at(params);
        let aim = (cap / remaining).min(1.0);
        let pulls: Vec<Pull<'_>> = session
            .rows
            .iter()
            .zip(&errors)
            .map(|(row, &error)| Pull {
                row,
                error: error * aim,
            })
            .collect();
        let Some(step) = system::tangent_step(eqs, &group, params, fixed, &pulls, weight_of) else {
            break;
        };
        let saved: Vec<f64> = params.to_vec();
        for (p, delta) in step {
            params[p] += delta;
        }
        let solved = system::solve_planned(
            eqs,
            params,
            fixed,
            plans,
            Some(&session.components),
            options,
            |_, _| {},
        );
        iterations += solved.iterations;
        status = solved.status;
        if !solved.converged {
            // The hard constraints are what the sketch owes the host: this
            // substep did not keep them, so it never happened.
            params.copy_from_slice(&saved);
            break;
        }
        let now = distance(&errors_at(params));
        if now > remaining * (1.0 - GAIN_MIN) {
            // Out of reach, or the linearisation has stopped paying. Keep
            // whichever of the two configurations tracks better.
            if now > remaining {
                params.copy_from_slice(&saved);
            }
            break;
        }
        remaining = now;
    }

    let tracking_error = distance(&errors_at(params));
    DragFrame {
        status,
        residual: crate::numerics::rms(&system::eval_residuals(eqs, &group, params)),
        iterations,
        tracking_error,
        reached: tracking_error <= track_tol,
    }
}
