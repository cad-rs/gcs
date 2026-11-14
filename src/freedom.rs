//! What an under-defined entity can still do, in geometric terms.
//!
//! [`crate::Status`] says that an entity is under-defined; an
//! [`EntityFreedom`] says how: "slides along this direction", "turns about
//! this point", "its radius can change". It describes the motions of the
//! entity that no constraint resists (the null space of the constraint
//! Jacobian, restricted to the entity) twice:
//!
//! - [`EntityFreedom::motions`]: while the rest of the system may move
//!   along. It is empty exactly when the entity is not under-defined in
//!   [`crate::Status`]. A body of an assembly that floats free with its
//!   neighbour has all six here.
//! - [`EntityFreedom::own`]: with everything else held where it is. This is
//!   the joint a body is left with relative to its neighbours (a pin, a
//!   slider, a ball), and what a constraint on the entity alone would have
//!   to take away. An entity whose freedom it shares with others (a point
//!   in the middle of a bent open chain) moves in `motions` and not here.
//!
//! A space of motions of more than one dimension has many bases. The one
//! given is chosen in a fixed order: translations first (along the world
//! axes when all of them are free), then rigid turns (about one point when
//! they share one), then changes of one measure of the shape each
//! ([`Measure`]), and what is left as [`Motion2::Other`] /
//! [`Motion3::Other`] with its raw velocity. Each motion is a direction:
//! its sign means nothing.
//!
//! Like the diagnosis this is first order, at the current parameters: read
//! it at a solved configuration. At a singular one (an exact tangency, a
//! straight chain) the first-order motions can be more than, or other than,
//! what the geometry can do over a finite distance.

use glam::{DQuat, DVec2, DVec3};

/// A measure of an entity's shape that a motion can change on its own.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Measure {
    /// A line's length.
    Length,
    /// A circle's, an arc's, a sphere's or a cylinder's radius.
    Radius,
    /// A 2D arc's angular span (end angle − start angle).
    Sweep,
}

/// One way a 2D entity can move (see the [module docs](self)).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Motion2 {
    /// Slides along `direction` (unit length) without turning.
    Translate { direction: DVec2 },
    /// Turns rigidly about `center`.
    Rotate { center: DVec2 },
    /// The measure `what` changes, with the rest of the shape kept.
    Resize { what: Measure },
    /// A motion that is none of the above (several measures changing
    /// together, or a measure changing as the entity moves in a way it
    /// cannot alone): its velocity per stored parameter of the entity,
    /// unit length.
    Other { velocity: Vec<f64> },
}

/// One way a 3D entity or an assembly body can move (see the
/// [module docs](self)).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Motion3 {
    /// Slides along `direction` (unit length) without turning.
    Translate { direction: DVec3 },
    /// Turns rigidly about the line through `point` along `axis` (unit
    /// length).
    Rotate { axis: DVec3, point: DVec3 },
    /// Turns about the line through `point` along `axis` (unit length)
    /// while sliding along it, `pitch` units of length per radian.
    Screw {
        axis: DVec3,
        point: DVec3,
        pitch: f64,
    },
    /// The measure `what` changes, with the rest of the shape kept.
    Resize { what: Measure },
    /// A motion that is none of the above: its velocity per stored
    /// parameter of the entity (for an assembly body its six coordinates:
    /// translation along world X, Y, Z, then rotation about its own x, y,
    /// z axes), unit length.
    Other { velocity: Vec<f64> },
}

/// How an entity can still move (see the [module docs](self)).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EntityFreedom<M> {
    /// The motions of the entity while the rest of the system may move
    /// along; empty exactly when the entity is not under-defined.
    pub motions: Vec<M>,
    /// The motions of the entity with everything else held.
    pub own: Vec<M>,
}

impl<M> Default for EntityFreedom<M> {
    fn default() -> Self {
        Self {
            motions: Vec::new(),
            own: Vec::new(),
        }
    }
}

// ── the description of one space of motions ─────────────────────────────────

/// How rigid motions and shape measures act on the stored parameters of
/// one entity, at its current parameters.
pub(crate) struct Frame {
    /// The point the turns are about (a line's midpoint, a circle's
    /// centre).
    pub(crate) reference: DVec3,
    /// The change of the parameters when the entity slides along each
    /// world axis (two in 2D, three in 3D).
    pub(crate) translations: Vec<Vec<f64>>,
    /// The change of the parameters when the entity turns about each world
    /// axis (`.0`; only Z in 2D) through `reference`.
    pub(crate) rotations: Vec<(DVec3, Vec<f64>)>,
    /// The rate of change of each measure of the shape (a row over the
    /// parameters). Rigid motions leave them all unchanged, and together
    /// with the rigid motions they account for every motion.
    pub(crate) measures: Vec<(Measure, Vec<f64>)>,
    /// Directions of the parameters that store the same entity (a plane's
    /// normal and offset scaled together): taken out before anything else.
    pub(crate) gauge: Vec<Vec<f64>>,
}

/// One motion in a dimension-neutral form.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Described {
    Translate(DVec3),
    Rotate {
        axis: DVec3,
        point: DVec3,
    },
    Screw {
        axis: DVec3,
        point: DVec3,
        pitch: f64,
    },
    Resize(Measure),
    Other(Vec<f64>),
}

/// Below this a vector of a basis is taken for zero, relative to entries of
/// order one.
const DIM_TOL: f64 = 1.0e-9;
/// Two spaces meet along a direction whose principal angle has a sine below
/// this.
const MEET_TOL: f64 = 1.0e-6;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// `v` less its components along the orthonormal `basis`.
fn reject(mut v: Vec<f64>, basis: &[Vec<f64>]) -> Vec<f64> {
    for _ in 0..2 {
        for u in basis {
            let d = dot(&v, u);
            for (a, b) in v.iter_mut().zip(u) {
                *a -= d * b;
            }
        }
    }
    v
}

/// An orthonormal basis of the span of `vs`, in their order.
fn orthonormal(vs: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut out: Vec<Vec<f64>> = Vec::new();
    for v in vs {
        let scale = norm(v);
        if scale <= DIM_TOL {
            continue;
        }
        let w = reject(v.iter().map(|x| x / scale).collect(), &out);
        let n = norm(&w);
        if n > DIM_TOL {
            out.push(w.iter().map(|x| x / n).collect());
        }
    }
    out
}

/// Eigenvalues and eigenvectors (columns of the second) of the symmetric
/// matrix `a`, by cyclic Jacobi rotations.
fn symmetric_eigen(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v = vec![vec![0.0; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let total: f64 = a.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum::<f64>()
            .sqrt();
        if off <= 1.0e-15 * total.max(f64::MIN_POSITIVE) {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vp, vq) = (row[p], row[q]);
                    row[p] = c * vp - s * vq;
                    row[q] = s * vp + c * vq;
                }
            }
        }
    }
    ((0..n).map(|i| a[i][i]).collect(), v)
}

/// The least-squares solution of smallest norm of `Σ c_k columns[k] = b`.
fn min_norm_solve(columns: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let m = columns.len();
    let normal: Vec<Vec<f64>> = (0..m)
        .map(|i| (0..m).map(|j| dot(&columns[i], &columns[j])).collect())
        .collect();
    let rhs: Vec<f64> = columns.iter().map(|c| dot(c, b)).collect();
    let (values, vectors) = symmetric_eigen(normal);
    let top = values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let mut c = vec![0.0; m];
    for (k, &lambda) in values.iter().enumerate() {
        if lambda <= 1.0e-12 * top || lambda <= 0.0 {
            continue;
        }
        let along: f64 = (0..m).map(|i| vectors[i][k] * rhs[i]).sum::<f64>() / lambda;
        for (i, ci) in c.iter_mut().enumerate() {
            *ci += along * vectors[i][k];
        }
    }
    c
}

/// An orthonormal basis of the part of the space with orthonormal basis
/// `s` that lies in the span of `u`.
fn meet(s: &[Vec<f64>], u: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let u = orthonormal(u);
    if s.is_empty() || u.is_empty() {
        return Vec::new();
    }
    let d = s.len();
    // K = I − (SᵀU)(UᵀS): its eigenvalues are the squared sines of the
    // principal angles between the two spaces.
    let c: Vec<Vec<f64>> = s
        .iter()
        .map(|si| u.iter().map(|uj| dot(si, uj)).collect())
        .collect();
    let k: Vec<Vec<f64>> = (0..d)
        .map(|i| {
            (0..d)
                .map(|j| f64::from(u8::from(i == j)) - dot(&c[i], &c[j]))
                .collect()
        })
        .collect();
    let (values, vectors) = symmetric_eigen(k);
    let mut order: Vec<usize> = (0..d).collect();
    order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
    let picked: Vec<Vec<f64>> = order
        .into_iter()
        .filter(|&k| values[k] <= MEET_TOL * MEET_TOL)
        .map(|k| {
            let mut w = vec![0.0; s[0].len()];
            for (i, si) in s.iter().enumerate() {
                for (a, b) in w.iter_mut().zip(si) {
                    *a += vectors[i][k] * b;
                }
            }
            w
        })
        .collect();
    orthonormal(&picked)
}

/// A unit vector of the line it spans, its largest entry made positive.
fn canonical_sign(v: DVec3) -> DVec3 {
    let v = v.normalize_or_zero();
    let big = [v.x, v.y, v.z]
        .into_iter()
        .fold(0.0_f64, |m, x| if x.abs() > m.abs() { x } else { m });
    if big < 0.0 {
        -v
    } else {
        v
    }
}

/// An orthonormal basis of the span of `ds` (3D vectors): the world axes
/// it contains first, as far as they reach into it.
fn world_basis(ds: &[DVec3]) -> Vec<DVec3> {
    let span = orthonormal(&ds.iter().map(|d| d.to_array().to_vec()).collect::<Vec<_>>());
    let project = |v: DVec3| -> DVec3 {
        span.iter().fold(DVec3::ZERO, |acc, u| {
            let u = DVec3::from_slice(u);
            acc + u * u.dot(v)
        })
    };
    let mut out: Vec<DVec3> = Vec::new();
    for axis in [DVec3::X, DVec3::Y, DVec3::Z] {
        if out.len() == span.len() {
            break;
        }
        let mut w = project(axis);
        for &u in &out {
            w -= u * u.dot(w);
        }
        if w.length() > 1.0e-6 {
            out.push(canonical_sign(w));
        }
    }
    out
}

fn skew(w: DVec3) -> [[f64; 3]; 3] {
    [[0.0, -w.z, w.y], [w.z, 0.0, -w.x], [-w.y, w.x, 0.0]]
}

/// The motions of the space with basis `space` (over the entity's stored
/// parameters), described with `frame` (see the [module docs](self)).
pub(crate) fn describe(frame: &Frame, space: &[Vec<f64>]) -> Vec<Described> {
    let gauge = orthonormal(&frame.gauge);
    let s = orthonormal(
        &space
            .iter()
            .map(|v| reject(v.clone(), &gauge))
            .collect::<Vec<_>>(),
    );
    if s.is_empty() {
        return Vec::new();
    }
    let translations: Vec<Vec<f64>> = frame
        .translations
        .iter()
        .map(|t| reject(t.clone(), &gauge))
        .collect();
    let rotations: Vec<Vec<f64>> = frame
        .rotations
        .iter()
        .map(|(_, r)| reject(r.clone(), &gauge))
        .collect();
    let rigid: Vec<Vec<f64>> = translations.iter().chain(&rotations).cloned().collect();
    let world = [DVec3::X, DVec3::Y, DVec3::Z];
    let mut out = Vec::new();

    // Translations.
    let st = meet(&s, &translations);
    let dirs: Vec<DVec3> = st
        .iter()
        .map(|v| {
            min_norm_solve(&translations, v)
                .iter()
                .zip(world)
                .fold(DVec3::ZERO, |acc, (&c, axis)| acc + axis * c)
        })
        .collect();
    let dirs = world_basis(&dirs);
    out.extend(dirs.iter().map(|&d| Described::Translate(d)));

    // Rigid turns beyond the translations, as twists (velocity of the
    // reference point, angular velocity) less the translations at hand.
    let sr = meet(&s, &rigid);
    let turns = orthonormal(
        &sr.iter()
            .map(|v| reject(v.clone(), &st))
            .collect::<Vec<_>>(),
    );
    let nt = translations.len();
    let mut twists: Vec<(DVec3, DVec3)> = turns
        .iter()
        .map(|r| {
            let c = min_norm_solve(&rigid, r);
            let v = c[..nt]
                .iter()
                .zip(world)
                .fold(DVec3::ZERO, |acc, (&c, axis)| acc + axis * c);
            let w = c[nt..]
                .iter()
                .zip(&frame.rotations)
                .fold(DVec3::ZERO, |acc, (&c, (axis, _))| acc + *axis * c);
            let v = dirs.iter().fold(v, |v, &d| v - d * d.dot(v));
            (v, w)
        })
        .filter(|(_, w)| w.length() > DIM_TOL)
        .collect();
    if !twists.is_empty() {
        // A point about which every turn is a pure rotation, up to the
        // translations at hand: P⊥(v_i + ω_i × x) = 0 for x = point −
        // reference, P⊥ taking out the free translation directions.
        let perp = |u: DVec3| dirs.iter().fold(u, |u, &d| u - d * d.dot(u));
        let mut columns = vec![Vec::new(); 3];
        let mut rhs = Vec::new();
        for &(v, w) in &twists {
            let m = skew(w);
            for k in 0..3 {
                let col = perp(DVec3::new(m[0][k], m[1][k], m[2][k]));
                columns[k].extend(col.to_array());
            }
            rhs.extend(perp(-v).to_array());
        }
        let x = min_norm_solve(&columns, &rhs);
        let residual: f64 = (0..rhs.len())
            .map(|i| {
                let fit: f64 = (0..3).map(|k| columns[k][i] * x[k]).sum();
                (fit - rhs[i]).powi(2)
            })
            .sum::<f64>()
            .sqrt();
        let size: f64 = twists
            .iter()
            .map(|(v, w)| v.length() + w.length() * (1.0 + frame.reference.length()))
            .sum();
        if residual <= MEET_TOL * size {
            let point = frame.reference + DVec3::new(x[0], x[1], x[2]);
            let axes = world_basis(&twists.iter().map(|&(_, w)| w).collect::<Vec<_>>());
            out.extend(
                axes.into_iter()
                    .map(|axis| Described::Rotate { axis, point }),
            );
        } else {
            // Screws, with the angular velocities made orthonormal.
            for i in 0..twists.len() {
                for j in 0..i {
                    let (vj, wj) = twists[j];
                    let along = twists[i].1.dot(wj);
                    twists[i].0 -= vj * along;
                    twists[i].1 -= wj * along;
                }
                let len = twists[i].1.length();
                twists[i].0 /= len;
                twists[i].1 /= len;
            }
            for (v, w) in twists {
                let pitch = w.dot(v);
                let point = frame.reference + w.cross(v);
                let axis = canonical_sign(w);
                if pitch.abs() <= MEET_TOL * (v.length() + 1.0 + frame.reference.length()) {
                    out.push(Described::Rotate { axis, point });
                } else {
                    out.push(Described::Screw { axis, point, pitch });
                }
            }
        }
    }

    // Changes of shape: what is left, one measure at a time where it can.
    let q = orthonormal(&s.iter().map(|v| reject(v.clone(), &sr)).collect::<Vec<_>>());
    if !q.is_empty() {
        // rates[k][i]: how fast measure k changes along q[i].
        let rates: Vec<Vec<f64>> = frame
            .measures
            .iter()
            .map(|(_, row)| q.iter().map(|qi| dot(row, qi)).collect())
            .collect();
        let mut resized: Vec<Vec<f64>> = Vec::new();
        for (k, (what, _)) in frame.measures.iter().enumerate() {
            // A combination of q that changes measure k alone.
            let target: Vec<f64> = (0..frame.measures.len())
                .map(|j| f64::from(u8::from(j == k)))
                .collect();
            // Columns of Mᵀ: one per q[i], the rates of every measure.
            let columns: Vec<Vec<f64>> = (0..q.len())
                .map(|i| rates.iter().map(|r| r[i]).collect())
                .collect();
            let c = min_norm_solve(&columns, &target);
            let fit: Vec<f64> = (0..target.len())
                .map(|j| (0..q.len()).map(|i| columns[i][j] * c[i]).sum())
                .collect();
            let miss = norm(
                &fit.iter()
                    .zip(&target)
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>(),
            );
            if miss <= MEET_TOL {
                out.push(Described::Resize(*what));
                let mut v = vec![0.0; q[0].len()];
                for (qi, &ci) in q.iter().zip(&c) {
                    for (a, b) in v.iter_mut().zip(qi) {
                        *a += ci * b;
                    }
                }
                resized.push(v);
            }
        }
        let resized = orthonormal(&resized);
        let rest = orthonormal(
            &q.iter()
                .map(|v| reject(v.clone(), &resized))
                .collect::<Vec<_>>(),
        );
        out.extend(rest.into_iter().map(Described::Other));
    }
    out
}

// ── frames of the entity kinds ──────────────────────────────────────────────

fn unit(n: usize, k: usize) -> Vec<f64> {
    let mut v = vec![0.0; n];
    v[k] = 1.0;
    v
}

/// The frame of a 2D entity of kind `kind` with parameters `p`.
pub(crate) fn frame_2d(kind: crate::EntityKind, p: &[f64]) -> Frame {
    use crate::EntityKind::*;
    let n = p.len();
    // ω × (x − c) for ω along +Z.
    let turn = |x: f64, y: f64, c: DVec2| [-(y - c.y), x - c.x];
    match kind {
        Point => Frame {
            reference: DVec3::new(p[0], p[1], 0.0),
            translations: vec![unit(n, 0), unit(n, 1)],
            rotations: vec![(DVec3::Z, vec![0.0; n])],
            measures: Vec::new(),
            gauge: Vec::new(),
        },
        Line => {
            let mid = DVec2::new(p[0] + p[2], p[1] + p[3]) * 0.5;
            let (a, b) = (turn(p[0], p[1], mid), turn(p[2], p[3], mid));
            let d = DVec2::new(p[2] - p[0], p[3] - p[1]).normalize_or_zero();
            Frame {
                reference: mid.extend(0.0),
                translations: vec![vec![1.0, 0.0, 1.0, 0.0], vec![0.0, 1.0, 0.0, 1.0]],
                rotations: vec![(DVec3::Z, vec![a[0], a[1], b[0], b[1]])],
                measures: vec![(Measure::Length, vec![-d.x, -d.y, d.x, d.y])],
                gauge: Vec::new(),
            }
        }
        Circle => Frame {
            reference: DVec3::new(p[0], p[1], 0.0),
            translations: vec![unit(n, 0), unit(n, 1)],
            rotations: vec![(DVec3::Z, vec![0.0; n])],
            measures: vec![(Measure::Radius, unit(n, 2))],
            gauge: Vec::new(),
        },
        Arc => Frame {
            reference: DVec3::new(p[0], p[1], 0.0),
            translations: vec![unit(n, 0), unit(n, 1)],
            rotations: vec![(DVec3::Z, vec![0.0, 0.0, 0.0, 1.0, 1.0])],
            measures: vec![
                (Measure::Radius, unit(n, 2)),
                (Measure::Sweep, vec![0.0, 0.0, 0.0, -1.0, 1.0]),
            ],
            gauge: Vec::new(),
        },
    }
}

/// The frame of a 3D entity of kind `kind` with parameters `p`.
pub(crate) fn frame_3d(kind: crate::space3d::SpaceEntityKind, p: &[f64]) -> Frame {
    use crate::space3d::SpaceEntityKind::*;
    let n = p.len();
    let world = [DVec3::X, DVec3::Y, DVec3::Z];
    let at = |k: usize| DVec3::new(p[k], p[k + 1], p[k + 2]);
    // Translations moving the points that start at the offsets `points`.
    let moving = |points: &[usize]| -> Vec<Vec<f64>> {
        world
            .iter()
            .map(|axis| {
                let mut v = vec![0.0; n];
                for &k in points {
                    v[k..k + 3].copy_from_slice(&axis.to_array());
                }
                v
            })
            .collect()
    };
    // Turns about `c` moving the points at `points` and turning the
    // directions at `directions`.
    let turning = |c: DVec3, points: &[usize], directions: &[usize]| -> Vec<(DVec3, Vec<f64>)> {
        world
            .iter()
            .map(|&axis| {
                let mut v = vec![0.0; n];
                for &k in points {
                    v[k..k + 3].copy_from_slice(&axis.cross(at(k) - c).to_array());
                }
                for &k in directions {
                    v[k..k + 3].copy_from_slice(&axis.cross(at(k)).to_array());
                }
                (axis, v)
            })
            .collect()
    };
    match kind {
        SpacePoint => Frame {
            reference: at(0),
            translations: moving(&[0]),
            rotations: turning(at(0), &[0], &[]),
            measures: Vec::new(),
            gauge: Vec::new(),
        },
        SpaceLine => {
            let mid = (at(0) + at(3)) * 0.5;
            let d = (at(3) - at(0)).normalize_or_zero();
            let mut length = vec![0.0; n];
            length[0..3].copy_from_slice(&(-d).to_array());
            length[3..6].copy_from_slice(&d.to_array());
            Frame {
                reference: mid,
                translations: moving(&[0, 3]),
                rotations: turning(mid, &[0, 3], &[]),
                measures: vec![(Measure::Length, length)],
                gauge: Vec::new(),
            }
        }
        Plane => {
            // n · x = d: sliding by t adds n · t to d; turning about the
            // foot point of the origin turns n and keeps d.
            let normal = at(0);
            let foot = if normal.length_squared() > 0.0 {
                normal * (p[3] / normal.length_squared())
            } else {
                DVec3::ZERO
            };
            let translations = world
                .iter()
                .map(|axis| vec![0.0, 0.0, 0.0, normal.dot(*axis)])
                .collect();
            Frame {
                reference: foot,
                translations,
                rotations: turning(foot, &[], &[0]),
                measures: Vec::new(),
                gauge: vec![p.to_vec()],
            }
        }
        Sphere | SpaceCircle | SpaceArc => Frame {
            reference: at(0),
            translations: moving(&[0]),
            rotations: turning(at(0), &[0], &[]),
            measures: vec![(Measure::Radius, unit(n, 3))],
            gauge: Vec::new(),
        },
        Cylinder => {
            // Axis point, axis direction, radius; the point sliding along
            // the axis and the direction rescaled are the same cylinder.
            let axis = at(3).normalize_or_zero();
            let mut slide = vec![0.0; n];
            slide[0..3].copy_from_slice(&axis.to_array());
            let mut scale = vec![0.0; n];
            scale[3..6].copy_from_slice(&axis.to_array());
            Frame {
                reference: at(0),
                translations: moving(&[0]),
                rotations: turning(at(0), &[0], &[3]),
                measures: vec![(Measure::Radius, unit(n, 6))],
                gauge: vec![slide, scale],
            }
        }
        SpaceSpline => {
            let points: Vec<usize> = (0..n / 3).map(|i| 3 * i).collect();
            let centroid = if points.is_empty() {
                DVec3::ZERO
            } else {
                points.iter().fold(DVec3::ZERO, |acc, &k| acc + at(k)) / points.len() as f64
            };
            Frame {
                reference: centroid,
                translations: moving(&points),
                rotations: turning(centroid, &points, &[]),
                measures: Vec::new(),
                gauge: Vec::new(),
            }
        }
    }
}

/// The frame of an assembly body at `position`, `rotation`: six columns,
/// translation along world X, Y, Z and rotation about its own x, y, z
/// axes.
pub(crate) fn frame_body(position: DVec3, rotation: DQuat) -> Frame {
    let world = [DVec3::X, DVec3::Y, DVec3::Z];
    let local = rotation.inverse();
    Frame {
        reference: position,
        translations: (0..3).map(|k| unit(6, k)).collect(),
        rotations: world
            .iter()
            .map(|&axis| {
                let w = local * axis;
                (axis, vec![0.0, 0.0, 0.0, w.x, w.y, w.z])
            })
            .collect(),
        measures: Vec::new(),
        gauge: Vec::new(),
    }
}

/// The 2D form of a description.
pub(crate) fn to_2d(d: Described) -> Motion2 {
    match d {
        Described::Translate(v) => Motion2::Translate {
            direction: v.truncate(),
        },
        Described::Rotate { point, .. } | Described::Screw { point, .. } => Motion2::Rotate {
            center: point.truncate(),
        },
        Described::Resize(what) => Motion2::Resize { what },
        Described::Other(velocity) => Motion2::Other { velocity },
    }
}

/// The 3D form of a description.
pub(crate) fn to_3d(d: Described) -> Motion3 {
    match d {
        Described::Translate(direction) => Motion3::Translate { direction },
        Described::Rotate { axis, point } => Motion3::Rotate { axis, point },
        Described::Screw { axis, point, pitch } => Motion3::Screw { axis, point, pitch },
        Described::Resize(what) => Motion3::Resize { what },
        Described::Other(velocity) => Motion3::Other { velocity },
    }
}

/// The freedom of every entity from its motion spaces: `frame(e)` is the
/// frame of entity `e`, `None` for a removed one.
pub(crate) fn freedoms_of<M>(
    spaces: &[crate::status::MotionSpace],
    frame: &dyn Fn(usize) -> Option<Frame>,
    form: fn(Described) -> M,
) -> Vec<Option<EntityFreedom<M>>> {
    spaces
        .iter()
        .enumerate()
        .map(|(e, space)| {
            let frame = frame(e)?;
            Some(EntityFreedom {
                motions: describe(&frame, &space.follow)
                    .into_iter()
                    .map(form)
                    .collect(),
                own: describe(&frame, &space.own).into_iter().map(form).collect(),
            })
        })
        .collect()
}
