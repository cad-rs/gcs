//! Rank-revealing diagnosis of a 3D constraint system.
//!
//! Mirrors the 2D [`crate::diagnosis`] module over the analytic 3D
//! Jacobian: the row space is analyzed by a column-pivoted Householder QR of
//! the transposed Jacobian, the pivot columns are the linearly independent
//! constraint rows, and the complement is classified by whether its residual
//! agrees with the least-squares step of the independent rows:
//!
//! - **redundant** — dependent and consistent (removable without changing
//!   the solution);
//! - **conflicting** — dependent and inconsistent (no solution satisfies the
//!   whole system at this linearization);
//! - **partially_redundant** — the gray band between the two.
//!
//! The analysis runs PER COMPONENT of the 3D constraint graph (dependencies
//! can only occur inside a component), and reported indices are positions in
//! the `constraints` slice.
//!
//! # Plane gauge
//!
//! A plane is stored as `(nx, ny, nz, d)`; uniform scaling of the pair
//! leaves every residual invariant, so only three of the four coordinates
//! are physical. The plane's free columns are therefore replaced by the two
//! tangent directions of its unit normal plus the offset direction (the
//! radial scaling direction is dropped), and the DOF count subtracts one
//! gauge dimension per fully free plane. Without this, every free plane
//! would contribute a phantom DOF and a fully constrained plane could show a
//! spurious dependency.
//!
//! # Cylinder gauge
//!
//! A cylinder is stored as an axis point, an axis direction and a radius:
//! seven parameters for five physical freedoms (a line in space has four,
//! plus the radius). Sliding the axis point along the axis and rescaling the
//! direction change the storage but not the cylinder. Unlike the plane
//! gauge these directions cannot simply be dropped, because `Concentric3D`
//! reads the axis point as a centre and thereby pins the slide. They are
//! tested against the Jacobian instead: every gauge direction that no row
//! pins is subtracted from the DOF count. The rank and the redundant /
//! conflicting classification are unaffected (they live in row space).

use crate::diagnosis::Diagnosis;
use crate::numerics::gt;
use crate::space3d::constraint::SpaceConstraint;
use crate::space3d::entity::{SpaceEntity, SpaceEntityKind};
use crate::space3d::solver::{SpaceEquations, SOLVER_NORM_TOL};
use crate::system::{self, Gauge, SparseRow};

pub use crate::diagnosis::{CONSIST_TOL, PARTIAL_TOL, RANK_REL_TOL};

/// Structured diagnosis of a 3D constraint system: the same report as the
/// 2D [`Diagnosis`], with `dof` counting physical freedoms (effective free
/// parameters − Jacobian rank; see the gauge notes in the module docs).
pub type SpaceDiagnosis = Diagnosis;

/// Diagnose a 3D constraint system at the current parameters.
///
/// Constraint indices are the positions in `constraints` (equations stack in
/// order, 1–3 rows per constraint). Reported indices are global.
pub fn diagnose(
    params: &[f64],
    fixed: &[bool],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
) -> SpaceDiagnosis {
    diagnose_equations(&SpaceEquations::new(entities, constraints), params, fixed)
}

/// [`diagnose`] over an equation view, which is how a sketch leaves its
/// suppressed and removed constraints out.
pub(crate) fn diagnose_equations(
    eqs: &SpaceEquations<'_>,
    params: &[f64],
    fixed: &[bool],
) -> SpaceDiagnosis {
    system::diagnose(
        eqs,
        params,
        fixed,
        &SpaceGauge {
            entities: eqs.entities,
        },
    )
}

/// The representation gauge of the 3D entities (see the module docs).
struct SpaceGauge<'a> {
    entities: &'a [SpaceEntity],
}

impl Gauge for SpaceGauge<'_> {
    /// One gauge dimension per fully free plane with a usable normal.
    fn removed_dimensions(&self, params: &[f64], fixed: &[bool]) -> usize {
        self.entities
            .iter()
            .filter(|e| e.kind == SpaceEntityKind::Plane)
            .filter(|e| {
                let idx = [e.param(0), e.param(1), e.param(2), e.param(3)];
                idx.iter().all(|&i| !fixed[i])
                    && plane_tangents(params[idx[0]], params[idx[1]], params[idx[2]]).is_some()
            })
            .count()
    }

    fn effective_columns(
        &self,
        params: &[f64],
        fixed: &[bool],
        touched: &[usize],
    ) -> Option<Vec<SparseRow>> {
        Some(effective_columns(params, fixed, self.entities, touched))
    }

    fn soft_directions(&self, params: &[f64], fixed: &[bool]) -> Vec<SparseRow> {
        cylinder_gauge_vectors(params, fixed, self.entities)
    }
}

/// Representation-gauge directions of the cylinders.
///
/// A cylinder is stored as an axis point `p`, an axis direction `a` and a
/// radius. Moving `p` along the axis and rescaling `a` both describe the
/// same cylinder, so of the six axis parameters only four are physical. Each
/// returned direction is one such move (unit length, in global parameter
/// indices); a direction is only produced when the three parameters it moves
/// are all free and the axis direction is usable.
///
/// Unlike the plane gauge these are not removed up front: `Concentric3D`
/// reads the axis point as a centre, which pins the slide. Whether a
/// direction is still free is decided against the Jacobian
/// (in [`crate::system::diagnose`]).
fn cylinder_gauge_vectors(
    params: &[f64],
    fixed: &[bool],
    entities: &[SpaceEntity],
) -> Vec<SparseRow> {
    let mut out = Vec::new();
    for e in entities {
        if e.kind != SpaceEntityKind::Cylinder {
            continue;
        }
        let point = [e.param(0), e.param(1), e.param(2)];
        let axis = [e.param(3), e.param(4), e.param(5)];
        let (ax, ay, az) = (params[axis[0]], params[axis[1]], params[axis[2]]);
        let len = (ax * ax + ay * ay + az * az).sqrt();
        if !gt(len, SOLVER_NORM_TOL) {
            continue;
        }
        let unit = [ax / len, ay / len, az / len];
        for moved in [point, axis] {
            if moved.iter().all(|&i| !fixed[i]) {
                out.push(moved.iter().copied().zip(unit).collect());
            }
        }
    }
    out
}

/// Effective free directions over the component-local `touched` columns.
///
/// Fully free planes contribute two normal-tangent directions and the offset
/// direction instead of four raw coordinate columns (the radial scaling
/// direction is gauge). Everything else contributes its unit column. The
/// returned vectors are sparse combinations of local column indices.
fn effective_columns(
    params: &[f64],
    fixed: &[bool],
    entities: &[SpaceEntity],
    touched: &[usize],
) -> Vec<Vec<(usize, f64)>> {
    let mut consumed = vec![false; touched.len()];
    let mut columns: Vec<Vec<(usize, f64)>> = Vec::new();
    for e in entities {
        if e.kind != SpaceEntityKind::Plane {
            continue;
        }
        let idx = [e.param(0), e.param(1), e.param(2), e.param(3)];
        let local: Vec<(usize, usize)> = touched
            .iter()
            .enumerate()
            .filter_map(|(lc, &p)| idx.iter().position(|&i| i == p).map(|off| (lc, off)))
            .collect();
        if local.is_empty() {
            continue;
        }
        if idx.iter().any(|&i| fixed[i]) {
            // Partially fixed plane: keep the raw free coordinates.
            continue;
        }
        let Some((t1, t2)) = plane_tangents(params[idx[0]], params[idx[1]], params[idx[2]]) else {
            continue;
        };
        let offset = (0.0, 0.0, 0.0, 1.0);
        let t1 = (t1.0, t1.1, t1.2, 0.0);
        let t2 = (t2.0, t2.1, t2.2, 0.0);
        for basis in [t1, t2, offset] {
            let mut col: Vec<(usize, f64)> = Vec::new();
            for &(lc, off) in &local {
                let coeff = match off {
                    0 => basis.0,
                    1 => basis.1,
                    2 => basis.2,
                    _ => basis.3,
                };
                if coeff != 0.0 {
                    col.push((lc, coeff));
                }
            }
            columns.push(col);
        }
        for &(lc, _) in &local {
            consumed[lc] = true;
        }
    }
    for (lc, used) in consumed.iter().enumerate() {
        if !used {
            columns.push(vec![(lc, 1.0)]);
        }
    }
    columns
}

/// A 3D vector as a plain tuple.
type Vec3 = (f64, f64, f64);

/// Two unit vectors spanning the tangent plane of `n = (nx, ny, nz)`.
/// Deterministic: the helper axis is the coordinate axis least aligned with
/// the normal, so the basis does not jitter between equivalent normals.
fn plane_tangents(nx: f64, ny: f64, nz: f64) -> Option<(Vec3, Vec3)> {
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if !gt(len, SOLVER_NORM_TOL) {
        return None;
    }
    let u = (nx / len, ny / len, nz / len);
    let ax = u.0.abs();
    let ay = u.1.abs();
    let az = u.2.abs();
    let helper = if ax <= ay && ax <= az {
        (1.0, 0.0, 0.0)
    } else if ay <= az {
        (0.0, 1.0, 0.0)
    } else {
        (0.0, 0.0, 1.0)
    };
    let c1 = cross(u, helper);
    let l1 = (c1.0 * c1.0 + c1.1 * c1.1 + c1.2 * c1.2).sqrt();
    if !gt(l1, 0.0) {
        return None;
    }
    let t1 = (c1.0 / l1, c1.1 / l1, c1.2 / l1);
    let t2 = cross(u, t1);
    Some((t1, t2))
}

#[inline]
fn cross(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (
        a.1 * b.2 - a.2 * b.1,
        a.2 * b.0 - a.0 * b.2,
        a.0 * b.1 - a.1 * b.0,
    )
}

/// Per entity, whether a physical motion of it is left (see
/// [`crate::status::movable_entities`]; the gauge of planes and cylinders
/// is not counted as motion).
pub(crate) fn movable_entities(
    eqs: &SpaceEquations<'_>,
    params: &[f64],
    fixed: &[bool],
) -> Vec<bool> {
    let entities = eqs.entities;
    crate::status::movable_entities(
        eqs,
        params,
        fixed,
        &SpaceGauge { entities },
        entities.len(),
        &|e| entities[e].param(0)..entities[e].param(0) + entities[e].param_count(),
    )
}

#[cfg(test)]
mod tests {
    use crate::numerics::householder_qr_column_pivoting;
    use crate::space3d::constraint::SpaceConstraint;
    use crate::space3d::entity::SpacePointRef;
    use crate::space3d::sketch::SpaceSketch;

    #[test]
    fn qr_reveals_rank_and_pivots() {
        let mut a = vec![
            vec![1.0, 1.0, 0.0, 0.0],
            vec![0.0, 1.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0],
        ];
        let (pivots, rank) = householder_qr_column_pivoting(&mut a, 1.0e-10);
        assert_eq!(rank, 2);
        assert_eq!(pivots[..2], [1, 0]);
    }

    #[test]
    fn fully_constrained_plane_has_zero_dof() {
        // PlaneNormal (2 rows) + PointOnPlane (1 row) determine the plane's
        // three physical DOF; the fourth storage coordinate is gauge.
        let mut sk = SpaceSketch::new();
        let plane = sk.add_plane(0.3, -0.4, 0.5, 1.2);
        let p = sk.add_point(0.0, 0.0, 0.0);
        sk.fix_entity(p);
        sk.add_constraint(SpaceConstraint::PlaneNormal {
            plane,
            nx: 0.0,
            ny: 0.0,
            nz: 1.0,
        });
        sk.add_constraint(SpaceConstraint::PointOnPlane {
            point: SpacePointRef::Point(p),
            plane,
        });
        let d = sk.diagnose();
        assert_eq!(d.dof, 0, "plane must be fully constrained: {d:?}");
        assert!(d.conflicting.is_empty(), "{d:?}");
        assert!(
            d.redundant.is_empty(),
            "gauge must not surface as a dependency: {d:?}"
        );
    }

    #[test]
    fn duplicate_sphere_radius_is_redundant() {
        let mut sk = SpaceSketch::new();
        let s = sk.add_sphere(0.0, 0.0, 0.0, 2.0);
        sk.add_constraint(SpaceConstraint::SphereRadius {
            sphere: s,
            radius: 2.0,
        });
        sk.add_constraint(SpaceConstraint::SphereRadius {
            sphere: s,
            radius: 2.0,
        });
        let d = sk.diagnose();
        assert!(!d.redundant.is_empty(), "{d:?}");
        assert!(d.conflicting.is_empty(), "{d:?}");
        assert!(
            d.dependency_groups.iter().any(|g| g == &vec![0, 1]),
            "duplicate rows form one dependency group: {d:?}"
        );
    }

    #[test]
    fn contradiction_is_conflicting() {
        let mut sk = SpaceSketch::new();
        let p = sk.add_point(0.0, 0.0, 0.0);
        sk.add_constraint(SpaceConstraint::fix_point(p, 0.0, 0.0, 0.0));
        sk.add_constraint(SpaceConstraint::fix_point(p, 1.0, 1.0, 1.0));
        let d = sk.diagnose();
        assert!(!d.conflicting.is_empty(), "{d:?}");
    }
}
