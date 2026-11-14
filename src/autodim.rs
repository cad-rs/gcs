//! Picking dimensions that take away the freedoms left (automatic
//! dimensioning). See `Sketch::auto_dimensions` and
//! `SpaceSketch::auto_dimensions`.
//!
//! The null space of the constraint Jacobian at the current parameters —
//! over the effective columns of the gauge, the soft gauge directions that
//! nothing pins taken out, so its dimension is the `dof` of the diagnosis —
//! is kept as an orthonormal basis of sparse vectors over the stored
//! parameters. A candidate is taken when every one of its rows has a part
//! in that space no earlier row has: its rows are then independent of the
//! system, it takes away as many freedoms as it has equations, and the
//! directions it takes are removed from the basis (only the basis vectors
//! its rows touch change). Candidates are offered in a fixed order and set
//! to what they measure, so a set taken this way holds where the geometry
//! stands and leaves nothing redundant.

use std::collections::HashMap;

use crate::system::{self, Equations, Gauge, SparseRow};

/// A part of a row in the null space below this, relative to the row,
/// does not count: it would be a near-dependency to the diagnosis.
const TAKE_TOL: f64 = 1.0e-6;

/// A sparse vector over stored parameters, ascending.
type Sparse = Vec<(usize, f64)>;

fn sparse_dot(a: &Sparse, b: &[(usize, f64)]) -> f64 {
    // `b` is short (a row); look each of its entries up in `a`.
    b.iter()
        .map(|&(p, v)| {
            a.binary_search_by_key(&p, |&(q, _)| q)
                .map_or(0.0, |i| a[i].1 * v)
        })
        .sum()
}

/// The physical null space of a system, kept orthonormal.
pub(crate) struct NullBasis {
    basis: Vec<Sparse>,
}

/// Adds `c · v` to the accumulator `acc` (parameter → value).
fn accumulate(acc: &mut HashMap<usize, f64>, v: &Sparse, c: f64) {
    for &(p, x) in v {
        *acc.entry(p).or_insert(0.0) += c * x;
    }
}

fn to_sparse(acc: HashMap<usize, f64>) -> Sparse {
    let mut v: Sparse = acc.into_iter().filter(|&(_, x)| x != 0.0).collect();
    v.sort_unstable_by_key(|&(p, _)| p);
    v
}

impl NullBasis {
    /// The null space of `eqs` at `params` (see the module docs).
    pub(crate) fn of<E: Equations + ?Sized, G: Gauge>(
        eqs: &E,
        params: &[f64],
        fixed: &[bool],
        gauge: &G,
    ) -> Self {
        let tol = crate::diagnosis::RANK_REL_TOL;
        let n_cols = eqs.column_count(params);
        let mut touched_any = vec![false; n_cols];
        let mut basis: Vec<Sparse> = Vec::new();
        for group in system::components(eqs, params) {
            let mut touched = system::touched_params(eqs, &group, params);
            for &p in &touched {
                touched_any[p] = true;
            }
            touched.retain(|&p| !fixed[p]);
            if touched.is_empty() {
                continue;
            }
            let n_eq = system::equation_total(eqs, &group);
            let col_of: HashMap<usize, usize> =
                touched.iter().enumerate().map(|(c, &p)| (p, c)).collect();
            let mut j_raw = vec![vec![0.0_f64; touched.len()]; n_eq];
            system::analytic_jacobian(eqs, &group, params, &col_of, &mut j_raw);
            let columns: Vec<SparseRow> = gauge
                .effective_columns(params, fixed, &touched)
                .unwrap_or_else(|| (0..touched.len()).map(|c| vec![(c, 1.0)]).collect());
            let j: Vec<Vec<f64>> = j_raw
                .iter()
                .map(|row| {
                    columns
                        .iter()
                        .map(|col| col.iter().map(|&(lc, w)| row[lc] * w).sum())
                        .collect()
                })
                .collect();
            for b in crate::numerics::null_space_basis(&j, tol) {
                let mut acc = HashMap::new();
                for (k, &coeff) in b.iter().enumerate() {
                    if coeff != 0.0 {
                        for &(lc, w) in &columns[k] {
                            *acc.entry(touched[lc]).or_insert(0.0) += coeff * w;
                        }
                    }
                }
                basis.push(to_sparse(acc));
            }
        }
        // Free parameters no constraint touches, over the effective columns
        // (a free plane's normal length is not a motion).
        let untouched: Vec<usize> = (0..n_cols)
            .filter(|&p| !fixed[p] && !touched_any[p])
            .collect();
        if !untouched.is_empty() {
            let columns = gauge
                .effective_columns(params, fixed, &untouched)
                .unwrap_or_else(|| (0..untouched.len()).map(|c| vec![(c, 1.0)]).collect());
            for col in columns {
                let mut acc = HashMap::new();
                for (lc, w) in col {
                    *acc.entry(untouched[lc]).or_insert(0.0) += w;
                }
                basis.push(to_sparse(acc));
            }
        }
        let mut out = Self {
            basis: Self::orthonormal(basis),
        };
        // Soft gauge directions nothing pins are not freedoms.
        for g in gauge.soft_directions(params, fixed) {
            let w: Vec<f64> = out.basis.iter().map(|b| sparse_dot(b, &g)).collect();
            let inside = w.iter().map(|x| x * x).sum::<f64>().sqrt();
            if inside > 1.0 - TAKE_TOL {
                out.remove(&[w]);
            }
        }
        out
    }

    /// Orthonormalizes `vs` (Gram–Schmidt twice), dropping dependent ones.
    fn orthonormal(vs: Vec<Sparse>) -> Vec<Sparse> {
        let mut out: Vec<Sparse> = Vec::new();
        for v in vs {
            let mut acc: HashMap<usize, f64> = v.iter().copied().collect();
            let scale = v.iter().map(|&(_, x)| x * x).sum::<f64>().sqrt();
            if scale == 0.0 {
                continue;
            }
            for _ in 0..2 {
                let current = to_sparse(acc.clone());
                for u in &out {
                    let d = sparse_dot(u, &current);
                    if d != 0.0 {
                        accumulate(&mut acc, u, -d);
                    }
                }
            }
            let w = to_sparse(acc);
            let n = w.iter().map(|&(_, x)| x * x).sum::<f64>().sqrt();
            if n > crate::diagnosis::RANK_REL_TOL * scale {
                out.push(w.into_iter().map(|(p, x)| (p, x / n)).collect());
            }
        }
        out
    }

    /// The number of freedoms left.
    pub(crate) fn dimension(&self) -> usize {
        self.basis.len()
    }

    /// Whether stored parameter `p` takes part in some freedom.
    pub(crate) fn moves(&self, p: usize) -> bool {
        self.basis
            .iter()
            .any(|b| b.binary_search_by_key(&p, |&(q, _)| q).is_ok())
    }

    /// Takes the candidate with rows `rows` (over stored parameters, held
    /// ones left out) if every row has a part in the null space no earlier
    /// row has; then removes what they take. Returns whether it did.
    pub(crate) fn take(&mut self, rows: &[SparseRow]) -> bool {
        if rows.is_empty() || self.basis.is_empty() {
            return false;
        }
        // Coordinates of each row in the basis.
        let w: Vec<Vec<f64>> = rows
            .iter()
            .map(|r| self.basis.iter().map(|b| sparse_dot(b, r)).collect())
            .collect();
        // Every row independent of the earlier ones, within the basis.
        let mut ortho: Vec<Vec<f64>> = Vec::new();
        for (r, wr) in rows.iter().zip(&w) {
            let size = r.iter().map(|&(_, v)| v * v).sum::<f64>().sqrt();
            let mut v = wr.clone();
            for _ in 0..2 {
                for u in &ortho {
                    let d: f64 = v.iter().zip(u).map(|(a, b)| a * b).sum();
                    for (a, b) in v.iter_mut().zip(u) {
                        *a -= d * b;
                    }
                }
            }
            let n = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if !(size > 0.0 && n > TAKE_TOL * size) {
                return false;
            }
            ortho.push(v.iter().map(|x| x / n).collect());
        }
        self.remove(&ortho);
        true
    }

    /// Removes the directions `w` (coordinates in the basis, orthonormal)
    /// from the space: the basis vectors they involve are replaced by an
    /// orthonormal basis of what is left of their span.
    fn remove(&mut self, w: &[Vec<f64>]) {
        let mut ortho: Vec<Vec<f64>> = Vec::new();
        for v in w {
            let mut v = v.clone();
            for u in &ortho {
                let d: f64 = v.iter().zip(u).map(|(a, b)| a * b).sum();
                for (a, b) in v.iter_mut().zip(u) {
                    *a -= d * b;
                }
            }
            let n = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if n > 0.0 {
                ortho.push(v.iter().map(|x| x / n).collect());
            }
        }
        let involved: Vec<usize> = (0..self.basis.len())
            .filter(|&j| ortho.iter().any(|u| u[j].abs() > 1.0e-15))
            .collect();
        let m = involved.len();
        // The rest of span(e_j : j involved) after taking out `ortho`,
        // restricted to those coordinates.
        let restricted: Vec<Vec<f64>> = ortho
            .iter()
            .map(|u| involved.iter().map(|&j| u[j]).collect())
            .collect();
        let mut keep: Vec<Vec<f64>> = Vec::new();
        let want = m.saturating_sub(restricted.len());
        // Pivoted Gram–Schmidt over the unit vectors: take the one with the
        // most left each time.
        let mut candidates: Vec<Vec<f64>> = (0..m)
            .map(|i| {
                let mut e = vec![0.0; m];
                e[i] = 1.0;
                for _ in 0..2 {
                    for u in restricted.iter().chain(&keep) {
                        let d: f64 = e.iter().zip(u).map(|(a, b)| a * b).sum();
                        for (a, b) in e.iter_mut().zip(u) {
                            *a -= d * b;
                        }
                    }
                }
                e
            })
            .collect();
        while keep.len() < want {
            let (best, norm) = candidates
                .iter()
                .enumerate()
                .map(|(i, c)| (i, c.iter().map(|x| x * x).sum::<f64>().sqrt()))
                .fold((0, -1.0), |acc, x| if x.1 > acc.1 { x } else { acc });
            if norm <= 1.0e-12 {
                break;
            }
            let u: Vec<f64> = candidates[best].iter().map(|x| x / norm).collect();
            for c in &mut candidates {
                for _ in 0..2 {
                    let d: f64 = c.iter().zip(&u).map(|(a, b)| a * b).sum();
                    for (a, b) in c.iter_mut().zip(&u) {
                        *a -= d * b;
                    }
                }
            }
            keep.push(u);
        }
        let new: Vec<Sparse> = keep
            .iter()
            .map(|c| {
                let mut acc = HashMap::new();
                for (k, &j) in involved.iter().enumerate() {
                    if c[k] != 0.0 {
                        accumulate(&mut acc, &self.basis[j], c[k]);
                    }
                }
                to_sparse(acc)
            })
            .collect();
        let mut rest: Vec<Sparse> = (0..self.basis.len())
            .filter(|j| !involved.contains(j))
            .map(|j| std::mem::take(&mut self.basis[j]))
            .collect();
        rest.extend(new);
        self.basis = rest;
    }
}

/// The rows of constraint 0 of `eqs` (a view holding just the candidate)
/// over the free parameters.
pub(crate) fn rows_of<E: Equations + ?Sized>(
    eqs: &E,
    params: &[f64],
    fixed: &[bool],
) -> Vec<SparseRow> {
    let mut rows: Vec<SparseRow> = vec![Vec::new(); eqs.equation_count(0)];
    eqs.jacobian(0, params, &mut rows);
    rows.into_iter()
        .map(|row| {
            let mut acc: HashMap<usize, f64> = HashMap::new();
            for (p, v) in row {
                if !fixed[p] {
                    *acc.entry(p).or_insert(0.0) += v;
                }
            }
            to_sparse(acc)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space3d::diagnosis::SpaceGauge;
    use crate::space3d::solver::SpaceEquations;
    use crate::space3d::SpaceSketch;

    /// The basis's dimension against the diagnosis's dof.
    fn dimension_and_dof(sk: &SpaceSketch) -> (usize, i64) {
        let eqs = SpaceEquations::new(&sk.entities, &sk.constraints);
        let fixed = vec![false; sk.params.len()];
        let gauge = SpaceGauge {
            entities: &sk.entities,
        };
        let space = NullBasis::of(&eqs, &sk.params, &fixed, &gauge);
        (space.dimension(), sk.diagnose().dof)
    }

    #[test]
    fn the_null_space_counts_what_the_diagnosis_counts_under_a_gauge() {
        // A free cylinder: seven stored numbers, five freedoms (the axis
        // point's slide and the axis length are how it is stored). A free
        // plane: four stored, three freedoms.
        let mut sk = SpaceSketch::new();
        sk.add_cylinder(1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.5);
        assert_eq!(dimension_and_dof(&sk), (5, 5));
        let mut sk = SpaceSketch::new();
        sk.add_plane(0.3, -0.4, 0.5, 1.2);
        assert_eq!(dimension_and_dof(&sk), (3, 3));
        // A cylinder through a held point: still the slide and the length
        // are not freedoms.
        let mut sk = SpaceSketch::new();
        let c = sk.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
        let p = sk.add_point(1.0, 0.0, 0.0);
        sk.add_constraint(crate::space3d::SpaceConstraint::fix_point(p, 1.0, 0.0, 0.0));
        sk.add_constraint(crate::space3d::SpaceConstraint::PointOnCylinder {
            point: crate::space3d::SpacePointRef::Point(p),
            cylinder: c,
        });
        assert_eq!(dimension_and_dof(&sk), (4, 4));
    }
}
