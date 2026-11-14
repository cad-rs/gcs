//! What one more constraint would do, without diagnosing the whole system
//! again.
//!
//! A what-if probe (`Sketch::what_if`, `SpaceSketch::what_if`,
//! `rigid3d::assembly_what_if`) keeps, per connected component, the row space of
//! the equilibrated constraint Jacobian at the current parameters as the
//! diagnosis saw it: an orthonormal basis of the independent rows, their
//! least-squares step, and how far each dependent row is from it. It is
//! built for a component the first time a candidate touches it. Asked
//! about a candidate, it projects the candidate's rows onto that space: a
//! row with something left over (beyond the rank tolerance of the
//! diagnosis) is independent and extends the space; a row without is
//! dependent, and is classified the way the diagnosis classifies a
//! dependent row — by how far its residual is from the least-squares step
//! of the independent rows. The verdict is the one [`crate::Diagnosis`]
//! gives with the candidate added: [`Verdict::Independent`] when the rank
//! grows by its number of equations; otherwise [`Verdict::Conflicting`]
//! when the diagnosis would name a constraint conflicting that it does not
//! name now, [`Verdict::PartiallyRedundant`] likewise for the gray band,
//! and [`Verdict::Redundant`] when neither.
//!
//! Candidates that change more than rows — held columns that differ, a
//! plane whose effective columns change because the candidate reaches a
//! parameter of it the component did not, a component with nothing free —
//! are diagnosed exactly on the components they touch instead. Either way the system itself is not
//! touched. Like the diagnosis, the answer is at the current parameters:
//! probe a solved system.

use std::cell::OnceCell;
use std::collections::HashMap;

use crate::diagnosis::{CONSIST_TOL, PARTIAL_TOL, RANK_REL_TOL};
use crate::system::{self, ComponentDiagnosis, Equations, Gauge, SparseRow};

/// What adding a constraint would make of it (see the
/// [module docs](self)).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Verdict {
    /// Every equation of it is independent of the system: it takes away as
    /// many freedoms as it has equations.
    Independent,
    /// Some equation of it depends on the system and agrees with it.
    Redundant,
    /// Some equation of it depends on the system and only nearly agrees
    /// (the gray band of [`crate::Diagnosis::partially_redundant`]).
    PartiallyRedundant,
    /// Some equation of it depends on the system and disagrees: the
    /// diagnosis would name a conflict it does not name now.
    Conflicting,
}

/// The answer of a what-if probe about one candidate constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Addition {
    pub verdict: Verdict,
    /// How many degrees of freedom adding it would take away (the drop of
    /// [`crate::Diagnosis::dof`]), from none up to its number of
    /// equations.
    pub removes: usize,
}

/// The cached row spaces of a system, component by component.
pub(crate) struct RowSpaces<E, G> {
    eqs: E,
    params: Vec<f64>,
    fixed: Vec<bool>,
    gauge: G,
    components: Vec<Vec<usize>>,
    /// Per column, the component whose rows list it.
    component_of: Vec<usize>,
    parts: Vec<OnceCell<ComponentDiagnosis>>,
    soft: Vec<SparseRow>,
    /// Per soft gauge direction, whether some component touches it.
    soft_touched: Vec<bool>,
}

/// A row whose part outside the row space is within this factor of the
/// rank tolerance either way could fall on either side of the diagnosis's
/// decision: the exact test decides it.
const RANK_MARGIN: f64 = 1.0e3;
/// Below this smallest relative pivot the independent rows are too poorly
/// conditioned for the step through them to classify a dependent row the
/// way a re-pivoted factorisation would: the exact test decides.
const WELL_CONDITIONED: f64 = 1.0e-3;
/// A dependent row within this factor of a classification threshold
/// could be classified either way once the factorisation re-pivots (its
/// distance scales with the coefficient of the row it swaps with): the
/// exact test decides.
const CLASS_MARGIN: f64 = 1.0e2;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// The key of an effective column: its global parameters and weights,
/// bit for bit.
fn key(column: &SparseRow) -> Vec<(usize, u64)> {
    column.iter().map(|&(p, w)| (p, w.to_bits())).collect()
}

impl<E: Equations, G: Gauge> RowSpaces<E, G> {
    pub(crate) fn new(eqs: E, params: Vec<f64>, fixed: Vec<bool>, gauge: G) -> Self {
        let components = system::components(&eqs, &params);
        let mut component_of = vec![usize::MAX; eqs.column_count(&params)];
        for (k, group) in components.iter().enumerate() {
            for p in system::touched_params(&eqs, group, &params) {
                component_of[p] = k;
            }
        }
        let soft = gauge.soft_directions(&params, &fixed);
        let soft_touched = soft
            .iter()
            .map(|dir| {
                dir.iter()
                    .any(|&(p, _)| !fixed[p] && component_of[p] != usize::MAX)
            })
            .collect();
        let parts = components.iter().map(|_| OnceCell::new()).collect();
        Self {
            eqs,
            params,
            fixed,
            gauge,
            components,
            component_of,
            parts,
            soft,
            soft_touched,
        }
    }

    /// The diagnosis of component `k`, with its row space kept.
    fn part(&self, k: usize) -> &ComponentDiagnosis {
        self.parts[k].get_or_init(|| {
            let group = &self.components[k];
            let f_c = system::eval_residuals(&self.eqs, group, &self.params);
            let mut seen = vec![false; self.soft.len()];
            system::diagnose_component(
                &self.eqs,
                &self.params,
                &self.fixed,
                &self.gauge,
                group,
                f_c,
                &self.soft,
                &mut seen,
                true,
            )
        })
    }

    /// What adding constraint `candidate` of `with` — this system's
    /// equations with the candidate appended, over the same columns, with
    /// the held columns `fixed_with` — would do.
    pub(crate) fn would_add<W: Equations + ?Sized>(
        &self,
        with: &W,
        fixed_with: &[bool],
        candidate: usize,
    ) -> Addition {
        let k = with.equation_count(candidate);
        if k == 0 {
            return Addition {
                verdict: Verdict::Independent,
                removes: 0,
            };
        }
        let listed = system::touched_params(with, &[candidate], &self.params);
        let mut merged: Vec<usize> = listed
            .iter()
            .map(|&p| self.component_of[p])
            .filter(|&c| c != usize::MAX)
            .collect();
        merged.sort_unstable();
        merged.dedup();
        self.fast(with, fixed_with, candidate, &merged, &listed)
            .unwrap_or_else(|| self.exact(with, fixed_with, candidate, &merged))
    }

    /// The incremental test; `None` when it does not apply.
    fn fast<W: Equations + ?Sized>(
        &self,
        with: &W,
        fixed_with: &[bool],
        candidate: usize,
        merged: &[usize],
        listed: &[usize],
    ) -> Option<Addition> {
        if fixed_with != self.fixed.as_slice() || merged.is_empty() {
            return None;
        }
        let parts: Vec<&ComponentDiagnosis> = merged.iter().map(|&c| self.part(c)).collect();
        let bases: Vec<&system::RowBasis> = parts
            .iter()
            .map(|p| p.basis.as_ref())
            .collect::<Option<Vec<_>>>()?;
        // The effective columns of the merged component, and where each
        // component's columns sit among them.
        let mut touched: Vec<usize> = bases
            .iter()
            .flat_map(|b| b.touched.iter().copied())
            .chain(listed.iter().copied().filter(|&p| !self.fixed[p]))
            .collect();
        touched.sort_unstable();
        touched.dedup();
        let columns: Vec<SparseRow> =
            match self
                .gauge
                .effective_columns(&self.params, &self.fixed, &touched)
            {
                None => touched.iter().map(|&p| vec![(p, 1.0)]).collect(),
                Some(local) => local
                    .iter()
                    .map(|col| col.iter().map(|&(lc, w)| (touched[lc], w)).collect())
                    .collect(),
            };
        let index: HashMap<Vec<(usize, u64)>, usize> = columns
            .iter()
            .enumerate()
            .map(|(c, col)| (key(col), c))
            .collect();
        let places: Vec<Vec<usize>> = bases
            .iter()
            .map(|b| {
                b.columns
                    .iter()
                    .map(|col| index.get(&key(col)).copied())
                    .collect()
            })
            .collect::<Option<Vec<_>>>()?;
        let n = columns.len();

        // The candidate's rows over those columns, equilibrated.
        let mut raw: Vec<SparseRow> = vec![Vec::new(); with.equation_count(candidate)];
        with.jacobian(candidate, &self.params, &mut raw);
        let mut f = vec![0.0; raw.len()];
        with.residuals(candidate, &self.params, &mut f);
        let mut rows: Vec<Vec<f64>> = raw
            .iter()
            .map(|row| {
                let value: HashMap<usize, f64> = row.iter().filter(|(p, _)| !self.fixed[*p]).fold(
                    HashMap::new(),
                    |mut m, &(p, v)| {
                        *m.entry(p).or_insert(0.0) += v;
                        m
                    },
                );
                columns
                    .iter()
                    .map(|col| {
                        col.iter()
                            .map(|&(p, w)| value.get(&p).copied().unwrap_or(0.0) * w)
                            .sum()
                    })
                    .collect()
            })
            .collect();
        for (row, fi) in rows.iter_mut().zip(f.iter_mut()) {
            let s = row.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            if s > 0.0 {
                for v in row.iter_mut() {
                    *v /= s;
                }
                *fi /= s;
            }
        }

        // Project onto the components' row spaces and the rows already
        // taken; extend the step for every independent row.
        let mut z = vec![0.0; n];
        for (b, place) in bases.iter().zip(&places) {
            for (c, &at) in place.iter().enumerate() {
                z[at] += b.z[c];
            }
        }
        let max0 = rows
            .iter()
            .map(|r| dot(r, r).sqrt())
            .chain(bases.iter().map(|b| b.max0))
            .fold(0.0_f64, f64::max);
        let mut extra: Vec<Vec<f64>> = Vec::new();
        let mut gain = 0usize;
        let mut dependent: Vec<usize> = Vec::new();
        for (i, row) in rows.iter().enumerate() {
            let mut rest = row.clone();
            for _ in 0..2 {
                for (b, place) in bases.iter().zip(&places) {
                    for q in &b.q {
                        let along: f64 = place.iter().zip(q).map(|(&at, qc)| rest[at] * qc).sum();
                        for (&at, qc) in place.iter().zip(q) {
                            rest[at] -= along * qc;
                        }
                    }
                }
                for u in &extra {
                    let along = dot(&rest, u);
                    for (a, b) in rest.iter_mut().zip(u) {
                        *a -= along * b;
                    }
                }
            }
            let left = dot(&rest, &rest).sqrt();
            let tol = RANK_REL_TOL * max0;
            if left > tol / RANK_MARGIN && left <= tol * RANK_MARGIN {
                return None;
            }
            if left > tol {
                gain += 1;
                let u: Vec<f64> = rest.iter().map(|v| v / left).collect();
                let beta = (-f[i] - dot(row, &z)) / dot(row, &u);
                for (a, b) in z.iter_mut().zip(&u) {
                    *a += beta * b;
                }
                extra.push(u);
            } else {
                dependent.push(i);
            }
        }

        let mut group: Vec<usize> = merged
            .iter()
            .flat_map(|&c| self.components[c].iter().copied())
            .collect();
        group.push(candidate);
        let unit = system::equation_unit(with.length_scale(&group, &self.params));
        let scale = f
            .iter()
            .map(|v| v.abs())
            .chain(bases.iter().map(|b| b.f_max))
            .fold(unit, f64::max);
        // A dependent row is only classified here when the answer cannot
        // hang on which rows a re-pivoted factorisation takes: the
        // components are well conditioned and consistent, and the row is
        // clearly on one side of the thresholds.
        if gain < rows.len()
            && (bases.iter().any(|b| b.min_pivot < WELL_CONDITIONED)
                || parts
                    .iter()
                    .any(|p| !p.conflicting.is_empty() || !p.partially_redundant.is_empty()))
        {
            return None;
        }
        let verdict = if gain == rows.len() {
            Verdict::Independent
        } else {
            // Which rows are conflicting / partially redundant now, and
            // which would be.
            let class = |e: f64| {
                if e <= CONSIST_TOL * scale {
                    0
                } else if e <= PARTIAL_TOL * scale {
                    1
                } else {
                    2
                }
            };
            let mut after: Vec<(usize, u8)> = bases
                .iter()
                .flat_map(|b| b.dependent.iter().map(|&(c, e)| (c, class(e))))
                .collect();
            for &i in &dependent {
                let e = (f[i] + dot(&rows[i], &z)).abs();
                if e > CONSIST_TOL * scale / CLASS_MARGIN && e <= PARTIAL_TOL * scale * CLASS_MARGIN
                {
                    return None;
                }
                after.push((candidate, class(e)));
            }
            let named = |list: fn(&ComponentDiagnosis) -> &Vec<usize>, c: usize| {
                parts.iter().any(|p| list(p).contains(&c))
            };
            let new_in = |level: u8, list: fn(&ComponentDiagnosis) -> &Vec<usize>| {
                after.iter().any(|&(c, l)| l == level && !named(list, c))
            };
            if new_in(2, |p| &p.conflicting) {
                Verdict::Conflicting
            } else if new_in(1, |p| &p.partially_redundant) {
                Verdict::PartiallyRedundant
            } else {
                Verdict::Redundant
            }
        };

        // Soft gauge directions the candidate pins take a rank without
        // taking a freedom.
        let pinned = if self.soft.is_empty() {
            0
        } else {
            let col_of: HashMap<usize, usize> =
                touched.iter().enumerate().map(|(c, &p)| (p, c)).collect();
            let mut j = Vec::new();
            for row in bases.iter().flat_map(|b| &b.j_raw).chain(&raw) {
                let mut dense = vec![0.0; touched.len()];
                for &(p, v) in row {
                    if let Some(&c) = col_of.get(&p) {
                        dense[c] += v;
                    }
                }
                j.push(dense);
            }
            let mut seen = vec![false; self.soft.len()];
            let after = system::unpinned_gauge_dims(&j, &col_of, &self.soft, &mut seen);
            let before = parts.iter().map(|p| p.soft_unpinned).sum::<usize>()
                + seen
                    .iter()
                    .zip(&self.soft_touched)
                    .filter(|(now, before)| **now && !**before)
                    .count();
            before.saturating_sub(after)
        };
        Some(Addition {
            verdict,
            removes: gain.saturating_sub(pinned),
        })
    }

    /// The exact test: the components the candidate touches diagnosed
    /// with and without it.
    fn exact<W: Equations + ?Sized>(
        &self,
        with: &W,
        fixed_with: &[bool],
        candidate: usize,
        merged: &[usize],
    ) -> Addition {
        let mut active = vec![false; self.eqs.constraint_count()];
        for &c in merged {
            for &i in &self.components[c] {
                active[i] = true;
            }
        }
        let before = system::diagnose(
            &crate::status::Only {
                eqs: &self.eqs,
                active: active.clone(),
            },
            &self.params,
            &self.fixed,
            &self.gauge,
        );
        active.resize(with.constraint_count(), false);
        active[candidate] = true;
        let after = system::diagnose(
            &crate::status::Only { eqs: with, active },
            &self.params,
            fixed_with,
            &self.gauge,
        );
        let k = with.equation_count(candidate) as i64;
        let gain = after.rank as i64 - before.rank as i64;
        let grew = |now: &[usize], then: &[usize]| now.iter().any(|c| !then.contains(c));
        let verdict = if gain == k {
            Verdict::Independent
        } else if grew(&after.conflicting, &before.conflicting) {
            Verdict::Conflicting
        } else if grew(&after.partially_redundant, &before.partially_redundant) {
            Verdict::PartiallyRedundant
        } else {
            Verdict::Redundant
        };
        Addition {
            verdict,
            removes: (before.dof - after.dof).max(0) as usize,
        }
    }
}
