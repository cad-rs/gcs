//! Elimination of the trivial equations of a component before it is solved.
//!
//! Many equations only say that one parameter equals another plus a
//! constant (two points coincide, a line is horizontal, two circles share a
//! radius) or that a parameter has a value (a fixed point, a radius). Each
//! such equation lets one unknown be written in terms of another, so the
//! parameters fall into classes whose members differ by known offsets: a
//! union–find with offsets. A class reached by a value, or by a held
//! parameter, is known outright. The solve then has one column per class
//! that is still unknown, instead of one per parameter.
//!
//! An equation that would join a class to itself (a cycle of such
//! equations) or two known classes is not used to eliminate anything: it
//! stays an ordinary equation of the solve, which then has to meet it like
//! any other — a consistent duplicate holds to rounding, a contradiction
//! leaves a residual nothing can move. So nothing an equation says is ever
//! dropped, and the diagnosis, which reads the full system, still names
//! every constraint.

use crate::numerics::StepSpace;
use crate::system::{Alias, Equations};

/// No class: a parameter the component does not touch.
const NONE: usize = usize::MAX;

/// The classes of the parameters one component touches.
pub(crate) struct Aliases {
    /// The parameters the component touches, held or free, ascending; the
    /// local index of a parameter is its position here.
    touched: Vec<usize>,
    /// Union–find over the local indices: `value(i) = value(parent[i]) +
    /// offset[i]`.
    parent: Vec<usize>,
    offset: Vec<f64>,
    /// Per root: the value of the root, once the class is known.
    known: Vec<Option<f64>>,
    /// Per equation row of the component, in order: whether it was used to
    /// eliminate an unknown (and so holds by construction).
    absorbed: Vec<bool>,
    /// How many rows were absorbed.
    count: usize,
}

impl Aliases {
    /// The classes of the parameters `touched` (ascending) that the
    /// constraints `indices` imply, the parameters `fixed` holds being
    /// known at their values in `params` — or `None` if no equation can be
    /// absorbed, which is found without allocating when no equation is an
    /// alias at all (every family but the 2D one).
    pub(crate) fn classes_of<E: Equations + ?Sized>(
        eqs: &E,
        indices: &[usize],
        touched: &[usize],
        params: &[f64],
        fixed: &[bool],
    ) -> Option<Classes> {
        let some_alias = indices
            .iter()
            .any(|&i| (0..eqs.equation_count(i)).any(|r| eqs.alias(i, r).is_some()));
        if !some_alias {
            return None;
        }
        let aliases = Self::new(eqs, indices, touched.to_vec(), params, fixed);
        (aliases.count > 0).then(|| aliases.classes())
    }

    fn new<E: Equations + ?Sized>(
        eqs: &E,
        indices: &[usize],
        touched: Vec<usize>,
        params: &[f64],
        fixed: &[bool],
    ) -> Self {
        let n = touched.len();
        let known = touched
            .iter()
            .map(|&p| fixed[p].then_some(params[p]))
            .collect();
        let mut aliases = Self {
            touched,
            parent: (0..n).collect(),
            offset: vec![0.0; n],
            known,
            absorbed: Vec::new(),
            count: 0,
        };
        for &i in indices {
            for r in 0..eqs.equation_count(i) {
                let absorbed = match eqs.alias(i, r) {
                    Some(Alias::Offset {
                        param,
                        base,
                        offset,
                    }) => aliases.relate(param, base, offset),
                    Some(Alias::Value { param, value }) => aliases.pin(param, value),
                    None => false,
                };
                aliases.count += usize::from(absorbed);
                aliases.absorbed.push(absorbed);
            }
        }
        aliases
    }

    /// The local index of parameter `param`, which the component touches.
    fn local(&self, param: usize) -> usize {
        self.touched
            .binary_search(&param)
            .expect("an alias names a parameter its constraint lists")
    }

    /// The root of local index `i` and `value(i) − value(root)`.
    fn find(&mut self, i: usize) -> (usize, f64) {
        let mut root = i;
        let mut total = 0.0;
        while self.parent[root] != root {
            total += self.offset[root];
            root = self.parent[root];
        }
        // Point the path at the root, each with its offset to it.
        let (mut at, mut rest) = (i, total);
        while self.parent[at] != root && at != root {
            let next = self.parent[at];
            let step = self.offset[at];
            self.parent[at] = root;
            self.offset[at] = rest;
            rest -= step;
            at = next;
        }
        (root, total)
    }

    /// Absorb `value(param) − value(base) = offset`, if it joins two
    /// classes not both known.
    fn relate(&mut self, param: usize, base: usize, offset: f64) -> bool {
        let (a, off_a) = self.find(self.local(param));
        let (b, off_b) = self.find(self.local(base));
        if a == b || (self.known[a].is_some() && self.known[b].is_some()) {
            return false;
        }
        // value(a) = value(b) + link.
        let link = offset + off_b - off_a;
        // The smaller index becomes the root, so the classes do not depend
        // on the order the equations come in beyond which rows join them.
        if a < b {
            self.parent[b] = a;
            self.offset[b] = -link;
            if let Some(vb) = self.known[b] {
                self.known[a] = Some(vb + link);
            }
        } else {
            self.parent[a] = b;
            self.offset[a] = link;
            if let Some(va) = self.known[a] {
                self.known[b] = Some(va - link);
            }
        }
        true
    }

    /// Absorb `value(param) = value`, if its class is not known yet.
    fn pin(&mut self, param: usize, value: f64) -> bool {
        let (root, off) = self.find(self.local(param));
        if self.known[root].is_some() {
            return false;
        }
        self.known[root] = Some(value - off);
        true
    }

    /// Settle the classes: every parameter's root and offset, final.
    fn classes(mut self) -> Classes {
        let n = self.touched.len();
        let mut root = vec![NONE; n];
        let mut offset = vec![0.0; n];
        for i in 0..n {
            let (r, off) = self.find(i);
            root[i] = r;
            offset[i] = off;
        }
        Classes {
            touched: self.touched,
            root,
            offset,
            known: self.known,
            absorbed: self.absorbed,
        }
    }
}

/// The settled classes of [`Aliases`].
pub(crate) struct Classes {
    pub(crate) touched: Vec<usize>,
    /// Per local index: its root.
    root: Vec<usize>,
    /// Per local index: `value(i) − value(root)`.
    offset: Vec<f64>,
    known: Vec<Option<f64>>,
    pub(crate) absorbed: Vec<bool>,
}

impl Classes {
    /// Whether the class of parameter `param` is known.
    pub(crate) fn is_known(&self, param: usize) -> bool {
        let i = self.touched.binary_search(&param).expect("touched");
        self.known[self.root[i]].is_some()
    }

    /// The root parameter of the class of `param`.
    pub(crate) fn root_of(&self, param: usize) -> usize {
        let i = self.touched.binary_search(&param).expect("touched");
        self.touched[self.root[i]]
    }

    /// Put every free parameter where its class says: a known class at its
    /// value, an unknown one at the value that moves its members least (the
    /// mean of what each member says the root is: the nearest point where
    /// the absorbed equations hold), each member at the root plus its
    /// offset. Held parameters are not written.
    pub(crate) fn snap(&self, params: &mut [f64], fixed: &[bool]) {
        let n = self.touched.len();
        let mut sum = vec![0.0_f64; n];
        let mut members = vec![0usize; n];
        for i in 0..n {
            let r = self.root[i];
            sum[r] += params[self.touched[i]] - self.offset[i];
            members[r] += 1;
        }
        for i in 0..n {
            let p = self.touched[i];
            if fixed[p] {
                continue;
            }
            let r = self.root[i];
            let base = self.known[r].unwrap_or(sum[r] / members[r] as f64);
            params[p] = base + self.offset[i];
        }
    }
}

/// The step space of a component solved over its classes: column `k` moves
/// every member of one unknown class by the same amount, so the offsets
/// within the class are kept.
pub(crate) struct ClassSpace {
    /// The members that steps move, ascending.
    storage: Vec<usize>,
    /// Per member, its column.
    column: Vec<usize>,
    /// Number of columns.
    dim: usize,
}

impl ClassSpace {
    /// The space of the columns `roots` (root parameters, ascending), whose
    /// members are the touched free parameters of `classes` whose root is
    /// among them.
    pub(crate) fn new(classes: &Classes, roots: &[usize], fixed: &[bool]) -> Self {
        let mut storage = Vec::new();
        let mut column = Vec::new();
        for &p in &classes.touched {
            if fixed[p] || classes.is_known(p) {
                continue;
            }
            if let Ok(k) = roots.binary_search(&classes.root_of(p)) {
                storage.push(p);
                column.push(k);
            }
        }
        Self {
            storage,
            column,
            dim: roots.len(),
        }
    }

    /// [`ClassSpace::new`] for a few of the classes, found through
    /// `members` — every free member of an unknown class as `(root,
    /// member)`, sorted — rather than by going through every parameter.
    pub(crate) fn of_roots(members: &[(usize, usize)], roots: &[usize]) -> Self {
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        for (k, &root) in roots.iter().enumerate() {
            let from = members.partition_point(|&(r, _)| r < root);
            let to = members.partition_point(|&(r, _)| r <= root);
            pairs.extend(members[from..to].iter().map(|&(_, p)| (p, k)));
        }
        pairs.sort_unstable();
        Self {
            storage: pairs.iter().map(|&(p, _)| p).collect(),
            column: pairs.iter().map(|&(_, k)| k).collect(),
            dim: roots.len(),
        }
    }

    /// Every free member of an unknown class of `classes`, as `(root,
    /// member)`, sorted (for [`ClassSpace::of_roots`]).
    pub(crate) fn members(classes: &Classes, fixed: &[bool]) -> Vec<(usize, usize)> {
        let mut members: Vec<(usize, usize)> = classes
            .touched
            .iter()
            .filter(|&&p| !fixed[p] && !classes.is_known(p))
            .map(|&p| (classes.root_of(p), p))
            .collect();
        members.sort_unstable();
        members
    }

    /// Per member of [`StepSpace::storage`], its column.
    pub(crate) fn columns(&self) -> &[usize] {
        &self.column
    }
}

impl StepSpace for ClassSpace {
    fn dim(&self) -> usize {
        self.dim
    }

    fn storage(&self) -> &[usize] {
        &self.storage
    }

    fn apply(&self, params: &mut [f64], saved: &[f64], step: &[f64]) {
        for ((&p, &k), s) in self.storage.iter().zip(&self.column).zip(saved) {
            params[p] = s + step[k];
        }
    }
}
