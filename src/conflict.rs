//! Which constraints conflict, as small sets a host can act on.
//!
//! [`crate::Diagnosis`] names conflicting constraints by linear algebra at
//! the current geometry: a dependency among the rows whose values
//! disagree. That misses a conflict that is not linear (two distances that
//! cannot bridge two fixed points: the rows are independent wherever the
//! geometry stands) and names one representative, not what to do.
//! [`Conflict`] is found by solving instead: for each connected component
//! that does not solve from where it stands, a deletion filter takes
//! constraints out one at a time, newest first, and leaves out each one
//! the component still fails without; what is left fails, and taking out
//! any one of it lets the rest solve — a minimal conflicting set. Its
//! newest member is the suggestion; with it taken out the filter runs
//! again, until the component solves. A conflicting dependency group of
//! the diagnosis, when it fails on its own, is filtered instead of the
//! whole component (fewer solves).
//!
//! "Fails" means a solve from the current geometry does not converge, so
//! like any solve the answer depends on where the geometry stands: a set
//! that could be met far away but not reached from here counts as
//! conflicting. Each test is one solve of the constraints left on; the
//! host's geometry is not touched.

/// A minimal set of constraints that cannot all hold, as far as a solve
/// from the current geometry can tell: with all of them on the component
/// does not solve, with any one of them off it does.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Conflict {
    /// The constraint indices, newest (highest index, the order they were
    /// added in) first. The first is the one to suggest taking out.
    pub constraints: Vec<usize>,
}

impl Conflict {
    /// The constraint to suggest taking out: the newest of the set.
    pub fn suggested(&self) -> usize {
        self.constraints[0]
    }
}

/// The conflicts of a system whose active constraints fall into
/// `components`, where `fails(on)` says whether a solve with only the
/// constraints `on` active (ascending) fails from the current geometry,
/// and `groups` are the conflicting dependency groups of the diagnosis
/// there (a hint).
pub(crate) fn find_conflicts(
    components: &[Vec<usize>],
    groups: &[Vec<usize>],
    fails: &dyn Fn(&[usize]) -> bool,
) -> Vec<Conflict> {
    let mut out = Vec::new();
    for component in components {
        let mut active: Vec<usize> = component.clone();
        active.sort_unstable();
        if !fails(&active) {
            continue;
        }
        // At most one conflict per constraint: each round takes one out.
        for _ in 0..component.len() {
            let start = groups
                .iter()
                .map(|g| {
                    let mut g = g.clone();
                    g.sort_unstable();
                    g
                })
                .find(|g| g.iter().all(|c| active.binary_search(c).is_ok()) && fails(g))
                .unwrap_or_else(|| active.clone());
            let set = deletion_filter(start, fails);
            let suggested = *set.last().expect("a failing set is not empty");
            out.push(Conflict {
                constraints: set.into_iter().rev().collect(),
            });
            active.retain(|&c| c != suggested);
            if !fails(&active) {
                break;
            }
        }
    }
    out
}

/// The dependency groups of `diagnosis` that hold a conflicting
/// constraint.
pub(crate) fn conflicting_groups(diagnosis: &crate::Diagnosis) -> Vec<Vec<usize>> {
    diagnosis
        .dependency_groups
        .iter()
        .filter(|g| g.iter().any(|c| diagnosis.conflicting.contains(c)))
        .cloned()
        .collect()
}

/// A minimal failing subset of the failing set `set` (ascending): each
/// member, newest first, is left out when the rest still fails.
fn deletion_filter(mut set: Vec<usize>, fails: &dyn Fn(&[usize]) -> bool) -> Vec<usize> {
    for c in set.clone().into_iter().rev() {
        let without: Vec<usize> = set.iter().copied().filter(|&k| k != c).collect();
        if fails(&without) {
            set = without;
        }
    }
    set
}

#[cfg(test)]
mod tests {
    use super::*;

    /// "Fails" when one of `sets` is all on.
    fn failing_if_any(sets: Vec<Vec<usize>>) -> impl Fn(&[usize]) -> bool {
        move |on: &[usize]| sets.iter().any(|s| s.iter().all(|c| on.contains(c)))
    }

    #[test]
    fn the_filter_keeps_the_newest_and_the_oldest_it_needs() {
        // {0, 2} and {1, 2} both fail; the newest-first filter drops 1.
        let fails = failing_if_any(vec![vec![0, 2], vec![1, 2]]);
        assert_eq!(deletion_filter(vec![0, 1, 2], &fails), vec![0, 2]);
    }

    #[test]
    fn each_round_takes_out_one_suggestion_until_the_component_solves() {
        // Two separate conflicts in one component, and a constraint (3)
        // in neither. Without 4 the component still fails (on {1, 2}), so
        // that conflict is found first.
        let fails = failing_if_any(vec![vec![0, 4], vec![1, 2]]);
        let found = find_conflicts(&[vec![0, 1, 2, 3, 4]], &[], &fails);
        assert_eq!(
            found,
            vec![
                Conflict {
                    constraints: vec![2, 1]
                },
                Conflict {
                    constraints: vec![4, 0]
                },
            ]
        );
        assert_eq!(found[0].suggested(), 2);
    }

    #[test]
    fn a_failing_group_is_filtered_instead_of_the_component() {
        let fails = failing_if_any(vec![vec![5, 7]]);
        let calls = std::cell::Cell::new(0);
        let counting = |on: &[usize]| {
            calls.set(calls.get() + 1);
            fails(on)
        };
        let found = find_conflicts(&[(0..10).collect()], &[vec![7, 5, 6]], &counting);
        assert_eq!(
            found,
            vec![Conflict {
                constraints: vec![7, 5]
            }]
        );
        // Whole component, group, three filter steps, the rest.
        assert_eq!(calls.get(), 6);
    }

    #[test]
    fn a_component_that_solves_has_no_conflict() {
        let fails = failing_if_any(vec![vec![9]]);
        assert!(find_conflicts(&[vec![0, 1]], &[], &fails).is_empty());
    }
}
