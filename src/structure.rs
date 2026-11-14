//! The structure of a component's Jacobian: which equations can be solved
//! for which unknowns, and in what order.
//!
//! A maximum matching of equations to unknowns (Hopcroft–Karp) splits the
//! system into the Dulmage–Mendelsohn parts: an over-determined part (more
//! equations than unknowns, every unknown of it tied to its equations
//! alone), an under-determined one (more unknowns than equations), and a
//! square part between them. The square part falls apart further into its
//! strongly connected components (Tarjan), each a block of equations that
//! has to be solved together for as many unknowns, in an order where every
//! block only uses unknowns of itself and of blocks before it. The
//! over-determined part only uses its own unknowns and comes first; the
//! under-determined part may use any and comes last.
//!
//! All of this reads the pattern alone, not the values: a block that is
//! structurally square may still be singular where the solve stands.

use crate::numerics::SparseJacobian;

const NONE: usize = usize::MAX;

/// A set of equations solved together for a set of unknowns.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Block {
    /// Rows of the pattern, ascending.
    pub(crate) rows: Vec<usize>,
    /// Columns of the pattern, ascending.
    pub(crate) cols: Vec<usize>,
}

/// The blocks of `pattern` (with `n_cols` columns) in the order they can be
/// solved, leaving out the rows `skip` marks. Rows with no column (other
/// than skipped ones) go with the over-determined part: nothing can move
/// them, so they are checked first.
pub(crate) fn blocks(pattern: &SparseJacobian, n_cols: usize, skip: &[bool]) -> Vec<Block> {
    let n_rows = pattern.n_rows();
    let row_cols = |r: usize| pattern.row_cols(r);
    let live: Vec<usize> = (0..n_rows).filter(|&r| !skip[r]).collect();
    let (col_of, row_of) = matching(pattern, &live, n_cols);

    // The rows touching each column (live rows only).
    let mut starts = vec![0usize; n_cols + 1];
    for &r in &live {
        for &c in row_cols(r) {
            starts[c + 1] += 1;
        }
    }
    for c in 0..n_cols {
        starts[c + 1] += starts[c];
    }
    let mut at = starts.clone();
    let mut col_rows = vec![0usize; starts[n_cols]];
    for &r in &live {
        for &c in row_cols(r) {
            col_rows[at[c]] = r;
            at[c] += 1;
        }
    }

    // Under-determined: what an unmatched column reaches, column → any row
    // touching it → that row's matched column.
    let mut under_c = vec![false; n_cols];
    let mut under_r = vec![false; n_rows];
    let mut stack: Vec<usize> = (0..n_cols).filter(|&c| row_of[c] == NONE).collect();
    for &c in &stack {
        under_c[c] = true;
    }
    while let Some(c) = stack.pop() {
        for &r in &col_rows[starts[c]..starts[c + 1]] {
            if !under_r[r] {
                under_r[r] = true;
                let next = col_of[r];
                if next != NONE && !under_c[next] {
                    under_c[next] = true;
                    stack.push(next);
                }
            }
        }
    }
    // Over-determined: what an unmatched row reaches, row → any column it
    // touches → that column's matched row.
    let mut over_r = vec![false; n_rows];
    let mut over_c = vec![false; n_cols];
    let mut stack: Vec<usize> = live
        .iter()
        .copied()
        .filter(|&r| col_of[r] == NONE)
        .collect();
    for &r in &stack {
        over_r[r] = true;
    }
    while let Some(r) = stack.pop() {
        for &c in row_cols(r) {
            if !over_c[c] {
                over_c[c] = true;
                let next = row_of[c];
                if next != NONE && !over_r[next] {
                    over_r[next] = true;
                    stack.push(next);
                }
            }
        }
    }

    let mut out = Vec::new();
    let over = Block {
        rows: live.iter().copied().filter(|&r| over_r[r]).collect(),
        cols: (0..n_cols).filter(|&c| over_c[c]).collect(),
    };
    if !over.rows.is_empty() {
        out.push(over);
    }

    // The square part: a graph on its rows, `r → r'` when row `r` uses the
    // column matched to `r'`. Tarjan emits a strongly connected component
    // once every component it reaches has been emitted, which is the order
    // they can be solved in.
    let square: Vec<usize> = live
        .iter()
        .copied()
        .filter(|&r| !under_r[r] && !over_r[r])
        .collect();
    let mut local = vec![NONE; n_rows];
    for (k, &r) in square.iter().enumerate() {
        local[r] = k;
    }
    let (row_of, local, square) = (&row_of, &local, &square);
    let successors = |k: usize| {
        let r = square[k];
        row_cols(r).iter().filter_map(move |&c| {
            let to = row_of[c];
            (to != NONE && to != r && local[to] != NONE).then(|| local[to])
        })
    };
    for component in strongly_connected(square.len(), successors) {
        let mut rows: Vec<usize> = component.iter().map(|&k| square[k]).collect();
        rows.sort_unstable();
        let mut cols: Vec<usize> = rows.iter().map(|&r| col_of[r]).collect();
        cols.sort_unstable();
        out.push(Block { rows, cols });
    }

    let under = Block {
        rows: live.iter().copied().filter(|&r| under_r[r]).collect(),
        cols: (0..n_cols).filter(|&c| under_c[c]).collect(),
    };
    if !under.cols.is_empty() {
        out.push(under);
    }
    out
}

/// A maximum matching of the rows `live` of `pattern` to its `n_cols`
/// columns (Hopcroft–Karp): per row its column, per column its row, or
/// `NONE`.
fn matching(pattern: &SparseJacobian, live: &[usize], n_cols: usize) -> (Vec<usize>, Vec<usize>) {
    let n_rows = pattern.n_rows();
    let mut col_of = vec![NONE; n_rows];
    let mut row_of = vec![NONE; n_cols];
    // A greedy start matches most rows of a sparse system at once.
    for &r in live {
        if let Some(&c) = pattern.row_cols(r).iter().find(|&&c| row_of[c] == NONE) {
            col_of[r] = c;
            row_of[c] = r;
        }
    }
    let mut level = vec![NONE; n_rows];
    let mut queue: Vec<usize> = Vec::with_capacity(live.len());
    loop {
        // Breadth first from the free rows, along unmatched edges to a
        // column and back along its matching edge.
        queue.clear();
        for &r in live {
            if col_of[r] == NONE {
                level[r] = 0;
                queue.push(r);
            } else {
                level[r] = NONE;
            }
        }
        let mut found = false;
        let mut head = 0;
        while head < queue.len() {
            let r = queue[head];
            head += 1;
            for &c in pattern.row_cols(r) {
                let next = row_of[c];
                if next == NONE {
                    found = true;
                } else if level[next] == NONE {
                    level[next] = level[r] + 1;
                    queue.push(next);
                }
            }
        }
        if !found {
            break;
        }
        // Depth first along the levels, iteratively, one augmenting path
        // per free row at most.
        let mut edge = vec![0usize; n_rows];
        let mut path: Vec<usize> = Vec::new();
        for &start in live {
            if col_of[start] != NONE {
                continue;
            }
            path.clear();
            path.push(start);
            while let Some(&r) = path.last() {
                let cols = pattern.row_cols(r);
                if edge[r] == cols.len() {
                    // A dead end: no path through `r` in this phase.
                    level[r] = NONE;
                    path.pop();
                    continue;
                }
                let c = cols[edge[r]];
                edge[r] += 1;
                let next = row_of[c];
                if next == NONE {
                    // Augment along the path: each row takes the column it
                    // went through.
                    for &row in path.iter().rev() {
                        let through = pattern.row_cols(row)[edge[row] - 1];
                        col_of[row] = through;
                        row_of[through] = row;
                    }
                    break;
                }
                if level[next] != NONE && level[next] == level[r] + 1 {
                    path.push(next);
                }
            }
        }
    }
    (col_of, row_of)
}

/// The strongly connected components of a graph on `n` nodes, in the order
/// Tarjan's algorithm emits them (every component after all the components
/// it reaches). Iterative, so a long chain does not exhaust the stack.
fn strongly_connected<I, F>(n: usize, successors: F) -> Vec<Vec<usize>>
where
    F: Fn(usize) -> I,
    I: Iterator<Item = usize>,
{
    let mut index = vec![NONE; n];
    let mut low = vec![0usize; n];
    let mut on_stack = vec![false; n];
    let mut stack: Vec<usize> = Vec::new();
    let mut out = Vec::new();
    let mut counter = 0;
    for root in 0..n {
        if index[root] != NONE {
            continue;
        }
        let mut calls: Vec<(usize, I)> = Vec::new();
        index[root] = counter;
        low[root] = counter;
        counter += 1;
        stack.push(root);
        on_stack[root] = true;
        calls.push((root, successors(root)));
        while let Some((v, edges)) = calls.last_mut() {
            let v = *v;
            match edges.next() {
                Some(w) if index[w] == NONE => {
                    index[w] = counter;
                    low[w] = counter;
                    counter += 1;
                    stack.push(w);
                    on_stack[w] = true;
                    calls.push((w, successors(w)));
                }
                Some(w) => {
                    if on_stack[w] {
                        low[v] = low[v].min(index[w]);
                    }
                }
                None => {
                    calls.pop();
                    if let Some((u, _)) = calls.last() {
                        low[*u] = low[*u].min(low[v]);
                    }
                    if low[v] == index[v] {
                        let mut component = Vec::new();
                        loop {
                            let w = stack.pop().expect("on the stack");
                            on_stack[w] = false;
                            component.push(w);
                            if w == v {
                                break;
                            }
                        }
                        out.push(component);
                    }
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pattern(n_cols: usize, rows: &[&[usize]]) -> SparseJacobian {
        SparseJacobian::from_rows(n_cols, rows)
    }

    #[test]
    fn a_lower_triangle_falls_into_single_blocks_in_order() {
        // x0 = .., x1(x0), x2(x1, x0): three blocks, x0 first.
        let p = pattern(3, &[&[1, 2], &[0], &[0, 1]]);
        let b = blocks(&p, 3, &[false; 3]);
        let order: Vec<(Vec<usize>, Vec<usize>)> =
            b.into_iter().map(|b| (b.rows, b.cols)).collect();
        assert_eq!(
            order,
            vec![(vec![1], vec![0]), (vec![2], vec![1]), (vec![0], vec![2])]
        );
    }

    #[test]
    fn a_coupled_pair_is_one_block() {
        let p = pattern(3, &[&[0, 1], &[0, 1], &[1, 2]]);
        let b = blocks(&p, 3, &[false; 3]);
        assert_eq!(b.len(), 2);
        assert_eq!(
            (b[0].rows.clone(), b[0].cols.clone()),
            (vec![0, 1], vec![0, 1])
        );
        assert_eq!((b[1].rows.clone(), b[1].cols.clone()), (vec![2], vec![2]));
    }

    #[test]
    fn over_and_under_parts_come_first_and_last() {
        // Rows 0, 1 both only on column 0 (over); row 2 on columns 0, 1, 2
        // (under: two unknowns for one equation); row 3 on column 3 alone.
        let p = pattern(4, &[&[0], &[0], &[0, 1, 2], &[3]]);
        let b = blocks(&p, 4, &[false; 4]);
        assert_eq!(
            (b[0].rows.clone(), b[0].cols.clone()),
            (vec![0, 1], vec![0])
        );
        assert_eq!((b[1].rows.clone(), b[1].cols.clone()), (vec![3], vec![3]));
        let last = b.last().unwrap();
        assert_eq!(
            (last.rows.clone(), last.cols.clone()),
            (vec![2], vec![1, 2])
        );
    }

    #[test]
    fn skipped_rows_are_left_out_and_empty_rows_are_checked_first() {
        let p = pattern(2, &[&[0, 1], &[], &[1], &[]]);
        let b = blocks(&p, 2, &[false, false, false, true]);
        assert_eq!((b[0].rows.clone(), b[0].cols.clone()), (vec![1], vec![]));
        let rows: Vec<usize> = b.iter().flat_map(|b| b.rows.clone()).collect();
        assert!(!rows.contains(&3));
        assert_eq!(b.len(), 3);
    }

    #[test]
    fn the_matching_is_maximum_where_greedy_is_not() {
        // Greedy gives row 0 column 0, leaving row 1 (only column 0)
        // unmatched; augmenting moves row 0 to column 1.
        let p = pattern(2, &[&[0, 1], &[0]]);
        let (col_of, _) = matching(&p, &[0, 1], 2);
        assert_eq!(col_of, vec![1, 0]);
    }

    #[test]
    fn a_long_chain_does_not_overflow() {
        let n = 200_000;
        let rows: Vec<Vec<usize>> = (0..n)
            .map(|i| if i == 0 { vec![0] } else { vec![i - 1, i] })
            .collect();
        let p = SparseJacobian::from_rows(n, &rows);
        let b = blocks(&p, n, &vec![false; n]);
        assert_eq!(b.len(), n);
        assert_eq!(b[0].cols, vec![0]);
        assert_eq!(b[n - 1].cols, vec![n - 1]);
    }
}
