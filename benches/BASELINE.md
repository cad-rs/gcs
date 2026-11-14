# Benchmark baseline

Numbers from `cargo bench` (`benches/solve.rs`, criterion, 10 samples each,
median reported). They are a reference for spotting regressions on the same
machine, not a portable claim: re-measure before comparing.

- Machine: Windows 11 Pro, x86_64, `rustc 1.99.1`, bench profile (optimised).
- Date: 2026-11-13.
- "0.1.0" is commit `8094753`; "P0" is that commit plus the baseline fixes
  recorded under *Unreleased* in the changelog. The 0.1.0 column was measured
  by running this same bench file against that commit.

## 2D sketch

Fixtures are rectangles of 18 parameters and 18 equations each. The size in
the benchmark name is the number of parameters.

| Benchmark | What it measures | 0.1.0 | P0 | Speed-up |
|---|---|---:|---:|---:|
| `sketch/solve_chain/54` | cold solve, one connected component | 1.60 ms | 0.364 ms | 4.4× |
| `sketch/solve_chain/198` | | 68.9 ms | 7.54 ms | 9.1× |
| `sketch/solve_chain/792` | | 4.35 s | 357 ms | 12× |
| `sketch/resolve_chain/54` | warm re-solve, 1e-4 off the solution | 1.57 ms | 0.372 ms | 4.2× |
| `sketch/resolve_chain/198` | | 68.2 ms | 7.16 ms | 9.5× |
| `sketch/resolve_chain/792` | | 4.31 s | 362 ms | 12× |
| `sketch/solve_field/360` | cold solve, 20 independent rectangles | 2.11 ms | 1.27 ms | 1.7× |
| `sketch/solve_field/3600` | cold solve, 200 independent rectangles | 13.7 ms | 8.66 ms | 1.6× |
| `sketch/diagnose_chain/54` | rank-revealing diagnosis, one component | 220 µs | 226 µs | — |
| `sketch/diagnose_chain/198` | | 8.31 ms | 8.20 ms | — |
| `sketch/diagnose_chain/792` | | 737 ms | 767 ms | — |

## 3D geometric

A braced truss of points tied together by distances, one connected
component. The size in the benchmark name is the number of parameters. Each
point carries one redundant brace, so the diagnosis has dependencies to
resolve (75 redundant rows at 240 parameters, 255 at 780).

This fixture was added after the 0.1.0 comparison run, so the solve has no
0.1.0 number. The diagnosis "before" column is the same tree without the
change described below, which is the 0.1.0 code path for that function.

| Benchmark | What it measures | before | P0 | Speed-up |
|---|---|---:|---:|---:|
| `space/solve_truss/60` | cold solve, one connected component | — | 321 µs | — |
| `space/solve_truss/240` | | — | 14.3 ms | — |
| `space/solve_truss/780` | | — | 396 ms | — |
| `space/diagnose_truss/60` | diagnosis with redundant rows | 1.38 ms | 435 µs | 3.2× |
| `space/diagnose_truss/240` | | 779 ms | 35.4 ms | 22× |
| `space/diagnose_truss/780` | | about 113 s (criterion's estimate, not run to the end) | 2.43 s | about 47× |

The diagnosis speed-up is one change: the Gram matrix of the independent
rows used to be rebuilt and eliminated afresh for every dependent row; it is
now built and factored once and solved per row.

## Assembly (rigid pose)

A grounded body plus N rigid bodies in one connected cluster. The size in
the benchmark name is N.

| Benchmark | What it measures | 0.1.0 | P0 | Speed-up |
|---|---|---:|---:|---:|
| `assembly/solve_chain/5` | cold solve, each body mated to the previous | 287 µs | 66 µs | 4.3× |
| `assembly/solve_chain/20` | | 7.05 ms | 0.817 ms | 8.6× |
| `assembly/solve_chain/40` | | 47.6 ms | 4.22 ms | 11× |
| `assembly/solve_star/5` | cold solve, each body pinned to the ground | 316 µs | 75 µs | 4.2× |
| `assembly/solve_star/20` | | 8.58 ms | 0.884 ms | 9.7× |
| `assembly/solve_star/40` | | 61.1 ms | 4.28 ms | 14× |
| `assembly/diagnose_star/5` | rank-revealing diagnosis | 275 µs | 131 µs | 2.1× |
| `assembly/diagnose_star/20` | | 7.36 ms | 4.06 ms | 1.8× |
| `assembly/diagnose_star/40` | | 53.4 ms | 31.8 ms | 1.7× |

### After the move to the shared pipeline (P1.1)

The rigid-pose family now solves and diagnoses through the same pipeline as
the other two, and a body pinned in full by a `Fixed` constraint is
grounded (placed, not solved for). Same machine, same fixtures:

| Benchmark | P0 | P1.1 | Change |
|---|---:|---:|---:|
| `assembly/solve_chain/5` | 66 µs | 109 µs | 1.65× slower |
| `assembly/solve_chain/20` | 0.817 ms | 1.01 ms | 1.24× slower |
| `assembly/solve_chain/40` | 4.22 ms | 4.41 ms | level |
| `assembly/solve_star/5` | 75 µs | 88.5 µs | 1.18× slower |
| `assembly/solve_star/20` | 0.884 ms | 0.348 ms | 2.5× faster |
| `assembly/solve_star/40` | 4.28 ms | 0.695 ms | 6.2× faster |
| `assembly/diagnose_star/5` | 131 µs | 307 µs | 2.3× slower |
| `assembly/diagnose_star/20` | 4.06 ms | 6.50 ms | 1.6× slower |
| `assembly/diagnose_star/40` | 31.8 ms | 44.9 ms | 1.4× slower |

- The star is faster because its parts hang off the ground and now solve as
  independent six-column systems; it scales linearly in the number of parts.
- The chain is one coupled cluster either way, so only the fixed cost of
  the generic pipeline shows: noticeable at 5 bodies, gone by 40.
- The diagnosis pays that fixed cost too, and now also extracts dependency
  groups, which the assembly diagnosis did not report before.
- The 2D and 3D-geometric benchmarks moved by less than 15%, mostly down:
  the Newton loop no longer re-evaluates the residuals of a step it has
  just accepted.

Re-measured at the end of P1 (typed constraints, two-row `Parallel` and
four-row `Concentric`, restore on failure, the editing state): every
assembly benchmark is within 10% of the P1.1 column, in both directions.
That is no finding either way: these fixtures use coincident and angle
mates only, so the row change does not reach them, and a successful solve
never takes the restore path. A single run of `assembly/solve_chain/5`
read 125 µs; three repeats of it alone read 104 µs each, which is the
scatter to expect from one pass at 10 samples.

### After the step control became a trust region (P2.3)

Same machine. The dogleg takes the same steps as the old line search
wherever the full Gauss–Newton step is accepted, so the cost per iteration
is what matters here, not the iteration count.

| Benchmark | before | dogleg |
|---|---:|---:|
| `sketch/solve_chain/54` | 0.364 ms | 0.320 ms |
| `sketch/solve_chain/198` | 7.54 ms | 6.24 ms |
| `sketch/solve_chain/792` | 357 ms | 335 ms |
| `sketch/solve_field/3600` | 8.66 ms | 7.66 ms |
| `space/solve_truss/240` | 14.3 ms | 12.3 ms |
| `space/solve_truss/780` | 396 ms | 352 ms |
| `assembly/solve_chain/5` | 104 µs | 99.7 µs |
| `assembly/solve_chain/20` | 0.917 ms | 0.910 ms |
| `assembly/solve_chain/40` | 4.33 ms | 4.19 ms |
| `assembly/solve_star/40` | 0.692 ms | 0.677 ms |

The "before" column for the sketch and space rows is the P0 column above,
measured earlier, so the few percent in the dogleg's favour there are
within what separates two sessions; the assembly rows were measured the
same day. The reading is "level", not "faster".

A first version of the dogleg multiplied by the dense Jacobian three times
per iteration and made the 5-body chain 30% slower (132 µs). The products
now run over the non-zero entries of each row.

### After equations were weighed by their precision (P2.4)

Same machine. Near the origin every weight is one and the solve takes the
steps it took before, so this only asks what the bookkeeping costs. Two
mid-size cold solves, the previous commit and this one run alternately,
four rounds each (criterion medians, 6 s of measurement per run):

| Benchmark | before, four runs | mean | after, four runs | mean |
|---|---|---:|---|---:|
| `assembly/solve_chain/40` | 4.35 / 4.33 / 4.21 / 4.42 ms | 4.33 ms | 4.49 / 4.51 / 4.31 / 4.46 ms | 4.44 ms |
| `sketch/solve_chain/198` | 6.72 / 6.46 / 6.43 / 6.75 ms | 6.59 ms | 6.74 / 6.56 / 6.91 / 6.93 ms | 6.79 ms |

About 3% slower, against a spread of about 5% between runs of one build:
visible in the means, not in any single pair of runs. The small fixtures
are level. The cause was not found. It is not the weighting: both builds
do the same number of eliminations, normal-equation builds and residual
evaluations on these fixtures (6 / 6 / 8 for the sketch, 2 / 2 / 5 for the
assembly), and the difference stays when the weighting is switched off.
What is left is how the rewritten loop compiles. (Forcing the elimination
kernel inline made the solve 3.5× slower, so its code generation is
sensitive to the surroundings.) Worth a look when the factorisation is
replaced in P2.1 / P2.2, not before.

### After a small 2D component was held to its own size (P2.4)

Same machine. The fixtures are 4 × 3 rectangles, so every component is
larger than one unit and keeps the absolute tolerance; what is new on
them is the box around each component's entities, taken once per solve.
The `sketch/` group, the previous commit and this one run alternately,
three rounds each (criterion medians, 6 s of measurement per run):

| Benchmark | before, mean | after, mean | change |
|---|---:|---:|---:|
| `sketch/solve_chain/54` | 0.363 ms | 0.383 ms | +5.4% |
| `sketch/solve_chain/198` | 7.07 ms | 7.16 ms | +1.3% |
| `sketch/solve_chain/792` | 379 ms | 395 ms | +4.1% |
| `sketch/resolve_chain/54` | 0.362 ms | 0.385 ms | +6.5% |
| `sketch/resolve_chain/198` | 7.55 ms | 7.23 ms | −4.2% |
| `sketch/resolve_chain/792` | 388 ms | 390 ms | +0.6% |
| `sketch/solve_field/360` | 1.36 ms | 1.33 ms | −2.2% |
| `sketch/solve_field/3600` | 8.97 ms | 9.19 ms | +2.4% |
| `sketch/diagnose_chain/198` (code unchanged) | 9.17 ms | 9.01 ms | −1.8% |

No difference is resolved. The machine drifted during this run: one
build's own rounds differ by 10–20% (`solve_chain/54` before:
0.343 / 0.357 / 0.390 ms; `solve_chain/792` before: 349 / 370 / 420 ms),
the first round being the fastest for both builds; within that round the
pairs are within 3% either way (`solve_chain/54` 0.343 → 0.348 ms,
`solve_chain/792` 349 → 337 ms). The diagnosis, whose code did not change,
scatters by the same couple of percent.

### After the stored LU factorisation was made as fast as the elimination (P2.2)

Same machine. `LuFactors` (used by the diagnosis for its Gram system) did
the arithmetic of `gaussian_elimination` in two to three times the time.
Measured on dense symmetric positive definite matrices, one factorisation
and solve, release build:

| Unknowns | elimination | `LuFactors` before | `LuFactors` after |
|---:|---:|---:|---:|
| 54 | 0.029 ms | 0.056 ms | 0.026 ms |
| 198 | 0.60 ms | 2.00 ms | 0.74 ms |
| 240 | 1.15 ms | 3.54 ms | 1.23 ms |
| 780 | 41.1 ms | 112.8 ms | 44.0 ms |

(The “after” column is the slice-based loop with the factors still in a
vector per column; keeping them in the triangle removed those
allocations too.) The diagnosis benchmarks, the previous commit and this
one alternately, three rounds each (criterion medians):

| Benchmark | before | after | change |
|---|---:|---:|---:|
| `assembly/diagnose_star/5` | 0.314 ms | 0.319 ms | +1.5% |
| `assembly/diagnose_star/20` | 7.26 ms | 6.93 ms | −4.5% |
| `assembly/diagnose_star/40` | 52.2 ms | 50.3 ms | −3.7% |
| `space/diagnose_truss/60` | 0.489 ms | 0.461 ms | −5.8% |
| `space/diagnose_truss/240` | 42.2 ms | 39.4 ms | −6.5% |
| `sketch/diagnose_chain/54`, `/198` | | | +0.7%, −1.2% |
| `sketch/diagnose_chain/792`, `space/diagnose_truss/780` | | | +3.6%, +5.9% |

The last two spread by 10–15% between rounds of one build (792: 906 to
1020 ms before) and are not resolved; the first five are steady to a few
percent. (These were measured with the next commit's convergence polish
in place as well; the diagnosis does not run the polish.)

### After the convergence polish began reusing its factorisation (P2.2)

Same machine, all benchmarks, the previous commit (`dfe53b0`) and this
change alternately, three rounds each (criterion medians, 6 s of
measurement per run; rounds of one build agreed to a few percent):

| Benchmark | before | after | change |
|---|---:|---:|---:|
| `sketch/solve_chain/54` | 0.415 ms | 0.256 ms | −38% |
| `sketch/solve_chain/198` | 8.25 ms | 4.40 ms | −47% |
| `sketch/solve_chain/792` | 456 ms | 231 ms | −49% |
| `sketch/resolve_chain/54` | 0.419 ms | 0.258 ms | −38% |
| `sketch/resolve_chain/198` | 8.24 ms | 4.36 ms | −47% |
| `sketch/resolve_chain/792` | 444 ms | 226 ms | −49% |
| `sketch/solve_field/360` | 1.49 ms | 1.03 ms | −31% |
| `sketch/solve_field/3600` | 9.98 ms | 9.26 ms | −7% |
| `space/solve_truss/60` | 0.391 ms | 0.423 ms | +8% |
| `space/solve_truss/240` | 16.8 ms | 10.7 ms | −36% |
| `space/solve_truss/780` | 472 ms | 407 ms | −14% |
| `assembly/solve_chain/5` | 0.133 ms | 0.151 ms | +14% |
| `assembly/solve_chain/20` | 1.14 ms | 1.17 ms | +2.6% |
| `assembly/solve_chain/40` | 5.33 ms | 5.31 ms | −0.4% |
| `assembly/solve_star/5`, `/20`, `/40` | | | +1.6%, +2.1%, +2.4% |

Counted on the 198-parameter chain (temporary counters, not committed): a
cold solve from 0.05 off and a warm one from 1e-4 off both take two
iterations, and both then ran four polish steps, each with a new Jacobian
and a factorisation — four of the six factorisations of the solve. That
is why the warm re-solve was no cheaper; it still is not, since Newton
needs two iterations from either start, but both now cost three
factorisations' worth instead of six.

The two small regressions are not the polish. Counted the same way, the
assembly corpus (12 592 solves) and both small fixtures take exactly as
many iterations and polish steps before as after, and on the assembly
every polish is a single step. A variant whose first polish step was the
old code path verbatim (elimination, no stored factors) was as slow on
those two fixtures (+17%, +9%) and gained less elsewhere (−28% to −35% on
the 2D chains, as a second step then needs its own factorisation), so it
was not kept. What remains is how the larger function compiles, as with
the 3% noted for P2.4 above.

The first version of this change factored with `LuFactors` as it then
was, two to three times slower than the elimination: it made the assembly
chains 65–101% slower and led to the previous section.

### After ill-conditioned steps were taken by QR (P2.2)

Same machine, all benchmarks, the previous commit (`40d6f9b`) and this
change alternately, three rounds each. The fixtures are well conditioned:
the QR path is never taken on them, and every step is solved as before
(through `LuFactors`, bitwise what the elimination gives). 2D and 3D
solves are within ±5% (`sketch/solve_chain/54` −5.6%, `space/solve_truss`
60 / 240 / 780 +5.6% / −4.1% / +8.2%, the last from one round of 518 ms
against 380–393 ms in the others); the assembly solves are within ±2%
except the 5-body chain (+9%, 0.125 → 0.136 ms, one of three rounds level).
The diagnosis does not take this path; its scatter (`assembly/diagnose_star/20`
+25% from one round of 11.5 ms, `sketch/diagnose_chain/792` −16%) is the
machine.

QR on every step was measured first, against the same commit:

| Benchmark | normal equations | QR everywhere |
|---|---:|---:|
| `sketch/solve_field/360` (18 parameters per component) | 0.99 ms | +13% |
| `sketch/solve_chain/54` | 0.24 ms | +56% |
| `sketch/solve_chain/198` | 4.17 ms | +137% |
| `sketch/solve_chain/792` | 216 ms | +157% |
| `space/solve_truss/240` | 10.6 ms | +255% |
| `space/solve_truss/780` | 391 ms | +320% |
| `assembly/solve_chain/40` | 5.01 ms | +133% |

A Householder QR of `J` is four to six times the arithmetic of forming
`JᵀJ` over the non-zero entries and eliminating it, so it is kept for the
steps that need it. How often that is, counted with temporary counters:
75 of 18 317 steps in the 2D corpus, 16 of 37 781 in the assembly corpus,
17 of 41 in the 3D corpus.

### After the Jacobian was stored sparse (P2.1)

Same machine, all benchmarks, the previous commit (`02b6953`) and this
change alternately, three rounds each. The third round of this change ran
10–50% slower on every benchmark, the diagnosis ones included, which this
change does not touch; it is the machine, so the table gives the median of
the three rounds on each side.

| Benchmark | before | after | change |
|---|---:|---:|---:|
| `sketch/solve_chain/54` | 0.228 ms | 0.186 ms | −18% |
| `sketch/solve_chain/198` | 3.79 ms | 3.59 ms | −5.3% |
| `sketch/solve_chain/792` | 206 ms | 200 ms | −2.6% |
| `sketch/resolve_chain/54` | 0.225 ms | 0.186 ms | −17% |
| `sketch/resolve_chain/198` | 3.84 ms | 3.58 ms | −6.7% |
| `sketch/resolve_chain/792` | 207 ms | 197 ms | −4.6% |
| `sketch/solve_field/360` | 0.930 ms | 0.754 ms | −19% |
| `sketch/solve_field/3600` | 8.22 ms | 7.78 ms | −5.4% |
| `space/solve_truss/60` | 0.343 ms | 0.245 ms | −29% |
| `space/solve_truss/240` | 9.57 ms | 8.42 ms | −12% |
| `space/solve_truss/780` | 369 ms | 347 ms | −6.1% |
| `assembly/solve_chain/5`, `/20`, `/40` | | | −3.3%, −3.6%, −5.2% |
| `assembly/solve_star/5`, `/20`, `/40` | | | −4.6%, −3.9%, −6.3% |

The diagnosis benchmarks are +0.1% to +6.3% on the medians, with the
rounds of the two sides overlapping (`sketch/diagnose_chain/198` 8.78–9.18
ms before, 8.84–9.42 after, leaving out the slow round); their code is not
part of this change.

Where the time went before, from temporary timers in the Newton loop
(release build, per solve, not committed): on the 54-parameter chain the
dense Jacobian — zeroing it, filling it through the hash map, gathering its
non-zero entries three times, weighting it, and scanning every column of
every row to form the normal equations — was a third of the solve, and the
factorisation 28%; on the 240-parameter truss the factorisation was 72%
and the dense passes 21%. The factorisation is what is left: it is now 39%
of the 54-parameter chain and 83% of the 240-parameter truss, and it is
why the largest fixtures gain least.

### After the normal matrix and the layout stopped allocating per row (P2.1)

Same machine, all benchmarks, the previous commit (`9ecc6b4`) and this
change alternately, three rounds each. This run was noisier than the
previous one (rounds of one build differ by 15–20%), so the table gives
medians and the range of the rounds.

| Benchmark | before | after | change | rounds overlap |
|---|---:|---:|---:|---|
| `sketch/solve_chain/54` | 0.213 ms | 0.155 ms | −27% | no |
| `sketch/solve_chain/198` | 3.98 ms | 3.69 ms | −7.3% | yes |
| `sketch/solve_chain/792` | 224 ms | 229 ms | +2.6% | yes (195–242 / 194–252) |
| `sketch/resolve_chain/54` | 0.211 ms | 0.161 ms | −24% | no |
| `sketch/solve_field/360` | 0.833 ms | 0.549 ms | −34% | no |
| `sketch/solve_field/3600` | 7.87 ms | 5.52 ms | −30% | no |
| `space/solve_truss/60` | 0.258 ms | 0.172 ms | −34% | no |
| `space/solve_truss/240` | 9.58 ms | 8.76 ms | −8.6% | yes |
| `space/solve_truss/780` | 410 ms | 428 ms | +4.3% | yes (409–450 / 356–448) |
| `assembly/solve_chain/5`, `/20`, `/40` | | | −23%, −24%, +4.9% | no, no, yes |
| `assembly/solve_star/5`, `/20`, `/40` | | | −26%, −31%, −25% | no |
| diagnosis benchmarks | | | −2% to −25% (one +9.7%) | mostly yes |

Allocations per solve, counted with a temporary counting allocator: 644 →
359 on the 54-parameter chain, 3 068 → 733 on the 240-parameter truss,
613 → 311 on the 5-body assembly chain.

The large fixtures were settled by timing the factorisation alone, the
old row-of-vectors layout against the flat one on the same banded matrix,
alternately, five runs each in one binary (median, building the matrix
included): 1.91 → 0.83 µs at n = 18, 15.8 → 10.8 µs at 54, 637 → 476 µs
at 200, 37.8 → 36.5 ms at 790. The flat layout is no slower anywhere, so
the +3–5% medians above are the machine. The diagnosis gains because its
Gram system is factored by the same code.

A first quick run with only the flat matrix (the layout as before) gained
3–15% on the small fixtures; the rest of the gain at that size is the
layout, which for the rigid family had grown a fresh vector for every row
of every constraint, several reallocations each.

## Corpus gates

Wall time of the metric gates in `tests/`, same machine.

| Gate | 0.1.0 | P0 | P1.1 |
|---|---:|---:|---:|
| `corpus_metrics_within_thresholds` (2D, 320 cases), unoptimised | 6.2 s | 2.4 s | 2.4 s |
| `assembly_metrics_within_thresholds` (333 cases), unoptimised | 411 s | 69 s | 38 s |
| `assembly_metrics_within_thresholds`, optimised | 21.0 s | 1.7 s | 1.5 s |

## Reading the numbers

- **Where the P0 speed-up came from.** Almost all of it is one change: the
  normal equations `JᵀJ` are accumulated over the non-zero Jacobian entries
  instead of over every pair of columns. The results are bitwise the same.
- **The single-component solve is still cubic.** A roughly four-fold
  increase in parameters costs about 20× (54 → 198) and about 47×
  (198 → 792): the dense Gaussian elimination now dominates. One connected
  component of roughly 800 parameters takes a third of a second; that is the
  practical ceiling until the factorisation is sparse.
- **Independent components scale about linearly.** Ten times the rectangles
  (360 → 3600 parameters) costs about 7×.
- **A warm re-solve costs as much as a cold solve.** On the rectangle chain
  Newton needs two iterations from 1e-4 off as from 0.05 off, so there is
  nothing to save. What made both expensive was the convergence polish,
  four factorisations after the tolerance was met; it now factors once,
  and both take about half the time (P2.2 above).
- **Diagnosis is now the slow half.** Diagnosing the 792-parameter component
  (767 ms) takes twice as long as solving it, and with 255 redundant rows
  the 780-parameter truss takes 2.4 s. The dense column-pivoted QR was not
  touched in P0.
- **A hard case that is deliberately not a fixture.** Without its braces the
  truss is a serial chain where every point hangs off the three before it.
  From a 0.05 perturbation the solver needs 5 iterations at 20 points, 9 at
  80, 17 at 120, and does not converge within 100 at 160 or at 260 points
  (it does at 200, in 35). It is a conditioning problem for a future
  trust-region step to win, and a good corpus case for it.

## Not covered yet

- No automatic regression gate: nothing fails when these numbers get worse.
