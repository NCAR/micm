# Results — damped Newton for DAE constraint initialization

- **Date:** 2026-08-29
- **Branch:** `dae-init-damping`, off `fix/dae-constraint-tolerance-measure` (NCAR/micm PR #1083) @ `f72c25df`
- **Spec:** `docs/superpowers/specs/2026-08-28-dae-init-cold-start-spec.md`
- **Plan:** `docs/superpowers/plans/2026-08-29-dae-init-cold-start.md`
- **Baseline:** `baseline.md` in this directory

## Summary

PR #1083 replaced the constraint-initialization acceptance rule with an
affine-covariant weighted-correction test, which resolved the musica#956 Case-2
knife edge but left the Newton iteration undamped, so cold starts still failed.
This work adds the backtracking line search. Cold starts now converge, nothing
that already converged is perturbed, and the already-consistent path costs
nothing extra.

## Test results

| Configuration | Baseline (before) | After |
|---|---|---|
| double | 65/65 | 65/65 |
| Kokkos serial | 90/90 | 90/90 |
| float (6 constraint/DAE targets) | 6/6 | 6/6 |

`constraint_initialization` grew from 12 to 21 gtest cases. **No pre-existing
test was modified**: the diff to `test_constraint_initialization.cpp` is
additions only. In particular the exact-count assertions
`EXPECT_EQ(loose_stats.constraint_init_iterations_, 1)` and
`EXPECT_EQ(tight_stats.constraint_init_iterations_, 2)` pass unchanged, which is
what an early `return Converged` from inside the line search would have broken.

## D1 — Cold-start reach at the default update budget

> **Correction.** An earlier revision of this note called the measurement below a
> "basin map" and reported the improvement as a widened basin of attraction. That
> framing is wrong and is corrected here. The undamped iteration converges from
> *every* starting depth tested if it is given enough Newton updates: it descends
> by halving, needing 14 updates from `XM0/XM* = 1e-3`, 24 from `1e-6`, 34 from
> `1e-9` and 48 from `1e-13`. What it cannot do is get there within the shipped
> `constraint_init_max_iterations_ = 10`. So the quantity below is the reach of a
> fixed update budget, not a boundary of the method, and the improvement is one
> of convergence *rate*: the line search reaches the root in 4-5 updates where the
> undamped iteration needs 14-48. See "Where the line search is worse" below.

`dae_init_basin.csv`, `dae_init_basin.png`. The musica#956 Case-2 mechanism over
15 temperatures (278-292 K) x 9 initial guesses for the solved ion = 135
projections per variant.

| Variant | Converged | Warm `XM = 900` | Cold `XM = 1` |
|---|---:|---:|---:|
| `origin/main` | 84/135 | 12/15, failing at 284, 287, 290 K | 0/15 |
| PR #1083 (undamped) | 120/135 | 15/15 | 0/15 |
| PR #1083 + line search | **135/135** | 15/15 | **15/15** |

This reproduces the spec's §1.1 table independently, including `origin/main`'s
warm-start failures landing on exactly 284, 287 and 290 K.

Two properties beyond the counts:

- **Zero regressions.** No `(temperature, initial guess)` cell converges undamped
  and fails damped. The behavioral discontinuity the plan flagged as the main
  risk of enabling damping did not materialize anywhere in this sweep.
- **Zero perturbation.** The 120 cells that converge under both variants agree
  **bit for bit** in the post-initialization `AQ` column. This is spec acceptance
  criterion 2, demonstrated across 120 cells rather than one.

## D2 — Initialization cost

`dae_init_cost.csv`, `dae_init_cost.png`. Median of 200 projections, `T = 286 K`.

| Start | Variant | Newton updates | Jacobians | Factorizations | Solves | Median |
|---|---|---:|---:|---:|---:|---:|
| already consistent | undamped | 2 | 1 | 1 | 1 | 0.459 us |
| already consistent | + line search | 2 | 1 | 1 | 1 | 0.542 us |
| warm (`XM = 900`) | undamped | 4 | 4 | 4 | 4 | 1.125 us |
| warm (`XM = 900`) | + line search | 4 | 4 | 4 | 7 | 1.625 us |

This confirms spec R7 as a measurement rather than an assumption. The
already-consistent path does **identical** work: the exactly-zero-residual short
circuit returns before any factorization, so the line search is never entered.
A warm start pays one extra triangular solve per Newton update and **no extra
factorization**, because every trial reuses the factorization taken at the
current iterate. That is the point of the simplified-Newton merit function.

## D3 — Row-scale invariance retained

Pinned by `ConstraintInitialization.BacktrackingIsInvariantToConstraintRowScaling`
over the nonlinear projection system `row_scale * ([Z]^2 - [X]) = 0` at row
scales `1e-12`, `1` and `1e12`. The test asserts **exactly equal**
`constraint_init_iterations_` and `solves_` across all three, not merely equal
answers: the line search made the identical accept/reject decision at every
trial. Scaling a complete constraint row scales `G` and `dG/dy` together, so the
simplified Newton correction is unchanged, and the merit built from it inherits
that invariance.

The pre-existing `WeightedCorrectionIsInvariantToConstraintRowScaling` test does
not cover this: its constraint is linear in the algebraic variable, so the
simplified correction at the candidate is exactly zero and the search never
engages.

## D4 — Convergence trace

`dae_init_trace.csv`, `dae_init_trace.png`. One cold start, `T = 284 K`,
initial `XM = 1`. Weighted correction norm per Newton update:

| Update | Undamped | Damped |
|---:|---:|---:|
| 0 | 9.983e+05 | 9.983e+05 (accepts `lambda = 2^-9 = 1.95e-3`) |
| 1 | 4.983e+05 | 9.964e+05 |
| 2 | 4.967e+05 | 1.642e+04 |
| 3 | 4.934e+05 | 3.192e+02 |
| 4 | 4.868e+05 | 1.168e-01 |
| 5 | 4.738e+05 | 1.565e-08, converged |
| 10 | 5.829e+04, budget exhausted | - |

The undamped iteration is not diverging, it is creeping: it runs out of budget
still five orders of magnitude from the manifold. One heavily damped step puts
the iterate inside the convergence basin, after which every accepted step is a
full Newton step (`lambda = 1`) and the norm falls quadratically.

**Provenance:** `SolverStats` does not expose the per-update norm, and adding a
field for a one-off measurement was rejected as widening the reviewed surface.
This trace was produced from a scratch copy of the tree (`git archive` of this
branch) carrying a temporary `fprintf` in `InitializeConstraints`. The live tree
was never instrumented and carries no trace code.

The per-update *norm* above needs that instrumentation, but the per-update
*iterate* — the algebraic variable's value after each Newton update — does not.
`InitializeConstraints` reaches its `Converged` exit through a plain `return`
rather than through `restore_and_return`, so a projection whose acceptance
tolerance is loose enough to stop on pass `j` hands back exactly the iterate
after `j` applied updates and reports `constraint_init_iterations_ == j + 1`.
Bisecting on the tolerance therefore recovers the whole iterate path from a
stock build. That path, before and after the fix, is plotted in the companion
ODE repository (`src/dae_init_cold_start.cpp`,
`scripts/plot_dae_init_cold_start.py`), and the driver cross-checks its harvest
against a default-tolerance run so a divergence fails loudly rather than
producing a plausible figure.

## Where the line search is worse

Measured on the knife edge at `T = 286 K`, sweeping `constraint_init_max_iterations_`:

| `XM0/XM*` | undamped, budget 10 | undamped, budget 60 | damped, any budget |
|---:|---|---|---|
| 1e-3 | fail | converges, 14 updates | converges, 4 |
| 1e-6 | fail | converges, 24 | converges, 5 |
| 1e-9 | fail | converges, 34 | **fail** |
| 1e-13 | fail | converges, 48 | **fail** |

Below `XM0/XM* ~ 3e-8` the damped path fails on its **first** update
(`constraint_init_iterations_ == 1`): the line search exhausts all 24 backtracks
without finding sufficient decrease. That floor is `2^-24 ~ 6e-8`, set by
`constraint_init_max_backtracks_`, and it does not move with the iteration
budget. An undamped run given ~48 updates converges from there.

This is a real trade the change makes, and it should be stated in review rather
than discovered by a referee. It does not affect the shipped defaults, where the
damped path strictly dominates: at `max_iterations_ = 10` the undamped iteration
reaches nothing below `XM0/XM* ~ 8e-3`. Raising `constraint_init_max_backtracks_`
extends the damped floor if an application needs it.

## Known limitation: heterogeneous batches

Found while writing the multi-cell tests, and pinned by
`ConstraintInitialization.HeterogeneousBatchIsNotRescuedByTheLineSearch`.

The projection reduces one weighted correction norm across the whole batch and
damps it with one step length, so a step short enough for the worst cell is
applied to every cell. Measured on the square-root system:

| Batch | Undamped | Damped |
|---|---|---|
| 1 cell, cold | Failed | Converged (5 updates, 28 solves) |
| 3 cells, all cold | Failed | Converged (5 updates, 28 solves) |
| 3 cells, all warm | Converged | Converged |
| 2 cells, mixed cold + warm | **Failed** | **Failed** |

A batch whose cells share a cold start converges exactly as the single-cell case
does. A mixed batch fails **both** with and without the line search, which makes
this a pre-existing property of the batched formulation rather than a regression
introduced by damping. Per-cell step lengths would be the fix, and would be a
separate change to the reduction machinery.

This matters for MICM specifically, because batching grid cells is the normal
mode of operation. It should be raised with maintainers as follow-on work.

## Acceptance criteria (spec §6)

| # | Criterion | Status | Evidence |
|---|---|---|---|
| 1 | Cold start `XM = 1` converges at all 15 T | **Met** | `dae_init_basin.csv`, damped, 15/15 |
| 2 | Warm start still converges, post-init `AQ` bit-identical to #1083 | **Met** | 120 shared cells, zero differing bits |
| 3 | `constraint_init_max_backtracks_ = 0` reproduces #1083 exactly | **Met** | `DisablingBacktracksReproducesUndampedBehavior`; verbatim opt-out block |
| 4 | Suite green in double, Kokkos, float | **Met** | 65/65, 90/90, 6/6 |
| 5 | No change to constraint-free solves | **Met** | `PureODESystemUnaffected`; guarded by `if (has_constraints)` at `rosenbrock.inl:47` |
| 6 | D1-D4 committed as CSV and figures | **Met** | this directory |

## Spec open questions

- **Q1 (PR shape)** — deferred by decision until this evidence existed. See the
  recommendation in the final report.
- **Q2 (iteration budget)** — settled: unchanged at 10 Newton updates, with the
  cost measured in D2 rather than capped.
- **Q3 (does the reporter's Case 2 actually return `ConstraintInitializationFailed`?)**
  — **answered, and the answer is no.** The two screenshots in musica#956 were
  retrieved and read:

  | Reported case | Failure state | Time |
  |---|---|---|
  | Case 2 (`T = 286 K`, `P = 85000`, `LWC = 0.3e-3`) | `StepSizeTooSmall` | t = 0.0000 s |
  | Case 4 (`T = 285 K`, `P = 85000`, `LWC = 0.03e-3`) | `ConvergenceExceededMaxSteps` | t = 7170 s |

  Neither is `ConstraintInitializationFailed`. **This work therefore cannot be
  said to fix the reported cases**, and no such claim should be made in the PR
  description or the paper.

  What is established is narrower: the reproducer matches the *phenomenology* of
  Case 2 — a failure at t = 0 that flips with roughly 1 K of temperature at
  `P = 85000 Pa` — and PR #1083 plus this line search resolve that. It does not
  match the reported error state, and its failing temperatures (284, 287, 290 K)
  are not the reporter's (286 K fails, 285 K and 280 K pass).

  A plausible chain would connect them: on `main` the raw-residual acceptance
  rule can report `Converged` for a state that is off the manifold, after which
  the integrator cannot take a step and reports `StepSizeTooSmall` at t = 0. If
  that is the mechanism, the weighted-correction rule in #1083 is the fix.
  **This is untested.** Settling it needs the reporter's actual configuration
  (mechanism, species, LWC, initial conditions) run against `main`, #1083, and
  this branch, with the returned `SolverState` recorded in each.

  Case 4 fails mid-integration and is not an initialization failure at all.
  Nothing in this work addresses it.
