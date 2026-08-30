# DAE Rosenbrock Constraint Convergence

This document describes the consistent-initial-condition projection on
`dae-init-damping`, stacked on `fix/dae-constraint-tolerance-measure`. The work
makes projection independent of constraint units and row scaling, globalizes the
Newton iteration, and restores the caller's state after initialization failure.
It does not change integration error control or reproject after a post-solve
physical-bound clamp; those are separate changes.

## DAE form used by MICM

MICM treats constrained systems as a semi-explicit, index-1 DAE:

```text
x' = f(x, z)
0  = G(x, z)
```

- `x` contains differential variables.
- `z` contains algebraic variables.
- `G` contains the algebraic constraints.
- Differential rows have a mass-matrix diagonal of one.
- Algebraic rows have a mass-matrix diagonal of zero.

Before Rosenbrock time integration begins, MICM holds `x` fixed and projects
`z` onto `G(x,z)=0`. This projection is also performed at the start of every
public `Solve()` call because callers may alter the state between calls.

## Why the previous convergence test was inadequate

The previous initializer stopped when the largest raw constraint residual was
smaller than a fixed number:

```text
max |G_a(x, z)| < 1e-10
```

That test had no invariant numerical meaning:

- A residual has the units and scale chosen for its constraint equation.
- Multiplying a complete constraint row by a constant changed whether the same
  physical state was classified as converged.
- A small residual did not necessarily imply a small error in `z` when the
  algebraic Jacobian was poorly conditioned.
- An undamped Newton update could leave the local convergence basin.
- The residual was checked before each update, but not after the final permitted
  update. A solve that reached the exact solution on its last update could still
  report failure.
- Failed initialization left the caller with a partially modified state.
- Clamping a negative differential variable after integration could move the
  state away from the algebraic manifold.

In the manufactured projection sweep, the raw-residual rule produced 564 false
convergence classifications and 842 cases that reached the solution but reported
failure out of 5,400 cases.

## New convergence contract

The initializer now measures the Newton correction in the same weighted state
space used to define integration accuracy. For each algebraic variable `a`, it
forms

```text
scale_a = atol_a + rtol * max(|z_a|, |z_a + delta_a|)

q = max_a |delta_a| / scale_a
```

where `delta` solves the linearized algebraic system. Initialization converges
when

```text
q <= constraint_init_tolerance
```

The default `constraint_init_tolerance` is `0.1`. Thus, the remaining estimated
algebraic correction must be no more than one tenth of the state-error scale
defined by the configured absolute and relative tolerances.

This criterion has three useful properties:

- It is dimensionless.
- It is tied directly to the user's state tolerances.
- Multiplying a complete residual and Jacobian row by the same nonzero constant
  does not change the Newton correction or convergence decision.

An exactly zero algebraic residual is accepted immediately without a matrix
factorization. NaN and infinite values are reported through the corresponding
solver states.

### Parameter semantics

The constraint-initialization controls are members of
`RosenbrockSolverParameters`:

| Parameter | Default | Meaning |
|---|---:|---|
| `constraint_init_max_iterations_` | `10` | Maximum Newton updates |
| `constraint_init_tolerance_` | `0.1` | Maximum weighted correction `q` |
| `constraint_init_max_backtracks_` | `24` | Maximum line-search reductions per update |
| `constraint_init_backtrack_factor_` | `0.5` | Multiplier applied to a rejected step length |
| `constraint_init_sufficient_decrease_` | `1e-4` | Required fractional decrease in the line-search merit function |

`constraint_init_tolerance_` no longer represents a raw residual threshold.
Applications that explicitly set this field should review the new dimensionless
meaning rather than carrying forward values such as `1e-10` automatically.

## Globalized Newton iteration

For a current state `y=(x,z)`, one initialization iteration performs the
following operations:

1. Evaluate `G(y)` and return immediately if the residual is exactly zero.
2. Form and factor the algebraic Newton matrix while retaining identity rows for
   differential variables.
3. Solve for the full Newton correction `delta`; only algebraic variables are
   eligible for update.
4. Form a candidate `y(lambda)=y+lambda*delta`, initially with `lambda=1`.
5. Evaluate the candidate residual and solve for a simplified Newton correction
   using the already-factored current Jacobian.
6. Accept the candidate if its weighted correction is converged or satisfies
   the sufficient-decrease test.
7. Otherwise reduce `lambda` by the backtrack factor and retry.

Setting `constraint_init_max_backtracks_` to `0` disables steps 4-7 entirely and
applies the full Newton step unconditionally, reproducing the undamped update bit
for bit. This is an exact opt-out, which makes the line search bisectable, but
`0 -> 1` is a behavior discontinuity rather than a tightening: with the search
enabled, an update whose trials all fail the sufficient-decrease test is a
failure, where an undamped update would have applied the full step regardless.

The merit function is the norm of the simplified Newton correction, rather than
the raw residual norm. This makes the line search affine-covariant and preserves
complete-row-scaling invariance. Reusing the existing factorization keeps each
backtrack less expensive than a new Newton iteration.

The line search does not report convergence from inside itself. An accepted
candidate is taken and control falls through to the next Newton update, which
re-evaluates the residual at the new state and decides through the same
zero-residual and weighted-correction tests an undamped update uses. Reporting
convergence on acceptance would consume one fewer Newton update and so change
the meaning of `constraint_init_iterations_`, which callers and tests already
depend on.

## Transactional failure behavior

The complete state-variable matrix is saved before projection. If initialization
fails because of non-finite values, an unacceptable line search, or exhaustion
of the Newton iteration limit, all differential and algebraic variables are
restored before returning. Conditions and custom parameters are not modified by
the projection.

This provides a clear caller contract:

- `Converged` means the algebraic state passed the weighted-correction test.
- `ConstraintInitializationFailed`, `NaNDetected`, or `InfDetected` does not
  leave a partial Newton update in the caller's state.

An exactly singular algebraic Jacobian can still return `InfDetected`. For
example, the constraint `z^2-1=0` cannot define a Newton direction at `z=0`.
Finite initial values as close as `z=+/-1e-6` converge with backtracking.

## Not in this branch

Two related pieces of the DAE work are described in the design notes but are
**not** present here, and this document is scoped to what this branch actually
contains:

- **Post-clamp reprojection.** `Solver::PostSolveClamp` exists and correctly
  clamps only differential variables in a DAE state
  (`include/micm/solver/solver.hpp:275`), but nothing reprojects afterwards.
  Changing a differential variable generally changes the algebraic solution, so
  a clamp can still leave the algebraic variables off the manifold. There is no
  `ReprojectAfterClamp` in this source tree.
- **The algebraic-variable local-truncation-error fix.** The raw algebraic
  step-change error override is still applied: see
  `constraints_.SetAlgebraicErrors(Yerror, Y, Ynew)` at `rosenbrock.inl:223`
  and `constraint_set.hpp:278`.

Each is a separate change with its own review surface.

## What did not change

The work deliberately leaves the Rosenbrock integration controller unchanged:

- Accepted time steps are still governed by the embedded local-truncation-error
  estimate.
- There is no separate nonlinear constraint-residual rejection test after every
  accepted time step.
- The error norm remains a global weighted RMS over the state and all cells.
- Each public `Solve()` call still starts from its configured `h_start` policy.

The new `constraint_init_tolerance_` controls projection accuracy. It does not
directly control time-step acceptance.

## Numerical evidence

All numbers below were measured on this branch. The benchmark that produces them
is `benchmark/dae_init_cold_start.cpp`; its CSV output and figures are committed
under `docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/`.

### Cold-start reach at the default update budget

This is deliberately not described as a widened basin of attraction. The undamped
iteration converges from every depth tested here given enough Newton updates (48
from `XM0/XM* = 1e-13`); it descends by halving. What the line search changes is
the number of updates needed - 4 or 5 rather than 14 to 48 - and therefore what
is reachable inside `constraint_init_max_iterations_ = 10`. Below
`XM0/XM* ~ 2^-24` the damped path fails on its first update, where an undamped
run with a raised iteration budget would converge; that floor is set by
`constraint_init_max_backtracks_`.

The reduced Henry/dissociation reproducer inspired by musica#956 Case 2
(`P = 85000 Pa`), swept over 15 temperatures from 278 K to 292 K and 9 initial
guesses for the solved ion, is 135 projections per variant:

| Variant | Converged | Warm start `XM = 900` | Cold start `XM = 1` |
|---|---:|---:|---:|
| `origin/main` | 84/135 | 12/15, failing at 284, 287, 290 K | 0/15 |
| Weighted correction, undamped | 120/135 | 15/15 | 0/15 |
| Weighted correction + line search | **135/135** | 15/15 | **15/15** |

Two properties of the middle-to-bottom transition matter as much as the counts:

- **No cell regresses.** There is no `(temperature, initial guess)` cell that
  converges undamped and fails damped.
- **Nothing already working is perturbed.** The 120 cells that converge under
  both agree bit for bit in the post-initialization algebraic column.

### Convergence trace

One cold start, `T = 284 K`, initial `XM = 1`, weighted correction norm per
Newton update:

| Update | Undamped | Damped |
|---:|---:|---:|
| 0 | 9.98e+05 | 9.98e+05 (accepts `lambda = 2^-9`) |
| 1 | 4.98e+05 | 9.96e+05 |
| 2 | 4.97e+05 | 1.64e+04 |
| 3 | 4.93e+05 | 3.19e+02 |
| 4 | 4.87e+05 | 1.17e-01 |
| 5 | 4.74e+05 | 1.56e-08, converged |
| 10 | 5.83e+04, budget exhausted | - |

The undamped iteration is not diverging; it is creeping, and it runs out of
budget still five orders of magnitude away. One heavily damped step puts the
iterate inside the basin, after which every step is a full Newton step and the
norm falls quadratically.

### Initialization cost

Median of 200 projections at `T = 286 K`:

| Start | Variant | Residual evaluations | Jacobians | Factorizations | Solves | Median |
|---|---|---:|---:|---:|---:|---:|
| already consistent | undamped | 1 | 0 | 0 | 0 | 0.208 us |
| already consistent | + line search | 1 | 0 | 0 | 0 | 0.208 us |
| warm (`XM = 900`) | undamped | 4 | 4 | 4 | 4 | 1.125 us |
| warm (`XM = 900`) | + line search | 7 | 4 | 4 | 7 | 1.667 us |

An exactly consistent projection does **identical algorithmic work** with the
search enabled: one residual evaluation and no Jacobian, factorization, or
solve. The measured median is also identical at this resolution. For the warm
start, three applied updates add three candidate-residual evaluations and three
triangular solves, with no extra factorization, because every trial reuses the
factorization taken at the current iterate. Its measured median rises from
1.125 us to 1.667 us (about 48%).

### Known limitation: heterogeneous batches

The projection reduces one weighted correction norm across the whole batch and
damps it with one step length, so a step short enough for the worst cell is
applied to every cell. A batch mixing a cold cell with a warm one therefore fails
**both** with and without the line search: damping neither fixes nor worsens it.
A batch whose cells share a cold start converges exactly as the single-cell case
does. Per-cell step lengths would be the fix and would be a separate change to
the reduction machinery. This is pinned by
`ConstraintInitialization.HeterogeneousBatchIsNotRescuedByTheLineSearch`.

## Remaining numerical policy questions

### Global RMS dilution across cells

With one active cell and an increasing number of stationary cells, the existing
global RMS norm produced:

| Cells | Accepted steps | Active-cell relative error |
|---:|---:|---:|
| 1 | 44 | `3.73e-7` |
| 10 | 32 | `1.20e-6` |
| 100 | 24 | `3.81e-6` |
| 1000 | 18 | `1.22e-5` |

This is the expected square-root dilution of a global RMS norm, but it makes
per-cell accuracy depend on batch composition. Changing to a maximum cellwise
WRMS norm would be a separate integration-controller policy decision.

### Nearly singular algebraic blocks

For a manufactured coupled system with condition parameter `epsilon=1e-14`, the
computed residual was zero in working precision while the algebraic forward
error was `2.22e-2`. Neither a residual test nor a correction test can certify
forward accuracy once the algebraic block is numerically singular. A future
condition estimate or pivot-quality diagnostic would address this case.

### Repeated public solve calls

Splitting one interval into `1`, `10`, `100`, or `1000` public `Solve()` calls
increased accepted steps from `9` to `90`, `900`, and `9000`, respectively,
while the result remained accurate to roundoff. Carrying a step-size estimate
between calls would be a separate performance/API change.

## Implementation and validation files

- `include/micm/solver/rosenbrock.inl`: weighted-correction projection, the
  backtracking line search, and transactional rollback.
- `include/micm/solver/rosenbrock_solver_parameters.hpp`: projection controls
  and defaults.
- `include/micm/solver/rosenbrock_temporary_variables.hpp`: the projection's
  reused buffers and the line-search step-length scalar.
- `test/unit/solver/test_constraint_initialization.cpp`: projection, rollback,
  row scaling, and line-search tests.
- `test/integration/test_dae_constraint_overshoot.cpp`,
  `test/integration/test_dae_algebraic_error_insensitivity.cpp`:
  integration-level DAE controller tests.
- `benchmark/dae_init_cold_start.cpp`, `benchmark/plot_dae_init_cold_start.py`:
  the cold-start basin, cost and trace evidence.
- `docs/superpowers/specs/2026-08-28-dae-init-cold-start-spec.md` and
  `docs/superpowers/plans/2026-08-29-dae-init-cold-start.md`: the spec this work
  implements and the plan it followed.
- `FABLE_PLAN.md`: the updated Phase 0b roadmap status and upstream boundary.

## Reproducing the evidence

```bash
# Tests, in the three configurations this work is gated on
cmake -S . -B bld-double -DMICM_ENABLE_TESTS=ON -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER
cmake --build bld-double -j10 && (cd bld-double && ctest)

# Basin map, cost counters and their figures
cmake -S . -B bld-bench -DMICM_ENABLE_DAE_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build bld-bench --target dae_init_cold_start -j10
(cd bld-bench && ./dae_init_cold_start)
python3 benchmark/plot_dae_init_cold_start.py bld-bench <output-dir>
```

`dae_init_cold_start` writes `dae_init_basin.csv` and `dae_init_cost.csv` to its
working directory. The `origin/main` panel of the basin figure is produced by
building the same program against `origin/main` with the line-search knob
removed, and the per-update trace in `dae_init_trace.csv` comes from a scratch
build carrying a temporary `fprintf` in `InitializeConstraints`; neither is
reproduced by the command above. Both are recorded in the results note beside
the committed CSVs.
