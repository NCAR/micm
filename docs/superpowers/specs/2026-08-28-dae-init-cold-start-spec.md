# Spec — Damped Newton for DAE Constraint Initialization (PR #1083 follow-on)

- **Date:** 2026-08-28
- **Target branch:** `fix/dae-constraint-tolerance-measure` (NCAR/micm PR #1083, OPEN, approved, CI 29/29 green)
- **Base:** `main` @ `92f00b12` (2026-08-24) — already contains the cumulative-alpha fix (`c0290729`, PR #1047, merged 2026-08-07)
- **Related:** NCAR/musica#956 Case 2;
  `dae-rosenbrock-benchmark:docs/superpowers/notes/2026-07-26-musica956-knife-edge-repro.md`;
  ODE repo `paper_dae/` §4, `update_dae/update_dae.tex`

> **Outcome note (2026-08-29):** Q3 was answered from the issue screenshots.
> The reporter's Case 2 returns `StepSizeTooSmall` at t=0, not
> `ConstraintInitializationFailed`, and Case 4 fails during integration. The
> reduced reproducer below remains valid evidence for the initialization
> algorithm, but this work is not evidence that either reported case is fixed.

## 1. Problem

PR #1083 replaces the DAE constraint-initialization acceptance rule with an
affine-covariant weighted-correction test:

```text
scale_a = atol_a + rtol * max(|z_a|, |z_a + delta_a|)
q       = max_a |delta_a| / scale_a  <=  constraint_init_tolerance_   (default 0.1)
```

That change is correct and it resolves the reduced, musica#956-inspired knife
edge described below. It does **not** establish a fix for the reporter's actual
Case 2, and it does not globalize the Newton iteration. `InitializeConstraints`
still applies a full, undamped step every update (`rosenbrock.inl` step 9,
`apply_update(Y, delta)`).

For constraints nonlinear in the solved algebraic variable — the quadratic
dissociation `K2*[AQ] = [XM]^2` of real aqueous chemistry is the motivating case —
an undamped Newton step from a guess far below the root overshoots, and the
iteration then halves its way back down, exhausting
`constraint_init_max_iterations_` (default 10) before reaching the manifold.

### 1.1 Measured evidence

Reproducer: `benchmark/musica956_knife_edge.cpp`, ported to main's templated
constraint API. Henry's law (Van't Hoff, `K_ref = 0.5`, `delta_H = -60000 J/mol`)
feeding a quadratic aqueous dissociation, both as built-in
`EquilibriumConstraint`s, `P = 85000 Pa`, `T` swept 278–292 K.

Two starting guesses for the solved ion `XM`, whose root sits near 600–800:

| tree | warm start (`XM = 900`) | cold start (`XM = 1`) |
|---|---|---|
| `origin/main` @ `7c62ea76` | fails at 284, 287, 290 K | **fails at all 15 T** |
| PR #1083 @ `f72c25df` | converges at all 15 T | **fails at all 15 T** |
| `init-weighted-correction` @ `03dc6693` | converges at all 15 T | **converges at all 15 T** |

Both `main` and PR #1083 carry the cumulative-alpha fix, so the warm-start
difference is attributable to the initialization rule alone. The cold-start
column isolates the missing globalization: `init-weighted-correction` differs
from PR #1083 in exactly one respect that bears on it — it carries a damped
backtracking line search.

On failure PR #1083 rolls back transactionally (the caller's `AQ` is returned
unchanged at `1.0`), which is a genuine improvement over main's partially
modified state, but the projection still fails.

### 1.2 Root cause

`fix/dae-constraint-tolerance-measure` does not contain the line search. These
symbols exist on `init-weighted-correction` and `dae-rosenbrock-benchmark` and
are absent from the PR branch:

- `constraint_init_max_backtracks_`
- `constraint_init_backtrack_factor_`
- `constraint_init_sufficient_decrease_`

`DAE.md` documents all five initialization parameters as if they shipped
together. `update_dae/update_dae.tex` (ODE repo) correctly states that PR #1083
"does not add Newton damping or backtracking." The measurement above confirms
`update_dae.tex`; `DAE.md` describes the `dae-rosenbrock-benchmark` state, not
the PR.

## 2. Goal

Add affine-covariant damped Newton to PR #1083's `InitializeConstraints`, so
that consistent-initial-condition projection converges from starting guesses
outside the full-step convergence basin, and produce the numerical evidence
(CSV + figures) that demonstrates it.

## 3. Requirements

### R1 — Globalized Newton update

Replace the unconditional `apply_update(Y, delta)` with a backtracking line
search over step length `lambda`, starting at `lambda = 1` and multiplying by
`constraint_init_backtrack_factor_` on rejection.

### R2 — Affine-covariant merit function

The merit function MUST be the weighted norm of the **simplified** Newton
correction at the candidate — the correction obtained by solving with the
**already-factored** Jacobian from the current iterate, not a refactorization.
This is Deuflhard's natural monotonicity test (NLEQ). It is required, not
merely preferred: a merit function built from the raw residual `|G|` would
reintroduce exactly the row-scale dependence that PR #1083 exists to remove.

Accept a candidate when either holds:

```text
||delta_bar(candidate)||_w  <=  constraint_init_tolerance_
||delta_bar(candidate)||_w  <   (1 - sigma * lambda) * ||delta(Y)||_w
```

where `sigma = constraint_init_sufficient_decrease_` and `||.||_w` is the
existing weighted-correction norm.

### R3 — Parameters

Add to `RosenbrockSolverParameters`, with defaults matching
`dae-rosenbrock-benchmark`:

| Parameter | Type | Default |
|---|---|---|
| `constraint_init_max_backtracks_` | `Index` | `24` |
| `constraint_init_backtrack_factor_` | `Real` | `0.5` |
| `constraint_init_sufficient_decrease_` | `Real` | `1e-4` |

Setting `constraint_init_max_backtracks_ = 0` MUST reproduce PR #1083's current
undamped behavior exactly, so the change is opt-out and bisectable.

### R4 — Device-portable implementation

PR #1083's `InitializeConstraints` is device-abstracted: it uses
`DenseMatrixPolicy::MaxType<Real>` and `LOrType` reduction scalars with explicit
`CopyToDevice()` / `CopyToHost()`, and expresses element-wise work through
dispatched lambdas (`check_algebraic_values`, `check_weighted_correction`,
`apply_update`). The reference implementation on `dae-rosenbrock-benchmark` is
host-only and indexes matrices directly.

The port MUST be written in PR #1083's idiom. Direct transcription of the
reference loops will not compile for Kokkos/CUDA. The branch's tests currently
pass in double (65/65), Kokkos (90/90), and a float build (40 constraint/DAE
tests); all three MUST still pass.

### R5 — Buffers

Use existing dense temporaries; add no new dense-matrix allocation:

- `initial_forcing_` — residual / Newton correction `delta` (already used)
- `K_[0]` — rollback snapshot (already used)
- `Ynew_` — candidate iterate `Y + lambda * delta` (free at init time)
- `Yerror_` — candidate simplified correction `delta_bar` (free at init time)

Reuse of `Ynew_` and `Yerror_` is safe because `InitializeConstraints` runs
before the integration loop's first stage.

One `DenseMatrixPolicy::ScalarType<Real>` handle for the changing step length is
permitted. It is required for Kokkos because `MICM_LAMBDA` captures an ordinary
`Real` by value when the reusable `Function` is constructed. This is an
8-byte host scalar (and the corresponding scalar device view), not a per-cell
workspace allocation.

### R6 — Transactional failure preserved

Every failure exit MUST continue to route through `restore_and_return`, and the
line search MUST NOT write to `Y` until a candidate is accepted. A non-finite
candidate or a non-finite candidate correction backtracks rather than failing.

### R7 — Cost on the common path

For an already-consistent or near-consistent initial state the full step is
accepted on the first trial, so the added cost is bounded by one extra
constraint forcing evaluation plus one extra triangular solve per Newton update.
The exactly-zero-residual short circuit MUST still return before any
factorization. This cost MUST be measured, not assumed (see D2).

## 4. Non-goals

Out of scope for this change; each is a separate PR:

- The algebraic-variable LTE fix (removing `SetAlgebraicErrors`,
  `constraint_set.hpp:278` / `rosenbrock.inl:223`). Still present on main and
  on PR #1083.
- Schur-reduced stage solves.
- Post-clamp reprojection (`ReprojectAfterClamp`).
- The `constraint_init_min_pivot_ratio_` diagnostic. Optional; propose
  separately if reviewers want conditioning reporting.

## 5. Deliverables

- **D1 — Cold-start basin map.** A benchmark sweeping initial algebraic guess
  against temperature, emitting CSV, plus a three-panel figure (main /
  PR #1083 / PR #1083 + line search).
- **D2 — Initialization cost.** Work counters (function calls, factorizations,
  solves) for warm/consistent starts, PR #1083 vs PR #1083 + line search.
- **D3 — Row-scale invariance retained.** Confirmation that adding the line
  search does not reintroduce false convergences or false failures across row
  scales.
- **D4 — Convergence trace.** For one cold-start case, `||delta||_w` per Newton
  update, undamped vs damped.
- **D5 — Regression tests** in `test/unit/solver/test_constraint_initialization.cpp`.
- **D6 — Doc corrections** to `DAE.md` (which currently describes parameters the
  PR does not have).

## 6. Acceptance criteria

1. Cold start `XM = 1` converges at all 15 temperatures in the reproducer.
2. Warm start `XM = 900` still converges at all 15 temperatures, with the
   post-init `AQ` column bit-identical to PR #1083's.
3. `constraint_init_max_backtracks_ = 0` reproduces PR #1083 behavior exactly.
4. Full suite green in double, Kokkos, and float builds.
5. No change to any constraint-free (pure ODE) solve: zero extra work when
   `has_constraints` is false.
6. D1–D4 produced as committed CSV and figures.

## 7. Open questions for the follow-on session

- **Q1 — Ship inside PR #1083, or as a stacked follow-up PR?** #1083 is already
  approved with green CI; adding ~80 lines to `InitializeConstraints` will
  require re-review. A stacked PR keeps #1083 mergeable now. Recommendation:
  stacked follow-up, unless maintainers prefer one coherent change.
- **Q2 — Is `constraint_init_max_iterations_ = 10` still right with damping?**
  Each Newton update can now consume up to 25 residual evaluations. Consider
  whether the budget should be expressed in total forcing evaluations instead.
- **Q3 — Does the reporter's Case 2 actually return
  `ConstraintInitializationFailed`?** The issue screenshots were not transcribed
  into the thread. `StepSizeTooSmall` or `NaNDetected` would point at a
  different fix. Worth asking on the issue in parallel with this work.
