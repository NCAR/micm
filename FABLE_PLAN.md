# Fable Plan — DAE Solver Corrections and Improvements

- **Date:** 2026-07-26; Phase 0b status updated 2026-08-29
- **Branch base:** `dae-rosenbrock-benchmark` (post `5ef0ce01`, QDE-S stiff family)
- **Predecessor:** `dae-rosenbrock-benchmark:PLAN.md` (2026-07-20
  performance-and-rigor plan, all phases done)
- **Source:** 2026-07-26 code review of the DAE Rosenbrock path (`rosenbrock.inl`,
  `rosenbrock_solver_parameters.hpp`, `constraint_set.hpp`, `schur_stage_solver.hpp`,
  `solver_builder.inl`, `solver.hpp`, benchmark drivers) against the v5 paper
  (ODE repo `paper_dae/`).
- **Status:** Phase 0 is on `main` through PR #1047. Phase 0b's weighted
  initialization work is on `fix/dae-constraint-tolerance-measure`; its damped
  Newton follow-on and numerical evidence are complete on `dae-init-damping`.
  The line search remains a stacked follow-up; upstream review is pending.
  Phases 1–9 are not started.

## Branch status

| Phase | Branch | Status |
|---|---|---|
| 0 — rejection-alpha fix | `main` | merged through PR #1047; regression coverage retained |
| 0b — constraint initialization | `fix/dae-constraint-tolerance-measure`; `dae-init-damping` | weighted-correction projection, rollback, damped Newton, tests, and D1-D4 evidence complete locally; stacked follow-up chosen, review pending |
| 1–9 | — | not started |

## Purpose

Three threads, in strict priority order:

1. **Correctness** — one confirmed silent-wrong-answer bug in the rejection
   path (Phase 0). Everything else waits until it is fixed, because several
   later phases re-measure benchmarks that run through the affected builder.
2. **Performance** — the review found overheads that are artifacts of
   implementation, not of the DAE formulation: controller losses the
   literature already quantifies, a fill-reducing ordering computed on the
   wrong sparsity pattern, `std::pow` on unit stoichiometries in the batched
   hot loop, a full sparse factorization per initialization sweep, and
   alpha-independent work redone per step attempt.
3. **Capability** — solver-level features the paper currently works around in
   benchmark drivers: a supported QSSA constraint type, a runtime conditioning
   monitor to replace the manual twilight guard, and non-autonomous stages.

Phases are ordered so correctness lands first and measurement-affecting
changes land before the re-measurements that judge them. Dev branches fork
from `dae-rosenbrock-benchmark` and merge back independently. Phases 0, 1,
and 2 are upstream-relevant (they fix or improve ODE-path behavior shared
with NCAR/micm `main`) and should be shaped as cherry-pickable commits.

## Correctness thread

### Phase 0 — Fix the rejection-path alpha accumulation bug
- **Branch:** `fix-rejection-alpha` (cherry-pick to an upstream PR against
  NCAR/micm `main`; the bug is upstream, introduced by PR #691 commit
  `52ea1138` and still present on `main` as of 2026-07-26)
- **Evidence:** `rosenbrock.inl:167-177` (non-in-place LU path). The delta
  trick stores the *delta* where the cumulative shift belongs:
  `alpha -= last_alpha; last_alpha = alpha;`. Trace: attempt 1 applies and
  stores a1 (correct); attempt 2 applies a2−a1 (cumulative a2, correct) but
  stores a2−a1; attempt 3 applies a3−(a2−a1), factoring `(a3+a1)M − J`
  instead of `a3·M − J`. From the third attempt of a single step — exactly
  the two-successive-rejections case that also triggers
  `rejection_factor_decrease_` — the stage matrix disagrees with the `c/H`
  right-hand-side terms, which is not a valid Rosenbrock step at any H. The
  embedded error estimate is built from the same corrupted stages and the
  over-damped matrix biases it toward acceptance: silent wrong answers, no
  NaN. Affects `CpuSolverBuilder` (non-in-place LU — the builder every DAE
  benchmark uses, e.g. `ts1_dae.cpp:631`), ODE and DAE alike. In-place and
  CUDA builders are unaffected.
- **Tasks:**
  - One-line fix: `last_alpha += alpha;` (or compute `alpha_full` first and
    store it). Keep the delta trick itself — it is valid once the
    accumulation is right.
  - Regression test: non-in-place builder, force ≥2 successive rejections of
    one step (absurdly large `h_start_` on Robertson at tight rtol), assert
    the result matches the in-place builder to roundoff. Add a unit-level
    assertion that after k attempts the shifted diagonal equals
    `1/(H_k·gamma) − J_ii` for a 2×2 system.
  - Open the upstream PR (fix + test) independent of any DAE machinery.
  - Re-run the benchmarks whose runs can reject: work–precision sweeps at
    loose rtol (`robertson_dae`, `tropospheric_dae`), the 5400-case
    projection sweep and 200-mechanism structural suite
    (`dae_convergence_experiments`), and `qde_exact` adaptive runs. TS1 had
    zero rejections and is exempt. Record any numeric drift in a note and
    propagate to the paper if visible at reported precision.
- **Acceptance:** forced-rejection test green on both builders with
  bit-identical accepted trajectories between them; full suite green;
  re-measured numbers recorded (expected: no visible drift, since published
  runs had rare/no rejections — but that is a measurement, not an
  assumption).
- **Status (2026-07-26):** fix and regression test committed on
  `fix-rejection-alpha` (`a826c8d0`) and pushed to NCAR/micm. The fix tracks
  the cumulative shift (`cumulative_alpha`) explicitly; the test
  (`RosenbrockSolver.RejectedStepAlphaAccumulation` in `test_rosenbrock.cpp`)
  runs Robertson with `h_start_` set to the full interval — six successive
  first-step rejections — and pins the standard builder to the in-place
  builder on attempt/accepted/rejected counts and final state (1e-10
  relative). Verified on the dev machine via a standalone header-only mirror
  of the test (no cmake/gtest locally): **pre-fix, upstream `main` diverges
  by up to 2.8e-3 relative with 13 extra attempts and 7 spurious mid-run
  rejections; post-fix the builders agree to 3e-15 with identical step
  sequences.** The acceptance criterion's "bit-identical" is refined to that
  observed 1e-15-scale agreement: the delta-shift path accumulates alpha in a
  different floating-point order than fresh assembly, so last-ulp differences
  are expected and the committed test asserts 1e-10. The full test suite
  passes locally on the fix branch: 64/64 via ctest (Release, AppleClang,
  2026-07-26), including the committed regression test under gtest and the
  analytical_rosenbrock_integration accuracy tests. Still pending: the
  upstream PR (branch is pushed; open via
  github.com/NCAR/micm/pull/new/fix-rejection-alpha), GitHub CI across
  platforms, and the benchmark re-runs. Real-world relevance: NCAR/musica#956's
  mid-run cloud-chemistry failure (t = 7170 s at reduced LWC) matches this
  bug's rejection-cascade signature, and MUSICA's `src/micm/cpu_solver.cpp`
  builds all CPU solvers, including `RosenbrockDAE4/DAE6`, with the affected
  `CpuSolverBuilder`; a diagnosis comment is drafted but not yet posted.
  Re-testing musica#956 Case 4 against the fix belongs in the re-measurement
  task.
- **Merge + re-measurement (2026-07-26):** merged into
  `dae-rosenbrock-benchmark` as `19ed87d8` (one conflict, resolved by applying
  the fix inside the branch's Schur else-path; brings upstream #1045/#1046
  along). Full suite on the merged branch: 70/70. All 13 benchmark programs
  rerun clean; everything reproduces the paper at printed precision —
  tropospheric 678/556 with floor 8.7e-3, TS1 463/465/455/462 with errors
  1.3e-2/1.6e-2/1.8e-2 (zero rejections), all QDE table cells and the stiff
  221/218/... sweep, Van der Pol (DAE 1080 steps; 5.3e-2 → 2.1e-7 at first
  order), diurnal envelopes 1.8e-2 / 0.87 / 3.4e-2, Prothero–Robinson orders
  1.02 / 3.03, 5400-case projection sweep with 0 false convergences,
  conservation audit at ulp — **except two numbers, attributed to the fix by
  an in-place A/B (fix reverted → published values reproduce exactly):**
  (1) equilibrium family / Schur ceiling full-ODE accepted steps 503 → 496,
  DAE unchanged at 371 — the published 503 contained rejection-cascade
  fallout; (2) HIRES at rtol 1e-6: max rel. error 4.3e-6 → 1.2e-6 (accepted
  1680 → 1679). Two further deltas are NOT fix-related (identical with the
  fix reverted): the long-time Robertson tight-rtol roundoff-floor cells
  (rtol 1e-8/1e-10 now ODE 2.9e-8/9.7e-8, DAE 4.0e-8/3.1e-8) moved with the
  #1045/#1046 merge and/or toolchain and need re-canonicalizing on the
  paper's rig, and the mechanism-scale per-cell timing ratios (0.85–0.87 vs
  published 0.71–0.87, identical step counts) are timing scatter.

### Phase 0b — Correct and globalize DAE constraint initialization

- **Branches:** `fix/dae-constraint-tolerance-measure` (weighted correction
  and transactional rollback) followed by `dae-init-damping` (affine-covariant
  damped Newton).
- **Scope:** This phase changes only consistent-initial-condition projection.
  The algebraic-variable LTE policy, post-clamp reprojection, Schur-reduced
  stage solves, and pivot-quality diagnostics remain separate work.
- **Status (2026-08-29): complete locally.** The implementation, regression
  tests, benchmark, committed CSV/figure evidence, and three-configuration
  validation are recorded in
  `docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/results.md`.
- **Evidence:**
  - The weighted state-space correction replaced the unit-dependent raw
    residual acceptance test.
  - The follow-on line search uses the simplified Newton correction under the
    already-factored Jacobian, retaining complete-row-scale invariance.
  - At the default ten-update budget, the reduced Henry/dissociation sweep
    converges for all 135 damped starts, including all 15 `XM = 1` cold starts;
    no cell that converges undamped regresses with damping.
  - `constraint_init_max_backtracks_ = 0` retains the prior undamped update,
    and failure remains transactional.
  - The exact-zero common path performs one residual evaluation and no
    Jacobian update, factorization, or solve. Warm starts pay one additional
    residual evaluation and simplified solve per applied Newton update.
  - Current local gates are double 65/65, Kokkos serial 90/90, and float 6/6.
- **MUSICA issue boundary:** The reduced sweep reproduces a temperature-sensitive
  t=0 failure signature inspired by NCAR/musica#956; it is not the reporter's
  configuration. The issue's Case 2 reports `StepSizeTooSmall`, and Case 4 is a
  mid-integration `ConvergenceExceededMaxSteps`. Neither is a demonstrated
  `ConstraintInitializationFailed`, so this phase must not be claimed as a fix
  for either reported case without the original configuration and an A/B run.
- **Remaining upstream work:** submit the line search as a stacked follow-up to
  PR #1083, incorporate current `main` when the stack is ready to advance, and
  let the normal GitHub CI/review process validate supported platforms. Exact
  reproduction of musica#956 remains a separate investigation.

## Performance thread

### Phase 1 — Step-size controller options
- **Branch:** `dae-controller`
- **Evidence:** `rosenbrock.inl:282-287` implements the classic KPP
  controller (`safety/err^(1/(p+1))`, facmin/facmax, no growth after
  rejection, `rejection_factor_decrease_` on double rejection;
  defaults at `rosenbrock_solver_parameters.hpp:30-33`). Dreger et al. 2025
  (GMD 18:4273, `10.5194/gmd-18-4273-2025`) measure 13–43% fewer function
  evaluations in KPP-family Rosenbrock chemistry from Söderlind's H211b
  digital filter or a raised safety factor. The DAE's recommended operating
  point (loose rtol, flat model-error floor) is precisely where controller
  quality dominates cost.
- **Tasks:**
  - Add an opt-in controller policy to `RosenbrockSolverParameters`:
    `{classic (default), h211b, predictive}` with the two state variables
    (previous error, previous H) carried in the solver loop; expose
    `safety_factor_` prominently in docs.
  - H211b: `H_new = H · (tol/err_n)^{1/4·1/p̂} · (tol/err_{n-1})^{1/4·1/p̂} ·
    (H_n/H_{n-1})^{-1/4}` (Söderlind's standard coefficients, adapted to the
    existing `estimator_of_local_order_`); predictive (Gustafsson): the
    RADAU-style accepted-step formula. First accepted step falls back to
    classic.
  - Benchmark matrix: {classic, h211b, predictive} × {Robertson,
    tropospheric, equilibrium N=64, TS1} × rtol {1e-2, 1e-6}, reporting
    accepted/rejected steps, function calls, wall-clock, achieved error.
  - Keep default = classic until the matrix is in; flip the default in a
    separate follow-up decision if the evidence is one-sided.
- **Acceptance:** measurable step/function-call reduction at unchanged
  achieved error on at least the loose-rtol chemistry cases (literature
  expectation 10–40%); bitwise-identical behavior with `classic`; suite
  green.

### Phase 2 — DAE-aware fill-reducing ordering and an honest Schur guard
- **Branch:** `dae-ordering`
- **Evidence:** `solver_builder.inl:94-114` computes the diagonal Markowitz
  reorder from the *kinetic* Jacobian pattern only — before constraint rows
  are merged and before kinetic entries of algebraic-species rows are pruned
  (constraint wiring happens later, `solver_builder.inl:341+`). For TS1 with
  17 dense constraint rows the factored pattern differs materially from the
  one the ordering optimized. Separately, the Schur profitability guard
  (`rosenbrock.inl:88`, `schur_stage_solver.hpp:184-205`) compares raw
  pattern counts (`nnz(S) > nnz(J)`), not post-LU fill of the two
  alternatives, and the paper's "profitable regime may be empty" conclusion
  was measured against the mis-ordered full matrix.
- **Tasks:**
  - Build the reorder input as (kinetic pattern with algebraic-dependent
    rows pruned) ∪ (built-in constraint rows) ∪ (external constraint rows);
    requires resolving constraint sparsity before `GetSpeciesMap()` —
    restructure the build so the constraint row/column sets are available at
    reorder time.
  - Record LU fill (nnz of L+U) before/after for tropospheric-DAE and TS1
    4/9/17-family matrices; report factor/solve flops proxy.
  - Upgrade the Schur guard to compare symbolic post-LU fill:
    `nnz(L+U of S) + Σ nz_g²` vs `nnz(L+U of αM−J)` (both available from the
    existing symbolic machinery), keeping the cheap pattern-count pre-pass as
    a first filter.
  - Re-run `schur_ceiling` and the TS1 Schur-guard row after the ordering
    fix; the paper's structural claim gets either confirmation on a fair
    baseline or a correction.
- **Acceptance:** reordering demonstrably reflects the factored DAE pattern
  (fill no worse than status quo on every benchmark matrix, strictly better
  on TS1); Schur guard decision derived from post-LU fill; re-measured
  Schur/ODE ratios recorded in a note and propagated to the paper if
  changed.

### Phase 3 — Batched constraint evaluation: pow fast paths and vectorized support
- **Branch:** `dae-batch-fastpath`
- **Evidence:** `constraint_set.hpp:478-483` — the batched equilibrium
  residual calls `std::pow(max(0,x), stoich)` unconditionally, including
  `pow(x, 1.0)` for the unit stoichiometries that dominate real equilibria
  (the entire equilibrium family is stoich-1); the Jacobian path already has
  a `stoich == 1.0` fast path, the residual does not. And
  `constraint_set.hpp:57-58`: `kBatchedEvaluation` is disabled for
  vectorized matrix policies, so the production SIMD path still pays the
  per-constraint `std::function` dispatch the scalar batching removed.
- **Tasks:**
  - Residual fast paths: `stoich == 1.0` → multiply, `stoich == 2.0` →
    `x·x`, else `pow`; mirror in the Jacobian partial-product loops.
  - Port the batched loops to vector-ordered dense/sparse policies (the
    arithmetic is identical; only the cell indexing differs — group-major
    with `GroupVectorSize()` inner stride), and enable `kBatchedEvaluation`
    there. Bit-exactness requirement as in the scalar batching commit
    (`775426f3`): identical arithmetic order term-for-term.
  - Re-run `equilibrium_efficiency` (scalar and a new vectorized variant)
    and `constraint_cost_probe`; record the DAE/ODE ratio sweep.
- **Acceptance:** bit-identical results on the existing scalar benchmarks;
  equilibrium-family DAE/ODE ratio at N=256 improves from the current
  0.93–1.40 band (target: ≤1.0 across the sweep if pow was load-bearing —
  measure, don't assume); vectorized DAE path runs the batched loops.

### Phase 4 — Constraint initialization on the algebraic blocks
- **Branch:** `dae-init-blocks`
- **Evidence:** `InitializeConstraints` (`rosenbrock.inl:676-732`) factors
  the full N×N sparse matrix — identity everywhere except the algebraic
  rows — once per Newton sweep, on every `Solve()` call, and again after
  every clamp reprojection (`solver.hpp:289-304`). For TS1 that is a
  227-row sparse factorization to solve a ≤17-variable problem. The
  `SchurStageSolver` symbolic phase (`schur_stage_solver.hpp:110-175`)
  already builds exactly the needed structure: connected components of the
  G_z pattern as small dense blocks with partial pivoting. External
  validation case: NCAR/musica#956's Case 2 — cloud-chemistry DAE failing at
  t = 0 under a 1 K temperature change (Van't Hoff K_eq sensitivity, low-LWC
  row scaling) — fails before any step is attempted, i.e. in exactly this
  routine; the phase's robustness upgrades should be measured against it.
- **Tasks:**
  - Factor out the group-decomposition + dense-LU machinery from
    `SchurStageSolver` into a shared helper; initialization assembles only
    the constraint Jacobian entries (it already zero-fills and writes only
    constraint rows) and solves `G_z δz = −G` per component, per cell.
  - Keep the affine-covariant acceptance logic unchanged (weighted
    correction norm, backtracking, transactional restore) — only the linear
    algebra changes.
  - Pivot diagnostic comes from the dense partial-pivot factors (strictly
    better than the current pivotless-U diagonal read at
    `rosenbrock.inl:691-722`); keep reporting
    `constraint_init_min_pivot_ratio_` with the same semantics.
  - Two free robustness upgrades while in the routine:
    (a) apply the accepted candidate's already-computed simplified-Newton
    correction when it passes the tolerance (it is currently solved and then
    used only as a merit value; Deuflhard's codes apply it);
    (b) clip line-search candidates at zero for algebraic variables so the
    iteration cannot wander to a negative root — nothing currently keeps
    iterates physical, and the batched residual's `max(0,x)` clamp
    (`constraint_set.hpp:479`) makes G/J inconsistent exactly there. Both
    behind the existing parameters struct; (a) on by default, (b) on by
    default for chemistry (documented).
  - Measure init cost: `constraint_cost_probe` plus a new
    per-`Solve()`-cadence micro-benchmark (diurnal-style 288 segments) on
    tropospheric and TS1-17.
- **Acceptance:** identical convergence decisions on the 5400-case
  projection sweep (allowing the two opt-in upgrades to *reduce* failures,
  never add any); init wall-clock at TS1 scale reduced by ≥5× on the
  cadence micro-benchmark; pivoting handles a permuted-coupling G_z block
  that the current pivotless path cannot (new structural test).

### Phase 5 — Stop redoing alpha-independent work on rejections
- **Branch:** `dae-reject-cheap`
- **Evidence:** `SchurStageSolver::Factor` (`schur_stage_solver.hpp:271-320`)
  refactors A_zz and recomputes W = A_zz⁻¹A_zx on every step attempt, but
  both depend only on J, not on alpha; only S needs re-assembly per attempt.
  Symmetrically, the in-place path re-evaluates the full rates+constraints
  Jacobian on every rejection (`rosenbrock.inl:338-349`); a snapshot of the
  unshifted −J vector before the alpha shift turns that into a memcpy — the
  trick the Schur path already uses by keeping `state.jacobian_` intact.
- **Tasks:**
  - Split Schur `Factor` into `FactorAlgebraic(J)` (per Jacobian refresh)
    and `FormAndFactorSchur(alpha)` (per attempt); invalidate the algebraic
    factors whenever the outer loop refreshes `state.jacobian_`.
  - In-place path: snapshot the flat −J vector after assembly, restore on
    rejection instead of re-evaluating rates and constraints (one
    `FlatBlockSize × cells` copy; skip the snapshot when
    `max_number_of_steps_` guarantees no rejection handling is needed — not
    worth a flag if the copy is cheap, measure first).
  - Counters: `jacobian_updates_` must keep meaning "evaluations", so the
    restore path must not increment it.
- **Acceptance:** bit-identical trajectories and counters (minus the removed
  redundant evaluations) on all benchmarks; rejection-heavy stress case
  (loose-rtol Robertson) shows reduced per-rejection cost; suite green.

## Capability thread

### Phase 6 — First-class table-driven QSSA constraint type
- **Branch:** `dae-qssa-constraint`
- **Evidence:** the paper attributes TS1's 1.4–1.6× DAE wall-clock to "the
  price of the generic table-driven constraint callbacks", and the driver
  shows why: `Ts1QssaModel` recomputes `NetStoichD` — a scan over each
  reaction's product list — inside every residual and Jacobian evaluation
  (`benchmark/ts1_dae.cpp:194-203, 288-344`), through per-model
  `std::function` indirection. This is benchmark code standing in for a
  missing library feature: every QSSA user must hand-write an external
  model.
- **Tasks:**
  - Add `QssaConstraint` to `constraint/types/`: constructed from the
    process list plus a set of constrained species; build-time it derives
    the relevant reactions, precomputes net stoichiometries, flattens
    (rate-constant index, reactant indices, net-stoich, Jacobian flat slots)
    into packed arrays; runtime it evaluates in the batched tight-loop style
    of `BatchedEquilibrium`, reading rate constants from the state's rate
    arrays so time-varying photolysis flows through without a bespoke
    parameter path.
  - Reuse the manifold projection already in `InitializeConstraints` (no
    per-model `Project` needed); support per-row scaling as an optional
    conditioning aid, default off (the affine-covariant test made it
    optional).
  - Port `ts1_dae` (and `diurnal_dae`'s constraint) to the library type;
    keep the external-model variants compiled as regression references.
  - Re-run `ts1_dae` 4/9/17 families and the diurnal benchmark; record the
    DAE/ODE wall-clock ratio.
- **Acceptance:** results match the external-model implementation to
  roundoff at matched step sequences; TS1 DAE/ODE wall-clock ratio drops
  materially from 1.4–1.6 (target ≤1.2 — measure); a user can constrain a
  radical family in ~5 lines against a mechanism with no hand-derived
  Jacobian.

### Phase 7 — Runtime conditioning monitor for the algebraic block
- **Branch:** `dae-condition-monitor`
- **Evidence:** the pivot-quality diagnostic exists only inside
  initialization (`rosenbrock.inl:691-722`); during stepping, a degenerating
  G_z surfaces only as NaN or step-size collapse. The paper's diurnal recipe
  switches DAE↔ODE on the sun clock, with the τ_HO2 monitor recorded but not
  acted on (`benchmark/diurnal_dae.cpp:177-230`). Reading the algebraic-row
  U-diagonal ratio after each stage factorization is nearly free — the
  factors already exist (dense-block pivots once Phase 4 lands; U diagonal
  on the sparse path).
- **Tasks:**
  - Per-`Solve()` stat: worst algebraic pivot ratio across step attempts
    (`SolverStats::algebraic_pivot_ratio_`, min-tracked like the init one),
    computed from whichever factorization path ran (sparse U diagonal,
    Schur/dense-block pivots).
  - Optional threshold parameter: when the ratio degrades past it, finish
    the current step and return a named state
    (`AlgebraicConditioningDegraded`) with the partial-interval result and
    `final_time_`, so a host can hand the segment to an ODE solver — the
    guarded diurnal recipe as a signal-driven mechanism instead of a clock.
  - Build-time structural check: verify each constraint row structurally
    depends on its own algebraic species (maximum-transversal / bipartite
    matching on the constraint-row pattern); warn — or with a flag, permute
    the constraint-to-row assignment — instead of letting the pivotless
    sparse LU hit a structural zero pivot on a nonsingular-but-permuted
    block. New structural test alongside the index-2 clean-failure test.
  - Rework `diurnal_dae` to drive the switch from the monitor; compare
    switching times against the sun-clock recipe and the recorded τ_HO2
    signal.
- **Acceptance:** monitor tracks the known conditioning sweep
  (`1e-2 → 1e-10` harness case) during stepping as the init diagnostic does;
  diurnal benchmark reproduces the guarded recipe with monitor-driven
  switches within one segment of the clock-driven ones; zero measurable
  overhead when the threshold is unset (stat only).

### Phase 8 — Non-autonomous stages (scoped)
- **Branch:** `dae-nonautonomous`
- **Evidence:** the tableaus already carry the stage abscissae `alpha_` and
  the d-vector in `gamma_` (`rosenbrock_solver_parameters.hpp:71-74`); the
  stage loop never reads them — rates are frozen per `Solve()` call
  (`rosenbrock.inl:189-204`). The paper's largest errors sit exactly where
  this bites: twilight windows under 300-s frozen-J cadence, where
  "transition errors persist in the slow species".
- **Tasks (scoped to rate interpolation, not general F(t)):**
  - Opt-in per-call rate endpoints: host supplies rate constants at t_n and
    t_n+Δ; the solver evaluates stage rates linearly interpolated at
    `t + alpha_i·H` and adds the `H·gamma_i·F_t` term with
    `F_t = (F(k_end) − F(k_start))/Δ` computed from the same interpolation
    (two extra forcing evaluations per step, or analytic in the linear-rate
    case).
  - Constraint residuals get the same interpolated parameters (the
    time-dependent QSSA constraint currently re-reads frozen photolysis).
  - Validate on Prothero–Robinson in true non-autonomous form (the current
    benchmark autonomizes it) and on the diurnal cycle with 300-s hosting:
    measure the twilight error envelope against the tight-ODE reference vs
    the frozen-J recipe.
  - Explicitly out of scope: arbitrary user F(t) callbacks; the append-t
    trick remains the documented alternative.
- **Acceptance:** design orders reproduced on non-autonomous
  Prothero–Robinson; diurnal twilight error envelope reduced vs frozen-J at
  matched cadence and cost within 1.2× per step; default-off path bitwise
  unchanged.

### Phase 9 — Small robustness and bookkeeping items
- **Branch:** `dae-bookkeeping` (or fold into the nearest phase touching
  each file)
- **Tasks:**
  - `(error < 1) || (H < h_min)` force-accept (`rosenbrock.inl:302`): count
    it (`forced_accepts_`) and surface in `SolverStats`; consider a named
    non-fatal flag on the result — silent accuracy escape hatches should be
    visible.
  - Rejections before the first acceptance are not counted
    (`rosenbrock.inl:334-337`); count them separately or drop the guard,
    so benchmark rejection statistics are complete.
  - Dense output (stretch): RODAS continuous extension so hosts can sample
    without end-of-interval truncation; pairs with `h_persist_`. Only if a
    consumer materializes.
- **Acceptance:** stats additions covered by unit tests; no behavior change
  to accepted trajectories.

## Paper thread (ODE repo, after the phases above)

- **Applied (2026-07-26, paper v6, 33 pp, builds clean):** the RODAS3
  default-tableau sentence in the Results preamble (all adaptive benchmarks
  use `FourStageDifferentialAlgebraicRosenbrockParameters`; convergence
  studies exercise all three); equilibrium full-ODE step count 503 → 496 in
  §6.3 text and the fig:equil caption ("three-quarters the steps" survives:
  371/496 = 0.75); the HIRES rtol=1e-6 cell of tab:testset 4.3e-6 → 1.2e-6;
  and a Reproducibility note documenting the rejected-step correction and
  its exactly-two number changes.
- **Figures deliberately not regenerated in v6:** neither corrected quantity
  appears in any plotted series — the equilibrium/Schur figures plot
  wall-clock, where the A/B attributes no change beyond scatter (the
  published ODE curves timed a 503-step integration, flattering the DAE
  ratio by ~1.4%, well inside the acknowledged scatter). Regenerating from
  this session's venv toolchain would swap the figures' timing provenance
  and force the timing text (0.93–1.40, "1.22 at N=256", the Schur bands)
  to follow. Fold figure regeneration, that timing text, and the long-time
  Robertson tight-rtol cells (toolchain-sensitive roundoff floor, not a fix
  effect) into the next canonical benchmarking session.
- After Phases 1–3: re-run the affected sweeps again (equilibrium-family
  ratios after Phase 3; Schur ceiling/guard rows after Phase 2; loose-rtol
  work–precision points after Phase 1). The Schur "profitable regime may be
  empty" paragraph gets re-examined against the Phase 2 fair baseline.
- After Phase 6: update the TS1 wall-clock column and the "price of the
  generic constraint callbacks" attribution.
- Version-bump per round, per the paper repo's convention.

## Ordering and dependencies

```
Phase 0 ──► everything (re-measurement baseline)
Phase 2 ──► Schur re-measurement; strengthens Phase 4/5 numbers
Phase 3 ──► equilibrium-family claims
Phase 4 ──► Phase 7 (dense-block pivots feed the monitor)
Phase 1, 5, 6, 8, 9 independent after Phase 0
```

Suggested first burst: Phase 0 (with upstream PR), then Phase 1 and Phase 3
in parallel — all three are small, isolated, and immediately measurable with
the existing benchmark harness.
