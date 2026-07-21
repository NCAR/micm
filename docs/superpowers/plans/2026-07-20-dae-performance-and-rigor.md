# DAE Rosenbrock Performance and Rigor Plan

- **Date:** 2026-07-20
- **Branch base:** `dae-rosenbrock-benchmark` (post `d7bf2087` constraint-convergence revision, merged with `main` at `f1165aa8`)
- **Predecessor:** `docs/superpowers/notes/2026-07-15-dae-constraint-convergence-experiments.md` (weighted-correction initialization, validated)
- **Status:** Planning; Phase 1 and Phase 7a in progress on dev branches.

## Purpose

Two threads, driven by the July 2026 benchmark results:

1. **Performance** — the QSSA-DAE always reduces accepted steps but often loses
   wall-clock (Robertson −12–16%, equilibrium triples 1.5–2.9× slower). Remove
   the overheads that are artifacts of implementation rather than of the
   formulation.
2. **Rigor** — extend the evidence base from three in-house systems to
   community-standard test problems, externally validated references, and
   adversarial/structural cases, so the DAE path can be defended as a general
   algorithm rather than a benchmark demonstration.

Every phase lists the measured evidence motivating it, its dev branch, and an
acceptance criterion. Phases are ordered so measurement improvements land
before the optimizations they are needed to judge. Dev branches fork from
`dae-rosenbrock-benchmark` and merge back independently.

## Performance thread

### Phase 1 — Step-size continuity across `Solve()` calls
- **Branch:** `dae-step-persistence`
- **Evidence:** controller_cadence experiment: splitting one interval into
  1000 public `Solve()` calls inflates 9 accepted steps to 9000 because `H`
  restarts from `h_start` every call. Dense-output benchmarks (19–25 output
  times) pay this on every segment, for both ODE and DAE.
- **Tasks:**
  - Opt-in `h_persist_` flag on `RosenbrockSolverParameters`; suggestion field
    on `State` (0 = unset) written on `Converged` exit, consumed at entry.
  - Track the controller's pre-clamp suggestion so the end-of-interval
    truncation of the final step does not poison the carried value.
  - Integration test: k segmented solves ≈ max(k, single-solve steps) accepted
    steps with the flag on; bitwise-identical behavior with the flag off.
  - Re-run cadence experiment and dense-output benchmarks with the flag on;
    record before/after step counts.
- **Acceptance:** cadence 1000-call case drops from ~9000 to ~1000 accepted
  steps (per-call minimum of one step each); default-off behavior unchanged;
  full test suite green.
- **Note:** whether the flag should default on is a follow-up policy decision;
  constraint reinitialization per call stays (cheap under the weighted test).

### Phase 2 — Benchmark measurement rigor
- **Branch:** `dae-work-precision`
- **Evidence:** ODE wall-clock at `j_scale=1` (1276 µs) vs `j_scale=2`
  (684 µs) is measurement noise; accuracy references are self-referential
  (same solver, rtol 1e-10), so a shared systematic error would cancel.
- **Tasks:**
  - Work–precision mode for `robertson_dae` and `tropospheric_dae`: sweep
    rtol 1e-2…1e-10, emit (error, wall-clock, steps) per method — the standard
    stiff-solver presentation.
  - Interleave A/B wall-clock repetitions; raise repetition count; report
    median and spread.
  - External reference solutions: generate with an independent solver
    (scipy `Radau` / SUNDIALS) once per system, check in as data with the
    generation script; benchmarks compare against these, not themselves.
- **Acceptance:** work–precision CSVs and cross-solver reference agreement at
  tight tolerance (differential variables within 10× rtol_ref).

### Phase 3 — Constraint evaluation cost
- **Branch:** `dae-constraint-cost`
- **Evidence:** Robertson DAE: fewer steps, +12–16% wall-clock — pure
  per-step overhead from constraint residual/Jacobian callbacks
  (`std::function` per stage) and re-evaluated constant entries
  (e.g. dG/dA = k1).
- **Tasks:** profile stage cost; devirtualize via the constraint-set policy
  template; cache state-independent Jacobian entries; fuse constraint terms
  into the process-set sweeps where possible.
- **Acceptance:** Robertson DAE wall-clock within ±5% of the full ODE at
  equal steps advantage; no accuracy change (deterministic counters
  identical).

### Phase 4 — Algebraic row elimination (Schur complement)
- **Branch:** `dae-schur-reduction`
- **Evidence:** equilibrium triples: DAE takes 371 vs 503 steps yet runs
  1.5–2.9× slower; the factored system stays 3N with two-thirds algebraic
  rows. Index-1 guarantees ∂G/∂z is nonsingular, so algebraic rows can be
  eliminated from the [αM − J] factorization.
- **Tasks:** prototype reduced-system factorization behind a builder option;
  reuse the sparsity machinery; validate on equilibrium family (largest
  algebraic fraction), Robertson, tropospheric; measure fill/solve cost.
- **Acceptance:** equilibrium DAE/ODE wall-clock ratio < 1 for N ≥ 16 at
  unchanged accuracy; orders retained (re-run fixed-step harness).
- **Risk:** largest engineering item; interacts with vectorized orderings and
  the in-place LU path. Prototype on the standard CPU path first.

### Phase 5 — Norm policy and diagnostics
- **Branch:** `dae-norms-diagnostics`
- **Evidence:** cell-dilution experiment (44 steps / 3.7e-7 error at 1 cell →
  18 / 1.2e-5 at 1000 cells): global WRMS makes per-cell accuracy depend on
  batch composition. Separately, the ε=1e-14 coupled-linear case shows zero
  residual against 2.2e-2 forward error — correction convergence cannot
  certify a numerically singular algebraic block.
- **Tasks:** cellwise-max WRMS as a selectable policy + benchmark both;
  algebraic-block condition/pivot-quality estimate surfaced in
  `SolverStats`/diagnostics with a documented threshold.
- **Acceptance:** cellwise policy holds active-cell error batch-invariant in
  the dilution experiment; singular-block case flags instead of silently
  converging.

### Phase 6 — Stiffly-accurate P-variant tableaus
- **Branch:** `dae-rodas-p`
- **Evidence:** six-stage observed order 3.88/3.94 (x/z) vs design order 4 —
  the classic Prothero–Robinson-type reduction the RODAS-P family targets.
- **Tasks:** add Rodas4P (and Rodas5P if references check out) parameter
  sets; extend the fixed-step harness with a Prothero–Robinson-type stiff
  order test; compare observed orders.
- **Acceptance:** P-variant holds observed order ≥ 3.95 on the stiff order
  test where the current six-stage method degrades.

## Rigor thread

### Phase 7a — Singular-perturbation limit (Van der Pol ε-sweep)
- **Branch:** `dae-vdp-epsilon`
- **Motivation:** QSSA-as-DAE is the ε→0 limit of the fast-species ODEs.
  Demonstrating uniform-in-ε convergence of the stiff ODE to the DAE limit on
  Van der Pol (x' = z, ε z' = (1−x²)z − x) is the canonical way to make that
  rigorous.
- **Tasks:** add a VdP external model (ε > 0: stiff ODE rows; ε = 0:
  algebraic z row) to `dae_convergence_experiments`; sweep
  ε ∈ {1e0…1e-6, 0}; CSV of terminal-state differences ODE(ε) vs DAE and
  step counts; summary lines.
- **Acceptance:** |x_ODE(ε) − x_DAE| → 0 monotonically as ε → 0 (first-order
  in ε), with DAE step count ≈ ε-independent and ODE step count growing.

### Phase 7b — Community test problems with published references
- **Branch:** `dae-ivp-testset`
- **Tasks:** Chemical Akzo Nobel (index-1, chemistry-flavored), HIRES, and
  Robertson in conservation-DAE form (A+B+C = 1 replacing an ODE row,
  t → 1e11) with reference solutions from the IVP test set literature;
  wire into the benchmark suite.
- **Acceptance:** terminal-state agreement with published references at
  standard tolerances; long-time conservation drift ≤ tolerance for the
  conservation form.

### Phase 7c — Time-dependent constraints (diurnal photolysis)
- **Branch:** `dae-diurnal`
- **Tasks:** diurnal J(t) forcing for the tropospheric mechanism (sunrise /
  midday / sunset / night); QSSA validity monitor prototype: flag when a
  constrained species' loss timescale approaches H; compare DAE error
  through twilight against steady-J results.
- **Acceptance:** documented error-vs-time-of-day envelope; monitor flags
  twilight windows where QSSA error exceeds the daytime envelope.

### Phase 7d — Conservation audit
- **Branch:** `dae-conservation`
- **Tasks:** element budgets (N, O, H) tracked over long tropospheric
  integrations for ODE vs DAE; quantify QSSA-induced non-conservation.
- **Acceptance:** drift quantified and reported; if material, evaluate a
  conservation-projection option.

### Phase 7e — Structural and adversarial coverage
- **Branch:** `dae-structural`
- **Tasks:** index-2 guard (constraint row missing its own algebraic
  variable → clean diagnostic, not wrong answers); randomized sparse
  mechanisms with prescribed timescale separation (loop until no new failure
  mode for K rounds); fold the Phase 5 condition diagnostic into these tests.
- **Acceptance:** all structural rejections produce named SolverStates and
  rolled-back state; random sweep runs N ≥ 1e3 mechanisms without silent
  failure.

### Phase 7f — Mechanism scale-up
- **Branch:** `dae-mechanism-scale`
- **Tasks:** apply the QSSA-DAE to an in-repo real mechanism (Chapman first,
  then CB05/TS1-class) with a 10–30-member constrained radical family;
  benchmark with the Phase 2 work–precision rig; this is the configuration
  where multi-mode stiffness removal should outrun per-step overhead.
- **Acceptance:** work–precision curves for full ODE vs DAE on a real
  mechanism; go/no-go recommendation for production use.

## Phase 8 — Paper v2 (DAE analysis repo)

After Phases 1–2 and any completed later phases: add work–precision figures,
the VdP uniform-convergence figure, external-reference validation, and update
the efficiency narrative (step persistence changes both methods' dense-output
numbers). Lives in the `DAE` analysis repository, not this repo.

## Sequencing

- Now: Phase 1 (`dae-step-persistence`) and Phase 7a (`dae-vdp-epsilon`) —
  independent, small, highest evidence-per-effort.
- Next: Phase 2 (measurement rig) before Phases 3–4 so optimization claims
  are credible; Phase 7b alongside.
- Then: Phases 3 → 4 (cost, then dimension), 5–6 as policy/method work,
  7c–7f as the rigor suite matures, Phase 8 last.

## Follow-on results (2026-07-21, branch `dae-constraint-batch`, merged `8efafbc6`)

- **Batched built-in constraint evaluation** (`775426f3`): equilibrium and
  linear constraints compiled into type-packed flat arrays (stoichiometries,
  species ids, pre-resolved Jacobian flat positions) evaluated in single
  tight loops, bit-identical arithmetic, guarded to the non-vectorized CPU
  path. Equilibrium family closed to ODE parity: N=256 DAE/ODE 2.59 -> 1.28
  in the before/after run; merged-branch sweep 0.93-1.40 across N, DAE
  faster than ODE at N=1.
- **Schur re-measurement with batching active**: on the fill-free
  equilibrium matrix, row elimination saves nothing and per-step S assembly
  is exposed (Schur/ODE 0.99-3.6, above plain DAE at every N). Exact
  reduction ceiling unchanged at 0.16-0.42x ODE. Schur stays a guarded
  opt-in for the fill-heavy-J / small-cheap-algebraic-block regime, which
  neither measured extreme realizes.
- **MOZART-TS1 import + scale-up** (`34deacbc`): real TS1 from CAM
  (`pp_trop_strat_mam4_vbsext` chem_mech.in, TS1-fullVBS for CESM2.2) via
  `benchmark/import_ts1.py` -> 227 species + M/O2/N2 fixed, 501 reactions
  (Arrhenius/Troe exact, steady-J photolysis, 3 usr_ rates in JPL form, 46
  usr_/het_ excluded). `benchmark/ts1_dae.cpp`: boundary-layer box 1e-2..1e5 s,
  rtol 1e-6, generic table-driven QSSA external model, families of 4/9/17
  radicals. **Result: NO step advantage at any family size** (ODE 463 acc;
  DAE-4 465, DAE-9 455, DAE-17 462; wall-clock 1.4-1.6x ODE; QSSA model
  error 1.3e-2 -> 1.8e-2 with family width). TS1's steps are
  accuracy-limited on slow species; the L-stable method already absorbs the
  radical timescales. Boundary-of-applicability finding for the paper:
  QSSA-DAE pays where constraints remove step-limiting stiffness.
  Infrastructure holds at scale (17-constraint init + factorization robust).
- **Schur fill guard**: TS1's radicals couple ~200 species so nnz(S) >
  nnz(J); `SchurStageSolver` now counts the S pattern in a cheap pre-pass
  and skips building S/LU when unprofitable; `rosenbrock.inl` falls back to
  the full stage matrix (`SchurNonZeros()` check). TS1 DAE+Schur overhead
  2.9x -> 6%.
- Paper v3 in the ODE repo (`35427af`) carries all of the above.
