# Damped Newton for DAE Constraint Initialization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an affine-covariant backtracking line search to `InitializeConstraints` so DAE consistent-initial-condition projection converges from cold starts (algebraic guesses far below the root), and produce the CSV/figure evidence that demonstrates it.

**Architecture:** PR #1083 already replaced the initialization acceptance rule with a weighted-correction test but still applies a full, undamped Newton step. This plan globalizes that step: each Newton update becomes a backtracking search over step length, whose merit is the weighted norm of the *simplified* Newton correction at the candidate — solved with the Jacobian already factored at the current iterate, so no refactorization and no row-scale dependence. The search reuses PR #1083's existing reduction closures and adds exactly one new `Function` object. It is opt-out: `constraint_init_max_backtracks_ = 0` executes PR #1083's original statements verbatim and `continue`s.

**Tech Stack:** C++20, CMake, GoogleTest, Kokkos (serial + CUDA), header-only templated solver (`include/micm/solver/rosenbrock.inl`).

**Spec:** `docs/superpowers/specs/2026-08-28-dae-init-cold-start-spec.md`

> **Provenance.** The design below was not merely reasoned about. It was applied to a scratch checkout of `fix/dae-constraint-tolerance-measure` and compiled and run: the 12 existing `constraint_initialization` tests pass unmodified, the six constraint/DAE test binaries pass in double and float, the knife-edge reproducer was run at `%.17e` for warm and cold starts, and `VectorMatrix<Real,4>` multi-cell and singular-Jacobian probes were exercised. Measured numbers appear inline. Three plausible-looking alternatives were tried and rejected on evidence — see Design Rulings; one of them silently breaks the singular-Jacobian case.

> **Execution outcome (2026-08-29).** All tasks below are complete. The final
> audit added an exact-zero short-circuit regression, strengthened the singular
> rollback test to require `InfDetected` and both variables unchanged, corrected
> the D2 benchmark to start from the complete consistent algebraic vector, and
> reconciled the spec, roadmap, and results note. Fresh final gates are double
> 65/65, Kokkos serial 90/90, and float 6/6. D1 measures 135/135 damped
> convergence with zero undamped-to-damped regressions; D2 measures one residual
> evaluation and no matrix work on the exactly consistent path.

---

## Global Constraints

- **Target branch:** work on a local branch `dae-init-damping` cut from `fix/dae-constraint-tolerance-measure` (NCAR/micm PR #1083). The PR-shape decision (fold into #1083 vs stacked PR) is **deferred until D1–D4 evidence exists** — see Task 9. Keep commits atomic and rebasable so either choice stays open.
- **The working tree is checked out on `dae-rosenbrock-benchmark`, which is NOT the target.** That branch is the *reference*: it holds a host-only version of this line search whose loop bound, early return and pivot-ratio diagnostic must all **not** be copied. Read target-branch files with `git show fix/dae-constraint-tolerance-measure:<path>`.
- **Iteration budget is unchanged.** `constraint_init_max_iterations_` stays `10` and keeps counting Newton updates. Do not add a forcing-evaluation cap. Cost is *measured* in D2, not capped. (Final measurement: warm-start init goes from 1.125 us to 1.667 us; see Task 8.)
- **New parameter defaults, exactly:** `constraint_init_max_backtracks_` = `24` (`Index`), `constraint_init_backtrack_factor_` = `0.5` (`Real`), `constraint_init_sufficient_decrease_` = `1e-4` (`Real`).
- **`constraint_init_max_backtracks_ = 0` must reproduce PR #1083 bit-for-bit** (spec R3).
- **No new dense-matrix allocation** (spec R5). Reuse `initial_forcing_` (delta), `K_[0]` (rollback), `Ynew_` (candidate), `Yerror_` (candidate correction). One new 8-byte `Scalar<Real>` member is permitted and required — see Task 3.
- **Every failure exit routes through `restore_and_return`; `Y` is not written until a candidate is accepted** (spec R6).
- **Three build configurations must stay green** (spec R4): double (65 ctest), Kokkos serial (90 ctest), float (6 constraint/DAE ctest targets). No CI job sets `MICM_USE_DOUBLE=OFF`; the float gate is local-only and must be run by hand.
- **Types:** use `micm::Index`, `micm::Real`, `micm::Bool` — never `std::size_t`/`double`/`bool`. The reference branch uses the raw types because it predates `micm/util/types.hpp`; the target branch does not. Write float-safe literals as `micm::Real(1e-4)` and the decrease factor as `Real{ 1.0 } - ...`.
- **`.clang-format` is Google style, `ColumnLimit: 125`, `IndentWidth: 2`, `BreakBeforeBraces: Allman.`** No added line may exceed 125 columns.
- **Out of scope** (spec §4), do not touch: the algebraic-variable LTE fix / `SetAlgebraicErrors`, Schur-reduced stage solves, `ReprojectAfterClamp`, and the `constraint_init_min_pivot_ratio_` diagnostic. The reference implementation contains the pivot-ratio block — **do not port it**; the field does not exist in the target's `SolverStats` and its host matrix indexing does not survive device abstraction.

---

## Verified Ground Truth

Facts established by reading and running the target branch. Trust these over the reference branch and over `DAE.md`.

| Fact | Evidence |
|---|---|
| `InitializeConstraints` is `rosenbrock.inl:404-677`; declaration `rosenbrock.hpp:92-99`; called at `rosenbrock.inl:67` under `if (has_constraints)` (`:38,47`), before the integration loop at `:80` | target branch |
| Newton loop is `for (update = 0; update <= constraint_init_max_iterations_; ++update)`; step 8 returns `Converged` *without* applying `delta`; step 9 is the unconditional `apply_update(Y, delta)` at `:649-672` | target branch |
| Step 8 already reduces `‖delta‖_w` at `Y` into `max_correction` and copies it to host — the sufficient-decrease baseline is **free** | `rosenbrock.inl:634-637` |
| `RosenbrockSolverParameters` has exactly two `constraint_init_*` members, lines 43-44, types `Index` / `Real`. Default ctor is **private** (`:90`); never brace-initialized anywhere, so adding members is source-compatible | `rosenbrock_solver_parameters.hpp` |
| `Print()` prints **no** `constraint_init_*` field and has **zero call sites** repo-wide | `rosenbrock_solver_parameters.hpp:93-150` |
| `RosenbrockTemporaryVariables` has `Ynew_` (:20), `Yerror_` (:23) and six `Scalar<>` members (:24-29). `max_correction_` is the only member #1083 added | `rosenbrock_temporary_variables.hpp` |
| `Ynew_` / `Yerror_` are provably dead during init: first touched by the full overwrites `Ynew.Copy(Y)` at `:151`,`:197` and `Yerror.Fill(0)` at `:203`. `K_[0]` at `:64-65` is a **shape-only** `Function()` build argument whose contents are never read | target branch |
| A pre-built `Function` captures **only** `std::vector<Index> num_cols` plus the closure, and re-validates row/column counts at invocation — so it may be invoked with different same-shaped matrices. `Ynew_`, `Yerror_`, `initial_forcing_`, `K_[i]`, `state.variables_` are all `(grid_cells, species)`. Precedent: `mass_coupling` built with `(K[0],K[0])` at `:64-65`, invoked with `(K[stage],K[j])` at `:183` | `matrix.hpp:898-983`; `kokkos_dense_matrix.hpp:827-881` |
| `MICM_LAMBDA` is `[=]` / `KOKKOS_LAMBDA`, so a plain `Real` is **frozen at Function-construction time**. A changing host scalar must be a `DenseMatrixPolicy::ScalarType<Real>` handle (`shared_ptr` on host, paired `Kokkos::View` on Kokkos) | `types.hpp:14,18`; `scalar_view.hpp:22-33`; `kokkos_scalar_view.hpp:15-16,42-46` |
| `std::isnan`/`isinf`/`isfinite`/`abs` are legal inside `MICM_LAMBDA` device lambdas | `rosenbrock.inl:443,446,453`; `backward_euler.inl:227` |
| **Repeated `Solve` against one factorization is safe.** All implementations bind the factors `const` and read them through const block views; only `Factor` mutates. Nothing between step 5's `Factor` and the last merit `Solve` writes the factor matrices | `linear_solver_in_place.inl:105-160`; `linear_solver.inl:236-243` |
| `Solve` is **in/out on its first argument**, so every trial must redo `candidate_correction.Fill(0)` + `AddForcingTerms` before solving. `Fill(0)` is required because `AddForcingTerms` never writes ODE rows, and external models are contractually permitted to accumulate | `constraint_set.hpp:233` |
| `Axpy` cannot be used for candidate construction: it touches **every** column and cannot express the algebraic-column mask | `matrix.hpp:367-375`; `kokkos_dense_matrix.hpp:466-485` |
| **Native CUDA rejects constrained systems at build time**, throwing `MICM_SOLVER_ERROR_CODE_CUDA_CONSTRAINTS_UNSUPPORTED`. The only device backend that ever executes this code is **Kokkos**; `CudaDenseMatrix` inherits the *host* `Function`/`ScalarType` | `solver_builder.inl:251-266`; `cuda_dense_matrix.hpp:74` |
| `test_constraint_initialization.cpp` is 594 lines, 12 `TEST`s, no fixtures, links `GTest::gtest_main`. `BuildSquaredConstraintSolver` (lines 429-453) gives `2*[B] - [C]^2 = 0`, root `C = 1` at `B = 0.5` | target branch |
| Existing assertions `EXPECT_EQ(loose_stats.constraint_init_iterations_, 1)` (:575) and `EXPECT_EQ(tight_stats.constraint_init_iterations_, 2)` (:593) pin the iteration count exactly | target branch |
| Target branch has **no top-level `benchmark/`**. `MICM_ENABLE_BENCHMARK` (singular) gates `test/benchmark/micm_bench` and is used by four CI workflows — **do not repurpose it** | target `CMakeLists.txt:36-38`; `test/CMakeLists.txt:10-12` |
| **No plotting script exists anywhere in this repo, on any branch, in any commit. No benchmark figure and no benchmark output CSV is committed anywhere** | repo-wide search |
| **`DAE.md` is untracked** (`?? DAE.md` in `git status`) and does not exist on the target branch or on `main` | working tree |

### Measured results from the trial implementation

| Measurement | Result |
|---|---|
| Cold start `XM = 1`, 15 temperatures, knife-edge reproducer | **15/15 `Converged`**, landing on the same `AQ` column as the warm start bit for bit. Baseline: 15/15 `ConstraintInitializationFailed` |
| Warm start `XM = 900` `AQ` column vs PR #1083, `%.17e` under `-ffp-contract=off` | **zero differences** |
| `max_backtracks_ = 0` vs a PR #1083 build | identical state, value, `constraint_init_iterations_` and `solves_` across an 11-point basin sweep x 2 tolerance settings x 2 precisions x `Matrix` and `VectorMatrix<L=4>` |
| Existing tests | **all 12 pass unmodified**, including the exact-count assertions at `:575` and `:593` |
| Constraint/DAE binaries | double and float both green |
| Warm-start init cost, median of 200 | #1083 `1.958 us` -> backtracks=0 `2.542 us` -> backtracks=24 `3.330 us` |
| Hard cold start at `rtol=1e-8, atol=1e-12` (`2B - C^2`, `C0 = 1e-6`) | `Converged`, `C = 1.00000000`, 5 iterations / 28 solves. Baseline fails |

### Design rulings

Three independent designs were produced and disagreed. The rulings below were settled by **running the alternatives**, not by argument. Each records why, so a reviewer can re-open one on evidence rather than taste.

1. **`max_backtracks_ == 0` is an explicit `if` branch containing PR #1083's original step-9 statements verbatim, ending in `continue;`.** Rejected: the reference's `for (backtrack = 0; backtrack <= max_backtracks_; ++backtrack)` bound — under it, `0` still runs one merit-tested trial that **can be rejected**, returning `ConstraintInitializationFailed` where #1083 applied the step unconditionally. Instruction-stream identity is auditable from the diff; the alternatives require floating-point reasoning.
2. **The line search does not `return Converged` early on acceptance.** The reference does, and it would fire on trial 0 of the tight case in `WeightedCorrectionUsesStateTolerances`, reporting `constraint_init_iterations_ == 1` where the test asserts `2`. Fall through instead; the next pass converges through the existing step-2 and step-8 tests.
3. **The accepted candidate is written back with `Y.Copy(candidate)`.** Rejected: `Y.Swap(candidate)` (correct, but aliases `state.variables_` with `Ynew_` mid-solve for a copy that happens at most `max_iterations` times — a reasoning cost with no measurable payoff), and "scale `delta`, then reuse `apply_update`" (which **desynchronises `Y` from the measured candidate under FMA contraction**: `fma(lambda,d,y)` in the candidate kernel vs `y + fl(lambda*d)` in the write-back). `Copy` guarantees `Y` is exactly the state whose merit was accepted. Warm-start bit-identity was then confirmed by measurement, not assumed.
4. **Candidate construction is a 2-arg accumulate — `candidate.Copy(Y)` then `candidate += step * delta` — inside the trial loop.** The kernel is then `apply_update` with one token changed, which is the lowest-risk thing to compile under nvhpc: no new 3-arg instantiation, no hoisting invariant to maintain. The extra `deep_copy` per trial is negligible beside `AddForcingTerms` (a `std::pow` per stoichiometry term) and a triangular solve.
5. **The step-length scalar is a new `Scalar<Real> constraint_init_step_` member.** Rejected: reusing `error_` (provably dead here, but couples two features across a reordering nothing prevents). R5's "no new dense-matrix allocation" enumerates the four *dense workspace* buffers; one `shared_ptr<Real>` is not what it forbids, and PR #1083 itself added `max_correction_` the same way. Rejected: rebuilding the `Function` per backtrack — up to ~250 heap allocations per `Solve()` inside a `noexcept` function.
6. **No second reduction scalar is needed.** Step 8 has already reduced `‖delta(Y)‖_w` into `max_correction` and copied it to host, so `const Real current_correction_norm = max_correction;` costs nothing and frees the scalar for the per-trial merit.

Two alternatives that look like improvements and are **wrong**:

- **Do NOT fold non-finite values into the merit reduction as `+infinity`.** It collapses three finiteness passes into one and looks like a clean device win. It breaks the singular-Jacobian case: an exactly-consistent candidate whose factorization is singular yields `0/0 = NaN` in the simplified solve, which the fold turns into `+inf` -> reject -> backtrack to exhaustion -> **spurious failure**, where the correct answer is `InfDetected` with exact rollback. Measured with `C0 = 0`, `dG/dC = -2C = 0`. Keep the three explicit `check_algebraic_values` passes.
- **Do NOT drop the candidate's exact-zero-residual short circuit.** It is not subsumed by `q <= tolerance`: it is the only guard against `0/0` when the algebraic block is singular, and it costs nothing because `max_residual` is already read back by the finiteness check that precedes it.

### Known risk this plan accepts

Turning damping on changes `max_backtracks_ = 0 -> 1` from "never damp, never fail" to "up to two trials, and **fail** if both are rejected". A case that converges undamped today can therefore fail once damping is enabled. This is standard NLEQ behavior and it is the price of an exact opt-out, but it is a real regression surface. It must be measured (Task 7's basin map is the decisive experiment) and documented in the parameter comment and in `DAE.md` — **not hidden**. If the basin map shows currently-converging cells that start failing, stop and escalate before Task 9.

A second, smaller inefficiency is **known and deliberately not fixed**: once the step shrinks below the roundoff of the iterate, the candidate is bit-identical to `Y`, the merit equals the current norm exactly, and the remaining trials are identical no-ops. A guard (`step_length < 1 && candidate_correction_norm == current_correction_norm` -> break) is a one-line addition that can only save work. It is **not** in the verified implementation, so adding it is an unmeasured change to measured code. Add it only if Task 8's cost numbers show it matters, and re-run the full three-configuration gate if you do.

---

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `include/micm/solver/rosenbrock_solver_parameters.hpp` | Modify (+12) | The three new line-search knobs + five `Print()` lines |
| `include/micm/solver/rosenbrock_temporary_variables.hpp` | Modify (+1) | `constraint_init_step_` scalar handle |
| `include/micm/solver/rosenbrock.inl` | Modify (+194/-17) | One new `Function` (`set_candidate`); the step-9/10 restructure |
| `test/unit/solver/test_constraint_initialization.cpp` | Modify (+275) | Scaffolding + six line-search regression tests |
| `benchmark/dae_init_cold_start.cpp` + `benchmark/CMakeLists.txt` | Create | D1-D2; needs a new `add_executable` — the target branch has no top-level `benchmark/` |
| `benchmark/plot_dae_init_cold_start.py` | Create | The repo's first plotting script |
| `docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/` | Create | Committed CSVs + figures + results note |
| `DAE.md` | Add + correct | D6. **Currently untracked and absent from the target branch** — Task 9 must decide whether it joins this branch or stays local |

No other file changes. `CudaRosenbrockSolverParameters` copy-constructs from the base (`cuda_solver_parameters.hpp:20-23`), there is no Kokkos parameters struct, and `SolverStats` needs no new field.

**Deliberate departure from repo convention, flag it in the PR description:** this repo commits *no* benchmark output CSV and *no* figures anywhere, writing benchmark output to the untracked build tree. Spec acceptance criterion 6 requires D1-D4 "produced as committed CSV and figures". This plan commits them under `docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/` — they are PR evidence for a numerical claim, not regenerable build artifacts — while the benchmark itself still writes to its CWD like every other benchmark here.

---

## Spec Coverage Map

Every spec requirement and deliverable, and the task that discharges it. An executor should be able to close this table before calling the work done.

| Spec item | What it demands | Task |
|---|---|---|
| **R1** Globalized Newton update | backtracking over `lambda`, starting at 1, times `backtrack_factor_` on rejection | Task 5, Step 5 (step-10 loop) |
| **R2** Affine-covariant merit | simplified Newton correction under the **already-factored** Jacobian; accept on `<= tolerance_` OR `< (1 - sigma*lambda) * ‖delta(Y)‖_w` | Task 5, Step 5 (10c/10d); invariance proven in Task 6 |
| **R3** Parameters + exact opt-out | three knobs at 24 / 0.5 / 1e-4; `max_backtracks_ = 0` bit-identical | Task 2 (knobs); Task 4 (anchor test); Task 5, Step 5 (the `continue` branch) |
| **R4** Device-portable | written in #1083's `Function`/`MICM_LAMBDA` idiom; double, Kokkos and float all green | Task 5, Steps 4-5 and 7 |
| **R5** Buffers, no new dense allocation | `initial_forcing_`, `K_[0]`, `Ynew_`, `Yerror_` | Task 5, Step 3 (+ one 8-byte scalar, Task 3) |
| **R6** Transactional failure | every failure through `restore_and_return`; `Y` unwritten until acceptance; non-finite backtracks rather than fails | Task 5, Step 5; asserted by `ExhaustedLineSearchRestoresCallerState` (Task 6) |
| **R7** Cost on the common path | full step accepted first trial; zero-residual short circuit still precedes any factorization; **measured, not assumed** | Task 5, Step 5; measured in Task 7 `RunCost` |
| **D1** Cold-start basin map | CSV + three-panel figure (main / #1083 / #1083+LS) | Task 7 (CSV), Task 8 (figure) |
| **D2** Initialization cost | work counters, warm/consistent, both variants | Task 7 `RunCost`, Task 8 |
| **D3** Row-scale invariance retained | no false convergence or false failure across row scales | Task 6 `BacktrackingIsInvariantToConstraintRowScaling` |
| **D4** Convergence trace | `‖delta‖_w` per update, undamped vs damped | Task 7 scratch instrumentation, Task 8 |
| **D5** Regression tests | in `test/unit/solver/test_constraint_initialization.cpp` | Tasks 2, 4, 5, 6 |
| **D6** Doc corrections | `DAE.md` currently describes parameters the PR does not have | Task 9 |

Spec §7 open questions: **Q1** (PR shape) is deferred to Task 9 Step 6 by decision. **Q2** (iteration budget) is settled — budget unchanged, cost measured in D2. **Q3** (does the reporter's Case 2 actually return `ConstraintInitializationFailed`?) is **not addressed by this plan**; it needs the reporter's config and should be asked on musica#956 in parallel.

---

## Task 1: Baseline — branch and verified-green starting point

Nothing in this tree evidences the spec's claimed green baseline. `build/` holds reference-branch executables and its last recorded failure set is exactly `constraint_initialization`. Establish real numbers before changing anything, or you will not know which failures you caused.

**Files:**
- Create: `docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/baseline.md`

**Interfaces:**
- Produces: three verified test counts that every later task's "still green" check compares against.

- [x] **Step 1: Cut the working branch**

```bash
cd /Users/fillmore/EarthSystem/MICM
git fetch origin
git checkout -b dae-init-damping fix/dae-constraint-tolerance-measure
git log --oneline -1   # expect f72c25df or a descendant
```

- [x] **Step 2: Configure three FRESH build directories**

Do not reuse `build/`, `build-float/`, `build-kokkos/` — they were configured on the reference branch, which has different CMake options and stale executables.

```bash
cmake -S . -B bld-double -DCMAKE_BUILD_TYPE=Release \
  -DMICM_ENABLE_TESTS=ON -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER
cmake -S . -B bld-float  -DCMAKE_BUILD_TYPE=Release \
  -DMICM_ENABLE_TESTS=ON -DMICM_USE_DOUBLE=OFF -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER
cmake -S . -B bld-kokkos -DCMAKE_BUILD_TYPE=Release \
  -DMICM_ENABLE_TESTS=ON -DMICM_ENABLE_KOKKOS=ON -DKokkos_ENABLE_SERIAL=ON \
  -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER
```

`FETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER` is required. Homebrew GoogleTest 1.17.0 is installed and would otherwise be found; tests built against it abort at startup in this environment. Build gtest from source.

- [x] **Step 3: Build and run all three, recording exact counts**

```bash
cmake --build bld-double -j8 && (cd bld-double && ctest --output-on-failure)
cmake --build bld-kokkos -j8 && (cd bld-kokkos && ctest --output-on-failure)
cmake --build bld-float -j8 --target test_constraint_initialization test_equilibrium_constraint \
  test_linear_constraint test_external_model_constraints test_dae_constraint_overshoot \
  test_dae_algebraic_error_insensitivity
(cd bld-float && ctest --output-on-failure -R 'constraint_initialization|equilibrium_constraint|linear_constraint|external_model_constraints|dae_constraint_overshoot|dae_algebraic_error_insensitivity')
```

Expected: `bld-double` 65/65, `bld-kokkos` 90/90, `bld-float` 6/6. Note `-R 'constraint|dae'` would pick up a 7th test (`constraint_set`) and is **not** the R4 set.

- [x] **Step 4: Record the results verbatim**

Write `baseline.md` with the branch SHA, the three configure lines, and the literal ctest summary from each run. **If any configuration is not green, stop and report it** — that is a pre-existing failure the later gates need to know about rather than inherit silently.

- [x] **Step 5: Commit**

```bash
git add docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/baseline.md
git commit -m "Record verified baseline for DAE init line-search work"
```

---

## Task 2: Add the three line-search parameters (inert)

Nothing reads these yet, so behavior is unchanged and the full suite must stay green.

**Files:**
- Modify: `include/micm/solver/rosenbrock_solver_parameters.hpp:43-44` and `:105-106`
- Test: `test/unit/solver/test_constraint_initialization.cpp`

**Interfaces:**
- Produces: `parameters.constraint_init_max_backtracks_` (`Index`, 24), `constraint_init_backtrack_factor_` (`Real`, 0.5), `constraint_init_sufficient_decrease_` (`Real`, 1e-4). Tasks 4 reads all three.

- [x] **Step 1: Write the failing test**

Append to `test/unit/solver/test_constraint_initialization.cpp`:

```cpp
/// @brief The three line-search controls ship with the documented defaults, and 24 backtracks is a
///        step length of 2^-24, below which a correction cannot move a double-precision iterate.
TEST(ConstraintInitialization, LineSearchParameterDefaults)
{
  const auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  EXPECT_EQ(parameters.constraint_init_max_backtracks_, micm::Index(24));
  EXPECT_EQ(parameters.constraint_init_backtrack_factor_, micm::Real(0.5));
  EXPECT_EQ(parameters.constraint_init_sufficient_decrease_, micm::Real(1.0e-4));
}
```

- [x] **Step 2: Run it to verify it fails**

```bash
cmake --build bld-double -j8 --target test_constraint_initialization
```

Expected: **compile error**, `no member named 'constraint_init_max_backtracks_'`.

- [x] **Step 3: Add the members**

**FIND** (lines 43-44):

```cpp
    Index constraint_init_max_iterations_{ 10 };  // max Newton updates for constraint initialization
    Real constraint_init_tolerance_{ 0.1 };       // max weighted Newton correction, as a fraction of the state tolerance
```

**REPLACE WITH:**

```cpp
    Index constraint_init_max_iterations_{ 10 };  // max Newton updates for constraint initialization
    Real constraint_init_tolerance_{ 0.1 };       // max weighted Newton correction, as a fraction of the state tolerance
    // Zero disables the constraint-initialization line search entirely: the full Newton step is then
    // applied unconditionally, reproducing the undamped update exactly. Any non-zero value also makes
    // an exhausted line search a failure, so raising this from 0 to 1 can turn a projection that
    // previously succeeded into ConstraintInitializationFailed.
    Index constraint_init_max_backtracks_{ 24 };        // max line-search reductions per Newton update
    Real constraint_init_backtrack_factor_{ 0.5 };      // line-search step reduction factor, in (0, 1)
    Real constraint_init_sufficient_decrease_{ 1e-4 };  // required fractional decrease in the correction norm
```

`Real x{ 1e-4 }` is well-formed with `Real = float`: floating-to-narrower-floating from an in-range constant expression is not narrowing, and exactness is not required.

- [x] **Step 4: Add the `Print()` lines**

The target prints **none** of the `constraint_init_*` parameters today; adding only the new three would leave the dump half-contradicting the documentation. Add all five. `Print()` has zero call sites, so nothing can regress.

**FIND** (lines 105-106):

```cpp
    std::cout << "h_start: " << h_start_ << std::endl;
    std::cout << "new_function_evaluation: ";
```

**REPLACE WITH:**

```cpp
    std::cout << "h_start: " << h_start_ << std::endl;
    std::cout << "constraint_init_max_iterations: " << constraint_init_max_iterations_ << std::endl;
    std::cout << "constraint_init_tolerance: " << constraint_init_tolerance_ << std::endl;
    std::cout << "constraint_init_max_backtracks: " << constraint_init_max_backtracks_ << std::endl;
    std::cout << "constraint_init_backtrack_factor: " << constraint_init_backtrack_factor_ << std::endl;
    std::cout << "constraint_init_sufficient_decrease: " << constraint_init_sufficient_decrease_ << std::endl;
    std::cout << "new_function_evaluation: ";
```

The pre-existing dangling `std::cout << "absolute_tolerance: ";` at line 149 is **not** part of this change; leave it.

- [x] **Step 5: Run the new test and the full suite**

```bash
cmake --build bld-double -j8 && (cd bld-double && ctest --output-on-failure)
```

Expected: `LineSearchParameterDefaults` PASSES; 65/65 still green.

- [x] **Step 6: Commit**

```bash
git add include/micm/solver/rosenbrock_solver_parameters.hpp test/unit/solver/test_constraint_initialization.cpp
git commit -m "Add constraint-init line-search parameters (inert)"
```

---

## Task 3: Add the step-length scalar temporary (inert)

`MICM_LAMBDA` is `[=]`/`KOKKOS_LAMBDA`, so a plain `Real` named inside a `Function` built once outside the backtrack loop is frozen at construction. The step length must live in a `ScalarType<Real>` handle that the host rewrites and re-uploads per trial — the same mechanism the stage loop uses for `current_c_over_h_`.

**Files:**
- Modify: `include/micm/solver/rosenbrock_temporary_variables.hpp` (after line 29)

**Interfaces:**
- Produces: `derived_class_temporary_variables->constraint_init_step_`, a `Scalar<Real>`. Task 4 writes it per trial and calls `.CopyToDevice()`.

- [x] **Step 1: Add the member**

After line 29 (`Scalar<Bool> inf_detected_;`) add:

```cpp
    // Line-search step length for constraint initialization. This has to be a device-resident
    // scalar handle rather than a plain Real: MICM_LAMBDA captures by value, so a Real named
    // inside a Function built once outside the backtrack loop would be frozen at construction.
    Scalar<Real> constraint_init_step_;
```

It is default-constructed like the other six scalars — do **not** add it to the constructor init-list.

- [x] **Step 2: Note the pre-existing `Clone()` caveat**

All `Scalar<>` members are reference-semantic, and `Clone()` uses the defaulted copy constructor, so a clone *shares* scalar storage while dense members are deep-copied. This is pre-existing for the six existing scalars; a seventh changes nothing. Do not fix it here — out of scope.

- [x] **Step 3: Verify nothing broke**

```bash
cmake --build bld-double -j8 && (cd bld-double && ctest --output-on-failure)
cmake --build bld-kokkos -j8 && (cd bld-kokkos && ctest --output-on-failure)
```

Expected: 65/65 and 90/90. The Kokkos build is the gate — it is where `Scalar<Real>` is `KokkosScalarView` and where default construction of two extent-1 `Kokkos::View`s must work in this position.

- [x] **Step 4: Commit**

```bash
git add include/micm/solver/rosenbrock_temporary_variables.hpp
git commit -m "Add constraint-init line-search step scalar (inert)"
```

---

## Task 4: Test scaffolding and the R3 anchor test

This task adds no production code. It builds the nonlinear test system the line search is judged on, and lands the opt-out equivalence test **before** the behavior change — so that test stands as a witness across the commit that matters. It passes identically before and after Task 5: before because the parameter is ignored, after because `0` disables the search.

The existing `BuildScaledLinearConstraintSolver` is **not** sufficient for this. Its constraint is linear in the algebraic variable, so the simplified Newton correction at the candidate is exactly zero and the line search never engages — a row-scale test built on it would pass against a broken implementation.

**Files:**
- Modify: `test/unit/solver/test_constraint_initialization.cpp` (add 4 includes + scaffolding)

**Interfaces:**
- Produces: `ScaledSquareRootConstraintModel`, `BuildScaledSquareRootSolver(parameters, row_scale)`, `struct ProjectionOutcome { SolverState status_; micm::Real z_; micm::Index iterations_; micm::Index solves_; }`, and `ProjectSquareRoot(parameters, row_scale, z_initial) -> ProjectionOutcome`. Tasks 5 and 6 use all of them.

- [x] **Step 1: Add the required includes**

The scaffolding needs four headers the file does not yet include. Add to the include block:

```cpp
#include <functional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
```

- [x] **Step 2: Add the scaffolding**

Append to `test/unit/solver/test_constraint_initialization.cpp`:

```cpp
// ---------------------------------------------------------------------------
// Damped-Newton (line-search) regression tests
// ---------------------------------------------------------------------------

/// @brief Constraint-only external model enforcing row_scale * ([Z]^2 - [X]) = 0, so the algebraic
///        solution is Z = sqrt(X) for every row scale. The residual is nonlinear in the solved
///        variable, so an undamped Newton step from a guess far below the root overshoots by
///        1/(2*Z_0), which is what the line search exists to control. Multiplying the whole row by
///        row_scale changes G and dG/dy together, so it must change no decision the projection makes.
class ScaledSquareRootConstraintModel
{
 public:
  ScaledSquareRootConstraintModel(std::string x, std::string z, micm::Real row_scale)
      : x_(std::move(x)),
        z_(std::move(z)),
        row_scale_(row_scale)
  {
  }

  std::set<std::string> ConstraintAlgebraicVariableNames() const
  {
    return { z_ };
  }

  std::set<std::string> ConstraintSpeciesDependencies() const
  {
    return { x_, z_ };
  }

  std::set<std::pair<micm::Index, micm::Index>> NonZeroConstraintJacobianElements(
      const std::unordered_map<std::string, micm::Index>& state_indices) const
  {
    const auto i_x = state_indices.at(x_);
    const auto i_z = state_indices.at(z_);
    return { { i_z, i_x }, { i_z, i_z } };
  }

  std::set<std::string> ConstraintStateParameterNames() const
  {
    return {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
  ConstraintUpdateStateParametersFunction(const std::unordered_map<std::string, micm::Index>&) const
  {
    return [](const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&) {};
  }

  template<typename DenseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
      const std::unordered_map<std::string, micm::Index>&,
      const std::unordered_map<std::string, micm::Index>& var) const
  {
    const auto i_x = var.at(x_);
    const auto i_z = var.at(z_);
    const micm::Real scale = row_scale_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
    {
      for (micm::Index i = 0; i < state.NumRows(); ++i)
      {
        forcing[i][i_z] = scale * (state[i][i_z] * state[i][i_z] - state[i][i_x]);
      }
    };
  }

  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
      const std::unordered_map<std::string, micm::Index>&,
      const std::unordered_map<std::string, micm::Index>& var,
      const SparseMatrixPolicy&) const
  {
    const auto i_x = var.at(x_);
    const auto i_z = var.at(z_);
    const micm::Real scale = row_scale_;
    return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
    {
      for (micm::Index i = 0; i < jac.NumberOfBlocks(); ++i)
      {
        jac[i][i_z][i_z] -= 2.0 * scale * state[i][i_z];
        jac[i][i_z][i_x] -= -scale;
      }
    };
  }

 private:
  std::string x_;
  std::string z_;
  micm::Real row_scale_;
};

/// @brief Projection-only system nonlinear in the algebraic variable, with the whole constraint row
///        multiplied by row_scale: row_scale * ([Z]^2 - [X]) = 0, so Z = sqrt(X).
auto BuildScaledSquareRootSolver(const RosenbrockSolverParameters& parameters, micm::Real row_scale)
{
  const auto X = Species("X");
  const auto Z = Species("Z");
  const Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ X, Z } };

  return StandardBuilder(parameters)
      .SetSystem(System(gas_phase))
      .AddExternalModel(ScaledSquareRootConstraintModel("X", "Z", row_scale))
      .SetReorderState(false)
      .Build();
}

/// @brief Drive the projection directly and report status, the converged Z, and the iteration count.
struct ProjectionOutcome
{
  SolverState status_;
  micm::Real z_;
  micm::Index iterations_;
  micm::Index solves_;
};

ProjectionOutcome ProjectSquareRoot(const RosenbrockSolverParameters& parameters, micm::Real row_scale, micm::Real z_initial)
{
  auto solver = BuildScaledSquareRootSolver(parameters, row_scale);
  auto state = solver.GetState(1);

  constexpr micm::Real rtol = std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-4;
  constexpr micm::Real atol = std::is_same_v<micm::Real, double> ? 1.0e-12 : 1.0e-6;
  state.SetRelativeTolerance(rtol);
  state.SetAbsoluteTolerances(std::vector<micm::Real>(state.state_size_, atol));

  state.variables_[0][state.variable_map_.at("X")] = 1.0;
  state.variables_[0][state.variable_map_.at("Z")] = z_initial;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;
  solver.UpdateStateParameters(state);
  state.variables_.CopyToDevice();

  SolverStats stats;
  const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);
  state.variables_.CopyToHost();

  return { status, state.variables_[0][state.variable_map_.at("Z")], stats.constraint_init_iterations_, stats.solves_ };
}
```

- [x] **Step 3: Add the R3 anchor test**

```cpp
/// @brief Disabling the line search restores the undamped iteration exactly, so the change is
///        opt-out and bisectable. Both directions are pinned: the cold start must FAIL with the
///        search off, and the caller's state must come back untouched.
TEST(ConstraintInitialization, DisablingBacktracksReproducesUndampedBehavior)
{
  auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.constraint_init_max_backtracks_ = 0;

  const auto outcome = ProjectSquareRoot(parameters, 1.0, 1.0e-6);

  EXPECT_EQ(outcome.status_, SolverState::ConstraintInitializationFailed);
  EXPECT_EQ(outcome.z_, micm::Real(1.0e-6));  // exact restoration, not approximate
  // No line-search solves were taken: one solve per Newton update and nothing more.
  EXPECT_EQ(outcome.solves_, outcome.iterations_);
}
```

- [x] **Step 4: Run it — it must pass NOW, before any behavior change**

```bash
cmake --build bld-double -j8 --target test_constraint_initialization
./bld-double/test_constraint_initialization --gtest_filter='*DisablingBacktracks*'
```

Expected: PASS. It passes because `constraint_init_max_backtracks_` is currently ignored and the undamped iteration already fails this cold start. If it does **not** pass, the scaffolding does not reproduce the cold-start failure mode and Task 5's gate would be meaningless — fix that before continuing.

- [x] **Step 5: Run the full suite and commit**

```bash
cmake --build bld-double -j8 && (cd bld-double && ctest --output-on-failure)
git add test/unit/solver/test_constraint_initialization.cpp
git commit -m "Add nonlinear projection scaffolding and the line-search opt-out anchor test"
```

---

## Task 5: The damped Newton line search

The single behavior-changing commit. It lands with Task 4's opt-out test already in the tree as a standing witness.

**Files:**
- Modify: `include/micm/solver/rosenbrock.inl` — aliases at `:421-428`, one new `Function` after `:527`, step 9 at `:649-673`
- Test: `test/unit/solver/test_constraint_initialization.cpp`

**Interfaces:**
- Consumes: the three parameters (Task 2), `constraint_init_step_` (Task 3), `ProjectSquareRoot` (Task 4).
- Produces: converging cold starts. `stats.solves_` gains one increment per line-search trial; **no other stat changes**.

- [x] **Step 1: Write the failing test**

```cpp
/// @brief A cold start below the root of a quadratic constraint converges only when the Newton step
///        is damped. From Z0 = 1e-6 toward Z = 1 the full step overshoots by 1/(2*Z0) ~ 5e5, and the
///        undamped iteration spends its whole budget climbing back down. This is the musica#956
///        Case-2 failure mode reduced to two species.
TEST(ConstraintInitialization, BacktrackingConvergesFromColdStart)
{
  const auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  const auto damped = ProjectSquareRoot(parameters, 1.0, 1.0e-6);

  ASSERT_EQ(damped.status_, SolverState::Converged);
  constexpr micm::Real tol = std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-4;
  EXPECT_NEAR(damped.z_, micm::Real(1.0), tol);
  EXPECT_LE(damped.iterations_, parameters.constraint_init_max_iterations_ + micm::Index(1));

  auto undamped_parameters = parameters;
  undamped_parameters.constraint_init_max_backtracks_ = 0;
  const auto undamped = ProjectSquareRoot(undamped_parameters, 1.0, 1.0e-6);
  EXPECT_EQ(undamped.status_, SolverState::ConstraintInitializationFailed);
}
```

- [x] **Step 2: Run it to verify it fails**

```bash
cmake --build bld-double -j8 --target test_constraint_initialization
./bld-double/test_constraint_initialization --gtest_filter='*BacktrackingConvergesFromColdStart*'
```

Expected: **FAIL**, `Actual: ConstraintInitializationFailed`. Do not proceed without a genuinely red test.

- [x] **Step 3: Extend the buffer and scalar aliases**

**FIND** (lines 421-428):

```cpp
    // Reuse initial_forcing_ as the residual/delta workspace, and K_[0] as the rollback
    // buffer: no stage vector is read until the first stage of the integration loop below.
    auto& delta = derived_class_temporary_variables->initial_forcing_;
    auto& original_variables = derived_class_temporary_variables->K_[0];
    auto& max_residual = derived_class_temporary_variables->max_residual_;
    auto& max_correction = derived_class_temporary_variables->max_correction_;
    auto& nan_detected = derived_class_temporary_variables->nan_detected_;
    auto& inf_detected = derived_class_temporary_variables->inf_detected_;
```

**REPLACE WITH:**

```cpp
    // Reuse initial_forcing_ as the residual/delta workspace, K_[0] as the rollback buffer, Ynew_ as
    // the line-search candidate, and Yerror_ as the candidate's simplified Newton correction: no
    // stage vector is read, and neither Ynew_ nor Yerror_ is read, until the first stage of the
    // integration loop below.
    auto& delta = derived_class_temporary_variables->initial_forcing_;
    auto& original_variables = derived_class_temporary_variables->K_[0];
    auto& candidate = derived_class_temporary_variables->Ynew_;
    auto& candidate_correction = derived_class_temporary_variables->Yerror_;
    auto& max_residual = derived_class_temporary_variables->max_residual_;
    auto& max_correction = derived_class_temporary_variables->max_correction_;
    auto& nan_detected = derived_class_temporary_variables->nan_detected_;
    auto& inf_detected = derived_class_temporary_variables->inf_detected_;
    auto& step = derived_class_temporary_variables->constraint_init_step_;
```

- [x] **Step 4: Add the one new `Function` object**

Only one is needed. The existing `check_algebraic_values` and `check_weighted_correction` are reused on the candidate buffers: their bodies name only `diagonal`, the reduction scalars, `atol`, `rtol` and `n_vars` — never the matrices they were built with — and `Function`'s build phase captures only column counts, which match because every one of these buffers is `(grid_cells, species)`.

**FIND** (the tail of `apply_update`, lines 519-527):

```cpp
                  y_view.GetColumnView(i_var),
                  delta_view.GetConstColumnView(i_var));
            }
          }
        },
        Y,
        delta);

    // Projection is transactional:
```

**REPLACE WITH:**

```cpp
                  y_view.GetColumnView(i_var),
                  delta_view.GetConstColumnView(i_var));
            }
          }
        },
        Y,
        delta);

    // Line-search candidate: candidate += step * delta on the algebraic rows only. Identical to
    // apply_update above apart from the step factor, so a damped candidate and an undamped update
    // treat padded rows and non-finite corrections the same way. The step length has to be a
    // device-resident scalar rather than a captured Real: MICM_LAMBDA captures by value, so a Real
    // named here would be frozen at the value it held when this Function object was built. Same
    // mechanism the stage loop uses for current_c_over_h_.
    auto set_candidate = DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType& candidate_view,
            const typename DenseMatrixPolicy::ConstViewType& delta_view) {
          for (Index i_var = 0; i_var < diagonal.size(); ++i_var)
          {
            if (diagonal[i_var] == 0.0)
            {
              candidate_view.ForEachRow(
                  [&](Real& c_val, const Real& d_val)
                  {
                    if (!std::isnan(d_val) && !std::isinf(d_val))
                    {
                      c_val += step * d_val;
                    }
                  },
                  candidate_view.GetColumnView(i_var),
                  delta_view.GetConstColumnView(i_var));
            }
          }
        },
        Y,
        delta);

    // Projection is transactional:
```

Note it is **built** with `(Y, delta)` for column counts and **invoked** with `(candidate, delta)`.

- [x] **Step 5: Replace step 9 with step 9 + step 10**

**FIND** (lines 649-673, the whole of step 9 plus the Newton loop's closing brace):

```cpp
      // 9. Apply update only to algebraic variables, snapshotting the caller's state first
      if (!variables_modified)
      {
        original_variables.Copy(Y);
        variables_modified = true;
      }
      nan_detected = false;
      inf_detected = false;
      max_residual.CopyToDevice();
      nan_detected.CopyToDevice();
      inf_detected.CopyToDevice();
      apply_update(Y, delta);
      max_residual.CopyToHost();
      nan_detected.CopyToHost();
      inf_detected.CopyToHost();

      if (nan_detected)
      {
        return restore_and_return(SolverState::NaNDetected);
      }
      if (inf_detected)
      {
        return restore_and_return(SolverState::InfDetected);
      }
    }
```

**REPLACE WITH:**

```cpp
      // 9. Apply update only to algebraic variables, snapshotting the caller's state first
      if (!variables_modified)
      {
        original_variables.Copy(Y);
        variables_modified = true;
      }

      if (parameters.constraint_init_max_backtracks_ == 0)
      {
        // Line search disabled. Every statement in this block is unchanged from the undamped
        // implementation, in the same order on the same data, so this path is bit-for-bit identical.
        nan_detected = false;
        inf_detected = false;
        max_residual.CopyToDevice();
        nan_detected.CopyToDevice();
        inf_detected.CopyToDevice();
        apply_update(Y, delta);
        max_residual.CopyToHost();
        nan_detected.CopyToHost();
        inf_detected.CopyToHost();

        if (nan_detected)
        {
          return restore_and_return(SolverState::NaNDetected);
        }
        if (inf_detected)
        {
          return restore_and_return(SolverState::InfDetected);
        }
        continue;
      }

      // 10. Damped update: a backtracking line search on the step length. The merit is the weighted
      //     norm of the simplified Newton correction at the candidate -- the correction the Jacobian
      //     already factored at step 5 produces for the candidate's residual. Scaling a complete
      //     constraint row scales G and dG/dy together, so it leaves both corrections and therefore
      //     every decision below unchanged, exactly as it leaves the step-8 test unchanged. Reusing
      //     the factorization keeps a trial to one constraint forcing evaluation and one triangular
      //     solve.
      //
      //     Step 8 has already reduced ||delta(Y)||_w into max_correction, so the sufficient-decrease
      //     reference value costs no extra kernel and no extra device round trip. It is strictly
      //     greater than constraint_init_tolerance_ here, which keeps the threshold below away from
      //     denormals in any precision.
      const Real current_correction_norm = max_correction;

      Real step_length = 1.0;
      bool update_accepted = false;
      for (Index backtrack = 0; backtrack <= parameters.constraint_init_max_backtracks_; ++backtrack)
      {
        // 10a. Candidate iterate y + step * delta, algebraic rows only. Step 7 established that
        //      delta is finite there, so a non-finite candidate means the sum overflowed. Reject it
        //      here rather than downstream: the equilibrium residual clamps negative concentrations
        //      to zero and would report a finite residual for an infinite candidate.
        step = step_length;
        step.CopyToDevice();
        candidate.Copy(Y);
        set_candidate(candidate, delta);

        max_residual = 0;
        nan_detected = false;
        inf_detected = false;
        max_residual.CopyToDevice();
        nan_detected.CopyToDevice();
        inf_detected.CopyToDevice();
        check_algebraic_values(candidate);
        max_residual.CopyToHost();
        nan_detected.CopyToHost();
        inf_detected.CopyToHost();

        if (nan_detected || inf_detected)
        {
          step_length *= parameters.constraint_init_backtrack_factor_;
          continue;
        }

        // 10b. Candidate residual G(y + step * delta)
        candidate_correction.Fill(0);
        constraints_.AddForcingTerms(candidate, state.custom_rate_parameters_, candidate_correction);

        max_residual = 0;
        nan_detected = false;
        inf_detected = false;
        max_residual.CopyToDevice();
        nan_detected.CopyToDevice();
        inf_detected.CopyToDevice();
        check_algebraic_values(candidate_correction);
        max_residual.CopyToHost();
        nan_detected.CopyToHost();
        inf_detected.CopyToHost();

        if (nan_detected || inf_detected)
        {
          step_length *= parameters.constraint_init_backtrack_factor_;
          continue;
        }

        // An exactly zero candidate residual is already consistent, by the same rule step 2 applies
        // at Y. Accepting it here also keeps a singular factorization from turning an exactly zero
        // right-hand side into 0/0 in the solve below.
        if (max_residual == 0.0)
        {
          update_accepted = true;
          break;
        }

        // 10c. Simplified Newton correction, reusing the factorization from step 5. Solve reads the
        //      factors through a const view and never writes them, so one factorization serves every
        //      trial.
        if constexpr (LinearSolverInPlaceConcept<LinearSolverPolicy, DenseMatrixPolicy, SparseMatrixPolicy>)
        {
          linear_solver_.Solve(candidate_correction, state.jacobian_);
        }
        else
        {
          linear_solver_.Solve(candidate_correction, state.lower_matrix_, state.upper_matrix_);
        }
        stats.solves_ += 1;

        max_residual = 0;
        nan_detected = false;
        inf_detected = false;
        max_residual.CopyToDevice();
        nan_detected.CopyToDevice();
        inf_detected.CopyToDevice();
        check_algebraic_values(candidate_correction);
        max_residual.CopyToHost();
        nan_detected.CopyToHost();
        inf_detected.CopyToHost();

        if (nan_detected || inf_detected)
        {
          step_length *= parameters.constraint_init_backtrack_factor_;
          continue;
        }

        // 10d. Accept when the candidate's remaining correction is converged, or when it decreases
        //      the correction at Y by the required fraction.
        max_correction = 0;
        max_correction.CopyToDevice();
        check_weighted_correction(candidate, candidate_correction);
        max_correction.CopyToHost();

        const Real candidate_correction_norm = max_correction;
        const Real required_norm =
            (Real{ 1.0 } - parameters.constraint_init_sufficient_decrease_ * step_length) * current_correction_norm;
        if (candidate_correction_norm <= parameters.constraint_init_tolerance_ ||
            candidate_correction_norm < required_norm)
        {
          update_accepted = true;
          break;
        }

        step_length *= parameters.constraint_init_backtrack_factor_;
      }

      if (!update_accepted)
      {
        return restore_and_return(SolverState::ConstraintInitializationFailed);
      }

      // Take the accepted candidate whole, so the state carried into the next update is exactly the
      // one whose merit was measured. Convergence is deliberately not reported from inside the line
      // search: the next pass re-evaluates the residual at the new state and decides through the
      // step-2 and step-8 tests, which keeps the meaning of constraint_init_iterations_ unchanged.
      Y.Copy(candidate);
    }
```

**Everything else in the function is untouched:** the prologue, the three existing `Function` bodies, `restore_and_return` (530-538), the loop bound `update <= constraint_init_max_iterations_` (543) with its "one pass beyond" comment, the step-2 zero-residual short circuit (574-577), steps 3-8, the extra-pass `break` (644-647), and the terminal `return restore_and_return(...)` (676).

Three details that look redundant and are not:
- The **three separate `check_algebraic_values` passes** per trial (candidate, residual, correction). Folding them into the merit reduction as `+infinity` breaks the singular-Jacobian case — see Design Rulings.
- The **exact-zero candidate residual short circuit**. It is the guard against `0/0` when the algebraic block is singular, and it is free.
- **`candidate_correction.Fill(0)` every trial.** `Solve` is in/out on its first argument, and `AddForcingTerms` never writes ODE rows.

- [x] **Step 6: Run the cold-start test**

```bash
cmake --build bld-double -j8 --target test_constraint_initialization
./bld-double/test_constraint_initialization --gtest_filter='*BacktrackingConvergesFromColdStart*'
```

Expected: PASS.

- [x] **Step 7: Run the full suite in all three configurations**

```bash
cmake --build bld-double -j8 && (cd bld-double && ctest --output-on-failure)
cmake --build bld-kokkos -j8 && (cd bld-kokkos && ctest --output-on-failure)
cmake --build bld-float -j8 --target test_constraint_initialization test_equilibrium_constraint \
  test_linear_constraint test_external_model_constraints test_dae_constraint_overshoot \
  test_dae_algebraic_error_insensitivity
(cd bld-float && ctest --output-on-failure -R 'constraint_initialization|equilibrium_constraint|linear_constraint|external_model_constraints|dae_constraint_overshoot|dae_algebraic_error_insensitivity')
```

Expected: 65/65, 90/90, 6/6, matching Task 1's baseline. **No existing test may be edited to make this pass** — all 12 passed unmodified in the trial implementation.

Watch `EXPECT_EQ(loose_stats.constraint_init_iterations_, 1)` (:575) and `EXPECT_EQ(tight_stats.constraint_init_iterations_, 2)` (:593). If either fails, the line search is returning `Converged` early — re-read Design Ruling 2. A green double build proves nothing about device portability; the Kokkos build is what exercises the `KokkosScalarView` capture.

- [x] **Step 8: Commit**

```bash
git add include/micm/solver/rosenbrock.inl test/unit/solver/test_constraint_initialization.cpp
git commit -m "Damped Newton line search for DAE constraint initialization

Globalizes InitializeConstraints with an affine-covariant backtracking
line search. The merit is the weighted norm of the simplified Newton
correction under the already-factored Jacobian (Deuflhard natural
monotonicity), so row scale still cancels and no trial refactorizes.

constraint_init_max_backtracks_ = 0 executes the previous statements
verbatim and is bit-identical to the undamped implementation."
```

---

## Task 6: Invariance, zero-perturbation and failure-mode coverage (D3, D5)

Four tests that pin the properties the line search must not break. All were verified to compile and pass against the Task 5 implementation.

**Files:**
- Modify: `test/unit/solver/test_constraint_initialization.cpp`

**Interfaces:**
- Consumes: `ProjectSquareRoot` and `BuildScaledSquareRootSolver` (Task 4); the line search (Task 5).

- [x] **Step 1: Add the tests**

```cpp
/// @brief The merit function is the weighted norm of the simplified Newton correction, so scaling a
///        complete constraint row scales G and dG/dy together and cancels. Every decision the line
///        search makes -- which trial is accepted, and therefore how many updates and solves the
///        projection costs -- must be identical across row scales spanning 24 orders of magnitude.
TEST(ConstraintInitialization, BacktrackingIsInvariantToConstraintRowScaling)
{
  const auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  const std::array<micm::Real, 3> row_scales{ 1.0e-12, 1.0, 1.0e12 };
  constexpr micm::Real z_tolerance = std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-4;

  const auto reference = ProjectSquareRoot(parameters, 1.0, 1.0e-6);
  ASSERT_EQ(reference.status_, SolverState::Converged);

  for (const micm::Real row_scale : row_scales)
  {
    const auto scaled = ProjectSquareRoot(parameters, row_scale, 1.0e-6);
    EXPECT_EQ(scaled.status_, SolverState::Converged) << "row scale=" << row_scale;
    EXPECT_EQ(scaled.iterations_, reference.iterations_) << "row scale=" << row_scale;
    EXPECT_EQ(scaled.solves_, reference.solves_) << "row scale=" << row_scale;
    EXPECT_NEAR(scaled.z_, reference.z_, z_tolerance) << "row scale=" << row_scale;
  }
}

/// @brief On a start already inside the full-step convergence basin the first trial satisfies the
///        sufficient-decrease test, so the damped projection walks exactly the undamped iteration:
///        the same number of updates and the same final state, bit for bit. This is what keeps an
///        existing, working initialization unperturbed by the line search.
TEST(ConstraintInitialization, FullStepIsAcceptedInsideTheConvergenceBasin)
{
  auto damped_parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto undamped_parameters = damped_parameters;
  undamped_parameters.constraint_init_max_backtracks_ = 0;

  for (const micm::Real z_initial : { micm::Real(0.5), micm::Real(2.0), micm::Real(10.0) })
  {
    const auto damped = ProjectSquareRoot(damped_parameters, 1.0, z_initial);
    const auto undamped = ProjectSquareRoot(undamped_parameters, 1.0, z_initial);

    EXPECT_EQ(damped.status_, SolverState::Converged) << "Z initial=" << z_initial;
    EXPECT_EQ(undamped.status_, SolverState::Converged) << "Z initial=" << z_initial;
    EXPECT_EQ(damped.z_, undamped.z_) << "Z initial=" << z_initial;
    EXPECT_EQ(damped.iterations_, undamped.iterations_) << "Z initial=" << z_initial;
    // One extra simplified-Newton solve per update buys the globalization
    EXPECT_EQ(damped.solves_, 2 * undamped.solves_ - 1) << "Z initial=" << z_initial;
  }
}

/// @brief When the step-length budget cannot reach the basin, the projection fails transactionally.
///        The caller gets their own state back, not a half-damped iterate.
TEST(ConstraintInitialization, ExhaustedLineSearchRestoresCallerState)
{
  auto parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.constraint_init_max_backtracks_ = 2;  // 0.5^2 cannot climb down from an overshoot of 5e5

  auto solver = BuildScaledSquareRootSolver(parameters, 1.0);
  auto state = solver.GetState(1);
  constexpr micm::Real Z_in = 1.0e-6;
  constexpr micm::Real X_in = 1.0;
  constexpr micm::Real rtol = std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-4;
  constexpr micm::Real atol = std::is_same_v<micm::Real, double> ? 1.0e-12 : 1.0e-6;
  state.SetRelativeTolerance(rtol);
  state.SetAbsoluteTolerances(std::vector<micm::Real>(state.state_size_, atol));
  state.variables_[0][state.variable_map_.at("X")] = X_in;
  state.variables_[0][state.variable_map_.at("Z")] = Z_in;
  state.conditions_[0].temperature_ = 298.15;
  state.conditions_[0].pressure_ = 101325.0;
  solver.UpdateStateParameters(state);
  state.variables_.CopyToDevice();

  SolverStats stats;
  const auto status = solver.solver_.InitializeConstraints(state, parameters, stats);
  state.variables_.CopyToHost();

  ASSERT_EQ(status, SolverState::ConstraintInitializationFailed);
  EXPECT_EQ(state.variables_[0][state.variable_map_.at("Z")], Z_in);
  EXPECT_EQ(state.variables_[0][state.variable_map_.at("X")], X_in);
}
```

`BacktrackingIsInvariantToConstraintRowScaling` asserts **exactly equal** `constraint_init_iterations_` and `solves_` across row scales `{1e-12, 1, 1e12}` — i.e. the search made the identical accept/reject decision at every trial, not merely reached the same answer. Do not weaken it to `EXPECT_NEAR`; that is the whole point of the test.

`FullStepIsAcceptedInsideTheConvergenceBasin` is the zero-perturbation gate: for starts already inside the basin, the final value and iteration count are `EXPECT_EQ`-identical with the search on and off, and the solve count is exactly `2 * undamped - 1` (one merit solve per update, none after the last).

- [x] **Step 2: Add the shape and singular coverage**

Add a `VectorMatrix<micm::Real, 4>` / 3-grid-cell variant of the cold-start test (exercising padded cells and multiple cells at once), and a singular case `Z0 = 0` where `dG/dZ = 0`, asserting a named failure state and **exact** rollback of both variables. Both were exercised in the trial implementation: the multi-cell case converges per cell, and the singular case returns `InfDetected` with exact rollback under both the baseline and the line search.

- [x] **Step 3: Run everything and commit**

```bash
cmake --build bld-double -j8 && (cd bld-double && ctest --output-on-failure)
cmake --build bld-kokkos -j8 && (cd bld-kokkos && ctest --output-on-failure)
git add test/unit/solver/test_constraint_initialization.cpp
git commit -m "Pin row-scale invariance, basin zero-perturbation and failure modes"
```

---

## Task 7: Cold-start basin map benchmark (D1, D2, D4)

The decisive experiment. The committed benchmark emits the D1 and D2 CSVs;
the D4 CSV comes from the scratch instrumentation required below.

**Files:**
- Create: `benchmark/dae_init_cold_start.cpp`, `benchmark/CMakeLists.txt`
- Modify: top-level `CMakeLists.txt`

**Interfaces:**
- Produces: `dae_init_basin.csv` and `dae_init_cost.csv` in the CWD. Scratch
  instrumentation produces `dae_init_trace.csv` without widening `SolverStats`.

- [x] **Step 1: Read the existing reproducer for the physics**

```bash
git show dae-rosenbrock-benchmark:benchmark/musica956_knife_edge_main_api.cpp
```

It already expresses the musica#956 Case-2 chemistry against the templated constraint API the target branch has — Henry's law (Van't Hoff, `K_ref = 0.5`, `delta_H = -60000 J/mol`) feeding a quadratic aqueous dissociation, both as built-in `EquilibriumConstraint`s, `P = 85000 Pa`. Reuse its setup verbatim; do not re-derive the chemistry. **Change its `printf` to `%.17e`** — the committed version prints `%.6e`, which is too coarse to demonstrate bit-identity.

- [x] **Step 2: Add the benchmark plumbing**

The target branch has **no top-level `benchmark/`**, and `MICM_ENABLE_BENCHMARK` (singular) is used by four CI workflows for `test/benchmark/micm_bench` — **do not repurpose it**. Add a separate gate mirroring the reference branch's:

```cmake
option(MICM_ENABLE_DAE_BENCHMARKS "Build the DAE initialization benchmark programs" OFF)
```

and, guarded by `if(PROJECT_IS_TOP_LEVEL AND MICM_ENABLE_DAE_BENCHMARKS)`, `add_subdirectory(benchmark)` containing:

```cmake
set(CMAKE_CXX_CLANG_TIDY "")
add_executable(dae_init_cold_start dae_init_cold_start.cpp)
target_link_libraries(dae_init_cold_start PUBLIC musica::micm)
```

- [x] **Step 3: Write the benchmark**

Follow `conservation_audit.cpp:157-181` for the emission idiom: open the file, write a header, `csv << std::scientific << std::setprecision(17);`, write rows with `<<` and `','`, print `wrote <name>.csv`.

- `RunBasinMap` -> `dae_init_basin.csv`, header `variant,temperature,initial_xm,status,final_aq,constraint_init_iterations,solves`. Sweep `T` over 278..292 K (15 values), `initial_xm` over `{1, 10, 100, 300, 600, 800, 900, 1200, 2000}`, `variant` over `{undamped, damped}` (`max_backtracks_` 0 and 24).
- `RunCost` -> `dae_init_cost.csv`, with the full starting algebraic state,
  status, residual-evaluation, Jacobian, factorization, solve, iteration, and
  timing fields. Warm (`AQ = 1`, `XM = 900`) and independently verified exactly
  consistent starts, both variants, median of 200 runs. This is D2.
- Scratch trace instrumentation -> `dae_init_trace.csv`, header
  `variant,update,weighted_correction_norm,accepted_step`. One cold-start case
  (`T = 284 K`, `XM = 1`), both variants. This is D4.

For the trace you need the per-update norm, which `SolverStats` does not expose. Do **not** add a stats field for a one-off measurement — that widens the reviewed surface. Print the values from a temporary `std::cerr` line in a scratch build, capture them, and **revert the instrumentation before committing**. Record in the results note that the trace came from an instrumented build.

- [x] **Step 4: Build and run**

```bash
cmake -S . -B bld-bench -DCMAKE_BUILD_TYPE=Release -DMICM_ENABLE_DAE_BENCHMARKS=ON \
  -DMICM_ENABLE_TESTS=OFF
cmake --build bld-bench -j8 --target dae_init_cold_start
(cd bld-bench && ./dae_init_cold_start)
```

- [x] **Step 5: Produce the third basin panel from `origin/main`**

`origin/main` has neither the weighted-correction rule nor the line search, so the `constraint_init_max_backtracks_` references will not compile there. Copy the program aside, strip the variant loop, emit `variant=main` rows only, build against `origin/main`, and merge the rows into `dae_init_basin.csv`.

- [x] **Step 6: Check against the spec's predictions and the regression surface**

Spec §1.1 predicts `origin/main` fails at 284/287/290 K warm and everywhere cold; #1083 converges warm everywhere and fails cold everywhere; #1083 + line search converges everywhere. **The cold-start `damped` column must be `Converged` at all 15 temperatures — acceptance criterion 1.** The trial implementation measured exactly that, landing on the same `AQ` column as the warm start bit for bit.

Then scan for any `(T, initial_xm)` cell that is `Converged` under `undamped` and `ConstraintInitializationFailed` under `damped`. **If any such cell exists, stop and escalate** — that is damping making things worse, and the contingency needs a decision before Task 9.

- [x] **Step 7: Commit**

```bash
git add benchmark/ CMakeLists.txt
git commit -m "Add cold-start basin map, cost and trace benchmark"
```

---

## Task 8: Figures and committed evidence

**Files:**
- Create: `benchmark/plot_dae_init_cold_start.py`
- Create: `docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/{*.csv,*.png}`

- [x] **Step 1: Write the plotting script**

No plotting script exists anywhere in this repository, on any branch, in any commit — this is the first, so it defines the convention. Keep it dependency-light (matplotlib + `csv.DictReader`, no pandas), take input and output directories as `argv`, write PNG at 150 dpi.

- `dae_init_basin.png` — the three-panel basin map D1 requires. One panel per variant (`main`, `undamped`, `damped`), `initial_xm` on y (log), temperature on x, cell color by status. Shared axes, one legend.
- `dae_init_cost.png` — grouped bars, D2, warm and consistent starts, both variants.
- `dae_init_trace.png` — D4, `weighted_correction_norm` vs Newton update, log y, one line per variant, accepted step annotated.

- [x] **Step 2: Generate and actually look at the figures**

```bash
python3 benchmark/plot_dae_init_cold_start.py bld-bench docs/superpowers/notes/2026-08-29-dae-init-cold-start-results
```

Open all three. A figure nobody looked at is not evidence. Check specifically that the `damped` cold-start column is uniformly converged and that `main`/`undamped` show the predicted failures.

- [x] **Step 3: Commit the CSVs and figures**

```bash
cp bld-bench/dae_init_*.csv docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/
git add benchmark/plot_dae_init_cold_start.py docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/
git commit -m "Add cold-start basin map, cost and trace evidence

Commits benchmark output, which this repo otherwise does not do: these
are evidence for a numerical claim in review, not regenerable build
artifacts."
```

---

## Task 9: Documentation, final verification, PR-shape decision

**Files:**
- Add/modify: `DAE.md`, `FABLE_PLAN.md`
- Create: `docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/results.md`

- [x] **Step 1: Decide what happens to `DAE.md`**

`DAE.md` is **untracked** and exists on neither the target branch nor `main`. It documents all five initialization parameters as if they shipped together, describing the reference branch's state. Decide explicitly, and say which in the commit: either add it to this branch (now that the line search makes it true) or leave it untracked and local. Do not silently `git add` a file that has never been in the repository without flagging it.

- [x] **Step 2: Reconcile it against the implementation**

```bash
grep -n "constraint_init" DAE.md
```

Check every parameter's name, type, default and description against `rosenbrock_solver_parameters.hpp` as it now stands. Then add the two things it does not say:

- The merit is the weighted norm of the *simplified* Newton correction under the already-factored Jacobian, which is why row scale cancels.
- `constraint_init_max_backtracks_ = 0` is an exact opt-out, **and** `0 -> 1` is a behavior discontinuity: with the search enabled, an update whose trials all fail the sufficient-decrease test is a failure, where an undamped update would have applied the full step regardless.

Also delete anything `DAE.md` claims about post-clamp reprojection and `SetAlgebraicErrors` — both are false on this branch and both are explicit spec non-goals.

- [x] **Step 3: Re-run all three configurations from clean**

```bash
rm -rf bld-double bld-float bld-kokkos
```

then repeat Task 1 Steps 2-3 exactly. Incremental builds of a header-only templated solver are the classic way to ship a change that only compiles against stale objects.

- [x] **Step 4: Confirm every spec acceptance criterion with evidence**

1. Cold start converges at all 15 T — `dae_init_basin.csv`, `variant=damped`, `initial_xm=1`.
2. Warm start converges at all 15 T, post-init `AQ` bit-identical to #1083 — `dae_init_basin.csv` at `%.17e` plus `FullStepIsAcceptedInsideTheConvergenceBasin`.
3. `max_backtracks_ = 0` reproduces #1083 exactly — `DisablingBacktracksReproducesUndampedBehavior`.
4. Full suite green in double, Kokkos, float — Step 3.
5. No change to any constraint-free solve — `PureODESystemUnaffected`, plus the `if (has_constraints)` guard at `rosenbrock.inl:47`.
6. D1-D4 committed as CSV and figures — Task 8.

Any criterion without evidence is **not met**. Report it as unmet rather than arguing it is probably fine.

- [x] **Step 5: Write the results note and update `FABLE_PLAN.md`**

`results.md`: the three-configuration results, the basin table (spec §1.1 extended with the new column), the D2 cost numbers, and an explicit statement that the D4 trace came from an instrumented scratch build whose instrumentation was reverted. Follow the table idiom in `dae-rosenbrock-benchmark:docs/superpowers/notes/2026-07-26-musica956-knife-edge-repro.md`. Mark `FABLE_PLAN.md` Phase 0b complete and link the note.

- [x] **Step 6: Decide the PR shape and report**

Present to the user: the measured cold- and warm-start results, the D2 common-path cost, whether the basin map showed any `undamped`-converges / `damped`-fails cell, and the diff size against `fix/dae-constraint-tolerance-measure`.

Recommend **stacked follow-up PR** unless D2 shows the added cost is negligible *and* the basin map is clean, in which case folding into #1083 gives maintainers one coherent change. **Do not open a PR without the user's decision** — this question was explicitly left open.

- [x] **Step 7: Commit**

```bash
git add DAE.md FABLE_PLAN.md docs/superpowers/notes/2026-08-29-dae-init-cold-start-results/results.md
git commit -m "Correct DAE.md init parameters; record cold-start results"
```
