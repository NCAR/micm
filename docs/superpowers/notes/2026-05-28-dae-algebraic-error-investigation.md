# DAE Rosenbrock — Algebraic-Variable Error Term Inflates Step Counts

**Date:** 2026-05-28
**Branch:** `dae-rosenbrock-benchmark`
**Status:** Root cause confirmed; fix is a design decision (not yet implemented)
**Context:** While building a QSSA-DAE-vs-full-ODE benchmark on the Robertson problem (see `docs/superpowers/plans/2026-05-28-dae-rosenbrock-benchmark.md`), the DAE solve took far MORE Rosenbrock steps than the full ODE. Investigation showed this is caused by how the DAE solver estimates the error of algebraic variables, not by the chemistry or the benchmark.

---

## 1. Observation

Robertson (`r1: A->B` k1; `r2: 2B->B+C` k2; `r3: B+C->A+C` k3), QSSA reduction makes B algebraic via `G = k1*A - k2*B^2 - k3*B*C = 0`. Both solves use `FourStageDifferentialAlgebraicRosenbrockParameters`, identical rtol=1e-6, atol=rtol*1e-2=1e-8, identical initial condition and output grid. Only difference: the DAE marks B algebraic and adds the constraint.

Initial benchmark result (uniform atol = 1e-8 for all species):

| k2 | full-ODE steps | QSSA-DAE steps | DAE/ODE |
|----|---------------:|---------------:|--------:|
| 3e5 | 507 | 28,719 | 57× |
| 3e6 | 481 | 9,338 | 19× |
| 3e7 | 430 | 3,116 | 7.2× |
| 3e8 | 359 | 1,141 | 3.2× |
| 3e9 | 248 | 495 | 2.0× |

The DAE costs more steps everywhere and the ratio *narrows* with stiffness — the opposite of the naive "DAE removes stiffness → fewer steps" intuition. (That intuition is an *explicit*-solver intuition; Rosenbrock is linearly implicit and already handles the stiff fast mode cheaply — note ODE steps actually *fall* as k2 grows.)

The DAE solution is correct: it matches a tight-tolerance (rtol=1e-10) full-ODE reference to ~1e-8 in the long tail, and the unit test `benchmark/test_robertson_qssa.cpp` passes. So the issue is cost, not correctness.

---

## 2. Root cause

### The mechanism (in the solver code)

`include/micm/solver/rosenbrock.inl:198-213`: for DAE systems, the embedded error estimate for algebraic rows is ~0 (the mass-matrix diagonal `M_ii = 0` zeroes the inter-stage coupling that forms `Yerror`). To give the controller *some* signal, the code overwrites the algebraic-row error with the **raw step change**:

```
Yerror[a] = Ynew[a] - Y[a]      (constraint_set.hpp: SetAlgebraicErrors / BuildAlgebraicErrorFunction)
```

The code comment already flags the hazard:

> "This uses the step change as-is without order scaling. ... The step change is O(H) while the true error is O(H^(p+1)), so this is conservative — the solver may take more steps than strictly necessary."

`NormalizedError` (`rosenbrock.inl:386`) then forms, per variable:

```
errors_over_scale = errors[i_var] / (atol[i_var] + rtol * ymax)
```

RMS-combined across variables; the step controller (`rosenbrock.inl:219-224`) drives `error ≈ 1`.

### Why this throttles the step size

For a **differential** variable, `errors[i]` is the embedded truncation error, `O(H^(p+1))` — vanishes fast as H shrinks, so the controller can take large H while keeping the norm ≈ 1.

For the **algebraic** variable B, `errors[B]` is the *actual change* `ΔB`, which is `O(H)` — it does NOT vanish like a truncation error. The controller therefore picks H so that

```
|ΔB| / (atol_B + rtol*|B|) ≈ 1   →   |ΔB per step| ≈ atol_B   (since atol_B dominates here)
```

i.e., it forces B to move by only ~`atol_B` per step, regardless of how large a step A and C could tolerate. The step count becomes

```
#steps ≈ total_variation(B) / atol_B
```

This is *not* an error control — it conflates "the slaved variable changed a lot" (legitimate, accurate) with "the integration is inaccurate" (false here).

### Quantitative confirmation

- **Magnitude.** At k2=3e7: B0 = sqrt(k1/k2) ≈ 3.65e-5, decaying ∝ sqrt(A) to ~5e-9. Total variation ≈ 3.6e-5. Predicted steps ≈ 3.6e-5 / atol_B(1e-8) ≈ 3.6e3. **Observed: 3116.** (~17%, good for a hand estimate over a log-time grid.)
- **Stiffness scaling.** B ∝ 1/sqrt(k2) ⇒ total_variation(B) ∝ 1/sqrt(k2) ⇒ steps ∝ 1/sqrt(k2). Observed step ratios per 10× k2: 28719→9338 = ×3.07, 9338→3116 = ×3.00 — almost exactly sqrt(10)=3.16. The narrowing-with-stiffness trend is the 1/sqrt(k2) decay of B's range, nothing to do with the ODE.

---

## 3. Decisive experiment (no solver change — public per-variable atol API)

Loosen `atol` for **B only**, leaving atol_A = atol_C = 1e-8. Prediction: DAE step count collapses to ODE levels with A/C accuracy unchanged. (Harness: scratch version of `benchmark/robertson_dae.cpp`, uncommitted; accuracy = max post-transient relative error of A and C vs tight-tol ODE reference, t_skip = 1.0 s.)

### Experiment 1 — across k2, atol_B 1e-8 (tight) vs 1e-2 (loose)

| k2 | ODE steps | DAE steps (atol_B=1e-8) | DAE steps (atol_B=1e-2) | err tight | err loose |
|----|----------:|------------------------:|------------------------:|----------:|----------:|
| 3e5 | 507 | 28,719 | 480 | 1.78e-3 | 1.78e-3 |
| 3e6 | 481 | 9,338 | 461 | 1.32e-3 | 1.32e-3 |
| 3e7 | 430 | 3,116 | 417 | 6.70e-4 | 6.70e-4 |
| 3e8 | 359 | 1,141 | 350 | 2.59e-4 | 2.59e-4 |
| 3e9 | 248 | 495 | 241 | 8.81e-5 | 8.82e-5 |

Loosening atol_B by **six orders of magnitude** drops DAE steps to ≈ ODE levels with **identical A/C accuracy** (differences in the 4th significant figure). The 60× step inflation bought zero accuracy.

### Experiment 2 — DAE steps vs atol_B at fixed k2=3e7

| atol_B | DAE steps | post-transient err (A&C) |
|-------:|----------:|-------------------------:|
| 1e-8 | 3116 | 6.70e-4 |
| 1e-7 | 590 | 6.70e-4 |
| 1e-6 | 420 | 6.70e-4 |
| 1e-5 | 417 | 6.70e-4 |
| 1e-4 | 417 | 6.70e-4 |
| 1e-3 | 417 | 6.70e-4 |
| 1e-2 | 417 | 6.70e-4 |
| 1e-1 | 417 | 6.70e-4 |

Step count is controlled *entirely* by atol_B until it stops binding (~1e-6), then floors at 417 — the genuine A/C-driven cost (≈ the ODE's 430). Accuracy is flat at 6.70e-4 across the whole sweep. The algebraic error term is the sole driver of the excess steps.

---

## 4. Conclusion

The QSSA-DAE does **not** intrinsically lose to the full ODE on Robertson. With a physically appropriate tolerance for the algebraic variable it *matches* the ODE step count at the same accuracy. The apparent blow-up is an artifact of the DAE error controller treating the algebraic variable's raw step change `ΔB = O(H)` as if it were a local truncation error `O(H^(p+1))`, then dividing by a tight absolute tolerance — throttling H to `total_variation(B)/atol_B` steps for no accuracy benefit.

This is a real limitation in the solver's algebraic-error handling, exactly as suspected. It penalizes any algebraic variable whose slaved value is small in magnitude but evolves substantially (the common case for QSS species, trace radicals, partitioning balances).

---

## 5. Candidate fixes (design decision — NOT yet implemented)

The step-change-as-error was introduced deliberately to prevent algebraic *overshoot* (see `test/integration/test_dae_constraint_overshoot.cpp`: an algebraic balance variable being driven negative). So we cannot simply drop algebraic variables from the error norm — that reintroduces the overshoot bug. The fix must keep overshoot protection while not penalizing accurate-but-large slaved evolution.

Options, roughly in order of increasing principled-ness:

1. **Order-scale the algebraic step change.** Multiply `ΔB` by a factor that converts the `O(H)` change into an `O(H^(p+1))`-comparable quantity before norming (e.g., scale by `(H/H_ref)^p`, or use the *second difference* of B across stages as a curvature/error proxy rather than the first difference). Keeps a step-rejection signal proportional to non-smoothness, not to magnitude.

2. **Use the constraint residual as the algebraic error signal.** After the trial step, evaluate `|G(Ynew)|` (how far off the manifold we landed) scaled by an appropriate norm, instead of `ΔB`. A slaved variable that lands on the manifold is accurate by construction; its error is the residual, not its motion. This decouples step size from B's magnitude entirely. (Overshoot then shows up as a large residual or via a separate physicality guard.)

3. **Relative-only / magnitude-aware scaling for algebraic rows.** Norm `ΔB` against `|B|` (pure relative) rather than `atol_B + rtol*|B|`, optionally with a dedicated, looser algebraic atol. Cheap heuristic; removes the atol_B floor problem but is less principled than (1)/(2).

4. **Physicality-only guard for overshoot.** Keep algebraic variables out of the *accuracy* norm, but add an explicit step rejection when an algebraic variable violates a physical bound (e.g., goes negative, or |G| exceeds a constraint tolerance). This separates "accuracy control" (differential vars) from "feasibility control" (algebraic vars) — arguably the cleanest conceptually, and matches how index-1 DAE integrators like RODAS treat algebraic components.

**Recommended next step:** prototype option (2) or (4) behind the existing DAE code path, verify that (a) `test_dae_constraint_overshoot` still passes (overshoot still rejected) and (b) the Robertson DAE step count drops to ODE levels at tight atol without a manual atol_B workaround. Both are quick to check with the harness in this branch.

## 6. Follow-up: is there a "proper local truncation error" to use? (experiment)

Question raised: rather than the step-change hack, why not use a proper local truncation error (LTE) for the algebraic variable? A Rosenbrock method already computes an embedded error `Yerror = Σ e_i k_i` (`rosenbrock.inl:192-196`) — for differential variables that IS the LTE. Is it usable for algebraic rows, or is it really ≈0 as the code comment claims?

Tested three configurations of the algebraic-row error on Robertson (DAE) and on `test_dae_constraint_overshoot` (scratch solver edits, reverted after measuring):

| Config (algebraic-row error) | Robertson DAE steps (k2=3e5..3e9) | DAE accuracy | overshoot test |
|------------------------------|-----------------------------------|--------------|----------------|
| 1. step change `Ynew-Y` (current) | 28719, 9338, 3116, 1141, 495 | 1.78e-3 … 8.8e-5 | PASS |
| 2. embedded `Σ e_i k_i` (override removed) | 480, 461, 417, 350, 241 | 1.78e-3 … 8.8e-5 | PASS |
| 3. forced to exactly 0 | 480, 461, 417, 350, 241 | 1.78e-3 … 8.8e-5 | PASS |

**Findings:**

- **Configs 2 and 3 are identical in every measurable way.** Keeping the embedded estimate is indistinguishable from zeroing it ⇒ the embedded algebraic error is genuinely ≈0 for these DA-Rosenbrock coefficients (the code comment's premise is correct). So there is **no existing LTE to "just use"** for algebraic rows — one would have to be *constructed*.
- **Config 1 (step change) is the sole source of the step inflation.** It is `O(H)`, not an LTE.
- **The override is NOT load-bearing for the current overshoot test:** the test passes in all three configs, including with the algebraic error zeroed (≈ pre-#969 "excluded from norm" behavior). The override and the overshoot test were introduced together in commit `95107e6a` (#969); the test was presumably red without the override then, but a later change (candidates: #993/#994 explicit algebraic species, #982 merge fix) appears to have made the override unnecessary for it. **This must be reconciled before removing the override** — either the test no longer reproduces the original overshoot, or the protection is now provided elsewhere.

### Revised fix options (since the embedded estimate is ≈0)

A. **Re-examine necessity (cheapest, do first).** The override is not needed to pass the overshoot test today. Determine whether it protects any real scenario the test doesn't cover. If not, removing it recovers ODE-level cost for free. If it does, capture that scenario as a new failing test first.

B. **Construct a real LTE for algebraic components:**
   - B1. **Step-doubling / Richardson** on the algebraic component(s): compare one step H against two of H/2; the difference is a genuine `O(H^{p+1})` LTE. Rigorous and general; costs extra solves per step (but far cheaper than the current excess steps).
   - B2. **Stiffly-accurate embedded pair** (RODAS3/4-style coefficients) that yields a nonzero algebraic error estimate. Cheap per step; larger implementation (new method coefficients).
   - B3. **Constraint-residual-based estimate:** derive the algebraic error from the predicted-vs-corrected constraint residual `G` (for index-1 this is proportional to the algebraic LTE). Cheap; reuses existing `G` evaluations.

C. **Feasibility guard for overshoot** (decouples accuracy from feasibility): keep algebraic vars out of the accuracy norm, add an explicit step rejection on physical-bound / constraint-residual violation. Matches how index-1 DAE integrators (RODAS) separate the concerns.

**Recommended sequence:** (A) reconcile the overshoot test — if the override is removable, that alone fixes the cost problem; otherwise (C) a feasibility guard, with (B3) or (B1) if a true accuracy signal for algebraic variables is wanted. All quick to evaluate with the harness in this branch.

## 7. Repro

- Branch `dae-rosenbrock-benchmark`. Build: `cmake -S . -B build -DMICM_ENABLE_BENCHMARKS=ON -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER && cmake --build build --target robertson_dae`.
- Run `./build/robertson_dae` (currently the scratch experiment harness producing the two tables above; uncommitted).
- Solver code: `include/micm/solver/rosenbrock.inl:198-216` (algebraic error injection) and `:363-431` (NormalizedError); `include/micm/constraint/constraint_set.hpp` (`SetAlgebraicErrors`, `BuildAlgebraicErrorFunction`).
