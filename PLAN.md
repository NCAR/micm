# DAE Rosenbrock Constraint Tolerance and Convergence Experiment Plan

## Purpose

Determine a numerically defensible treatment of algebraic-constraint tolerance and convergence in the MICM Rosenbrock DAE solver.

## Implementation status — 2026-07-15

- Completed the baseline API repair and focused DAE/chemistry validation.
- Added the `dae_convergence_experiments` benchmark target and reproducible CSV output beneath the build tree.
- Implemented the experiment-supported weighted-correction convergence rule, affine-covariant backtracking, final-update convergence check, and failure rollback.
- Added post-clamp algebraic reprojection and focused regression coverage.
- Preserved the existing global WRMS integration-error policy. The cell-dilution experiment quantifies a real batch-size dependence; changing that established controller contract requires a separate policy decision.
- Recorded commands, numerical results, remaining limitations, and interpretation in `docs/superpowers/notes/2026-07-15-dae-constraint-convergence-experiments.md`.

The experiments will answer four separate questions:

1. How should consistent-initial-condition Newton convergence be measured?
2. Do the four-stage and six-stage DAE Rosenbrock methods achieve their expected convergence orders for both differential and algebraic variables?
3. How much constraint drift occurs during integration, and is it controlled by the existing LTE test?
4. Are solver accuracy, feasibility, and work invariant to constraint scaling, output cadence, and the number of algebraic variables?

No production solver behavior should change until the manufactured-problem experiments establish a clear convergence contract.

## Current implementation baseline

- `constraint_init_tolerance_` defaults to `1e-10`.
- Initialization convergence is the unscaled infinity norm of the raw algebraic residual over all cells and algebraic rows:

  ```text
  max |G_i(y)| < constraint_init_tolerance
  ```

- Initialization uses undamped Newton updates, changes only algebraic variables, and allows ten updates by default.
- The residual is checked before each update but is not rechecked after the final permitted update.
- A failed initialization leaves the state at its final attempted iterate.
- Rosenbrock stage solves use linearized constraint equations; accepted steps do not receive a separate nonlinear constraint-residual test.
- Step acceptance uses the embedded LTE norm. Algebraic LTE entries are generally near zero after removal of the old raw-step-change override.
- Algebraic rows still contribute to the RMS norm denominator even when their LTE entries are zero.
- Every public `Solve()` call reinitializes the algebraic variables.
- Post-solve nonnegative clamping changes differential variables without a subsequent constraint projection.

Relevant implementation files:

- `include/micm/solver/rosenbrock.inl`
- `include/micm/solver/rosenbrock_solver_parameters.hpp`
- `include/micm/solver/solver.hpp`
- `include/micm/constraint/constraint_set.hpp`

## Experimental principles

- Use manufactured systems before chemistry mechanisms so exact solutions and conditioning are known.
- Separate initialization convergence, integration accuracy, manifold feasibility, and physical bounds.
- Treat deterministic solver counters as primary performance evidence; wall-clock timings are secondary.
- Always report both the raw constraint residual and a dimensionless or state-error-equivalent measure.
- Repeat key experiments after multiplying complete constraint rows, residual and Jacobian together, by arbitrary constants.
- Preserve exact parameter sets, seeds, compiler configuration, and CSV schemas so results are reproducible.
- Keep experimental instrumentation isolated from production solver behavior until a stopping rule is selected.

## Phase 0: Restore and freeze the baseline

### Tasks

1. Repair API drift in the existing branch-only tests and benchmarks:

   - Replace obsolete `K_HLC_ref` and `delta_H` field names with their current forms.
   - Replace obsolete `SystemParameters` construction with the current `System` API.
   - Update stale comments that still describe the removed algebraic step-change error override.

2. Rebuild and run these baseline targets:

   - `test_constraint_initialization`
   - `test_dae_constraint_overshoot`
   - `test_dae_algebraic_error_insensitivity`
   - `test_dae_algebraic_error_step_economy`
   - `test_robertson_qssa`
   - `robertson_dae`
   - `tropospheric_dae`

3. Record the compiler, build type, MICM commit, and baseline test results.

4. Add experiment-only diagnostics for:

   - Initial and final residual by cell and constraint row
   - Newton correction by algebraic variable
   - Residual-reduction factor by iteration
   - Initialization termination reason
   - Residual immediately after integration and after post-solve clamping

### Exit criterion

All baseline targets build from current sources, existing active regressions pass, and diagnostic output does not change solver results.

## Experiment 1: Constraint-row scaling invariance

### Model

Hold `x` fixed during initialization and solve

```text
G_s(x, z) = s (z - x^2) = 0
```

for the algebraic variable `z`. Scaling both the residual and Jacobian row by `s` leaves the mathematical problem and exact Newton correction unchanged.

### Parameter matrix

- `x`: `1e-6`, `1`, `1e6`
- Initial relative error in `z`: `1e-12`, `1e-6`, `1e-2`, `1`, `1e6`
- Row scale `s`: `1e-14` through `1e14` in powers of `1e2`
- `constraint_init_tolerance_`: `1e-4` through `1e-14` in powers of `1e2`
- `constraint_init_max_iterations_`: `0`, `1`, `2`, `10`
- One cell and a heterogeneous multi-cell case

### Record

- Solver state
- Number of residual evaluations and Newton updates
- Corrected `z`
- Algebraic state error `z - x^2`
- Raw scaled residual `G_s`
- Unscaled residual `G_s / s`
- Weighted Newton-correction norm
- Whether the last update solved the equation despite a failure return

### Required result

- Corrected state and success/failure classification are invariant to `s`, apart from an identified floating-point limit.
- A linear problem solved by the final allowed Newton update is reported as converged.
- A small scaled residual cannot falsely accept a large algebraic state error.

## Experiment 2: Fixed-step DAE convergence order

### Manufactured DAE

```text
x' = -z
G(x, z) = z - x^2 = 0
x(0) = 1
z(0) = 1
```

The exact solution is

```text
x(t) = 1 / (1 + t)
z(t) = 1 / (1 + t)^2
```

This is nonlinear, index one, and couples the differential dynamics directly to the algebraic variable.

### Parameter matrix

- Methods:

  - Four-stage differential-algebraic Rosenbrock, expected order 3
  - Six-stage differential-algebraic Rosenbrock, expected order 4

- Final time: `T = 1`
- Fixed steps: `H = 2^-k`, with `k = 3 ... 11`
- Constraint row scales: `1e-12`, `1`, `1e12`
- State scales: multiply both exact variables by representative small and large concentration scales

Set `h_start_`, `h_min_`, and `h_max_` to `H`, use tolerances loose enough to avoid adaptive rejection, and assert that the requested fixed-step sequence was actually taken.

### Record

- Global error in `x(T)`
- Global error in `z(T)`
- Final and maximum sampled `|G|`
- Accepted and rejected steps
- Function calls, Jacobian evaluations, decompositions, and solves
- Observed orders

  ```text
  p_obs = log2(error(H) / error(H/2))
  ```

### Required result

- Establish the observed differential, algebraic, and constraint-residual orders for both methods.
- Errors and orders remain invariant to constraint-row scaling.
- No unexplained error plateau occurs before the expected floating-point or conditioning limit.

## Experiment 3: Initialization tolerance versus integration accuracy

Use the manufactured DAE from Experiment 2, but perturb the initial algebraic value away from the manifold.

### Parameter matrix

- `rtol`: `1e-3`, `1e-5`, `1e-7`, `1e-9`
- Component `atol`: chosen consistently for each state scale
- `constraint_init_tolerance_`: `1e-4` through `1e-14`
- Initial algebraic relative error: `1e-12` through `1e2`
- Candidate dimensionless Newton thresholds `eta`: `1`, `0.1`, `0.01`

### Compare convergence tests

1. Current raw residual:

   ```text
   max |G| < tolerance
   ```

2. Explicitly scaled residual supplied by the constraint.

3. Weighted algebraic Newton correction:

   ```text
   max_a |delta_y_a| /
     (atol_a + rtol * max(|y_a|, |y_a + delta_y_a|)) < eta
   ```

4. Weighted correction plus a required residual decrease.

### Record

- Projection error before integration
- Global solution error at `T`
- Constraint residual at `T`
- Initialization work and failure rate
- Point at which tightening initialization convergence no longer improves integration accuracy

### Required result

- Initialization error remains safely below the requested integration error.
- The selected criterion is invariant to row scaling and concentration units.
- The criterion does not over-solve initialization by many orders of magnitude.

## Experiment 4: Nonlinear basin and algebraic-block conditioning

### Scalar nonlinear constraint

Use

```text
G(x, z) = z^2 - x = 0
```

with initial guesses for `z` that are positive, negative, zero, near the root, and many orders of magnitude from it.

### Coupled algebraic constraints

Use two algebraic variables with a tunable nearly singular block, for example

```text
G1 = z1 + z2 - b1
G2 = z1 + (1 + epsilon) z2 - b2
```

with `epsilon = 1e-2 ... 1e-14`.

### Variants

- Current full Newton update
- Backtracking or line-search Newton prototype
- Row-equilibrated residual/Jacobian prototype
- One cell and heterogeneous multi-cell cases

### Record

- Residual and correction history
- Algebraic-block condition estimate
- Convergence rate and stagnation
- Backtracks
- NaN, Inf, singular solve, or maximum-iteration termination
- Whether a small residual corresponds to a small algebraic state error
- State contents returned after failure

### Required result

- Clearly define the supported initial-guess and conditioning regime.
- Failures are deterministic, diagnosable, and do not masquerade as convergence.
- Determine whether line-search globalization is necessary for supported chemistry constraints.

## Experiment 5: LTE dilution, output cadence, and clamping

### Algebraic-fraction sweep

Take one differential equation with known solution and add `m` harmless slaved algebraic variables:

```text
m = 0, 1, 9, 99
Gj = zj - phi_j(x) = 0
```

Measure differential error and adaptive step count while holding all differential tolerances fixed.

### Multi-cell sweep

Repeat with `1`, `10`, `100`, and `1000` cells, including one difficult outlier cell among otherwise easy cells.

### Output-cadence sweep

Integrate the same total interval using:

```text
1, 10, 100, and 1000 public Solve() calls
```

Every case must use the same requested tolerances and total final time.

### Clamp-activation case

Construct a constrained decay or conservation problem in which a differential variable reaches a small negative numerical value and is clamped after `Solve()`.

### Record

- Differential and algebraic solution errors
- LTE norm and accepted/rejected steps
- Constraint residual before and after each public `Solve()` boundary
- Constraint residual after clamping
- Initialization iterations and total projection cost
- Maximum cellwise normalized error, not only the global RMS

### Required result

- Adding zero-LTE algebraic rows does not silently loosen differential accuracy, or the dilution is quantified and explicitly accepted as policy.
- Splitting an interval into output calls changes results only within the requested tolerance and has a documented work cost.
- Post-solve clamping does not leave an unacceptable off-manifold state.

## Experiment 6: Chemistry validation

Run these only after the manufactured experiments identify a candidate convergence rule.

### Robertson QSSA DAE

- Existing stiffness sweep: `k2 = 3e5 ... 3e9`
- `rtol = 1e-4 ... 1e-10`
- Initialization-tolerance sweep around the candidate criterion
- Constraint row scaling: at least `1e-8`, `1`, `1e8`
- Output cadence: single interval, current 19-point grid, and a much denser grid
- Tight full-ODE reference for `A` and `C`

Record differential error, algebraic `B` error relative to the QSSA root, physical and scaled residuals, initialization work, solver counters, and feasibility.

### Tropospheric O3-NOx-HOx QSSA DAE

- Existing photolysis or solar-factor sweep
- Raw production residuals and currently scaled residuals for `O1D`, `O`, `OH`, and `HO2`
- Reference-concentration scaling varied independently of the mathematical constraint
- Initialization and integration tolerance sweeps
- Output-cadence sweep
- Tight full-ODE reference for long-lived species and a direct algebraic-root reference for radicals

Record per-radical relative algebraic error, physical production residual, normalized residual, differential-species accuracy, initialization work, step counts, and failures.

### Chemistry acceptance criteria

- Row scaling does not change the accepted physical solution or convergence classification.
- Tightening algebraic state `atol` alone does not recreate the old step-count inflation.
- Differential accuracy remains consistent with requested tolerances.
- Constraint initialization is neither a dominant cost nor a source of output-cadence dependence.
- Existing conservation, non-negativity, overshoot, and step-economy regressions pass.

## Instrumentation and output schema

Write generated data beneath the build tree, for example `build/dae_experiments/`, so experiment output is not tracked accidentally.

### Case-level CSV

Include at least:

```text
experiment,case_id,method,state_scale,row_scale,rtol,atol,
constraint_tolerance,max_init_iterations,output_calls,cells,
solver_state,init_iterations,steps,accepted,rejected,function_calls,
jacobian_updates,decompositions,solves,error_differential,
error_algebraic,max_raw_residual,max_normalized_residual,wall_us
```

### Initialization-iteration CSV

Include at least:

```text
case_id,cell,constraint_row,iteration,raw_residual,
normalized_residual,correction,weighted_correction,residual_ratio,
accepted_update,backtracks
```

### Result summary

For every experiment, produce:

- A machine-readable CSV
- A concise Markdown table
- Log-log convergence plots where applicable
- A written interpretation that distinguishes observation from inference

## Decision criteria

Select a convergence design only if it satisfies all of the following:

1. Invariant to algebraically equivalent row scaling over the usable floating-point range.
2. Tied to the requested state accuracy rather than arbitrary residual units.
3. Correctly reports convergence after the final allowed successful update.
4. Robust for the supported nonlinear and coupled constraints.
5. Does not over-solve initialization relative to integration accuracy.
6. Preserves the expected DAE method convergence orders.
7. Does not reintroduce algebraic-`atol` step inflation.
8. Keeps feasibility and physical-bound handling separate from LTE accuracy control.
9. Produces stable results under changes in output cadence, cell count, and algebraic-variable count.
10. Passes all existing DAE, constraint, analytical, overshoot, and chemistry regressions.

## Expected design decision

The leading candidate is a weighted Newton-correction convergence test tied to the state `rtol` and component `atol`, combined with finite-value checks, residual-decrease monitoring, and line-search globalization when a full Newton step is not productive.

The experiments must determine whether this is sufficient or whether constraints also need an explicit per-row residual scale. Physical inequality enforcement should remain a separate acceptance or feasibility mechanism rather than being encoded through an algebraic raw-step-change error.

## Proposed execution order

1. Restore the baseline build and diagnostics.
2. Run row-scaling invariance and maximum-iteration boundary tests.
3. Establish fixed-step convergence orders.
4. Determine the relationship between initialization and integration tolerances.
5. Characterize nonlinear and ill-conditioned convergence.
6. Test RMS dilution, output cadence, and clamping.
7. Prototype the best convergence criterion behind the experiment harness.
8. Repeat manufactured experiments with the prototype.
9. Run Robertson and tropospheric chemistry validation.
10. Convert the accepted numerical behavior into focused regression tests before changing the production default.
