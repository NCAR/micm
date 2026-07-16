# DAE Rosenbrock Constraint Tolerance and Convergence Experiments

- **Date:** 2026-07-15
- **Branch:** `dae-rosenbrock-benchmark`
- **Baseline commit:** `f11b5a7b`
- **Status:** Weighted-correction initialization implemented and validated; two controller/conditioning policy questions remain explicit.

## Implementation

- Added `benchmark/dae_convergence_experiments.cpp` and the CMake target `dae_convergence_experiments`.
- Generated experiment output is written below the build tree and is not tracked.
- Replaced the raw residual stopping rule

  ```text
  max |G_i(y)| < tolerance
  ```

  with the maximum weighted algebraic Newton correction

  ```text
  max_a |delta_y_a| /
    (atol_a + rtol * max(|y_a|, |y_a + delta_y_a|)) <= eta
  ```

- The default is `eta = 0.1`.
- Added affine-covariant backtracking based on a simplified Newton correction computed with the current factored Jacobian. Complete constraint-row scaling cancels from this merit function.
- A final permitted Newton update is evaluated and may return convergence.
- Failed initialization restores the complete caller-provided state.
- If nonnegative clamping changes a differential variable after a successful or usable partial DAE solve, algebraic variables are reprojected. Successful reprojection preserves the original solver status.

## Baseline versus candidate

| Observation | Baseline | Candidate |
|---|---:|---:|
| Projection cases | 5,400 | 5,400 |
| False convergence classifications | 564 | 0 |
| Exact solution reached but failure returned | 842 | 0 |
| Nonlinear square-root failures | 9 | 3 |
| Failures from finite near-zero guesses | 6 | 0 |
| Failures at exactly singular `z=0` | 3 | 3 |
| Clamp-activation residual | `5.0e-1` | `0` |

The three retained nonlinear failures are deterministic `InfDetected` results at an exactly zero algebraic derivative, one for each row scale. Failure rollback leaves the input state unchanged.

## Row scaling and tolerance coupling

- Projection success/failure classification had zero mismatches across row scales `1e-14 ... 1e14`.
- The nonlinear integration sweep had zero classification mismatches across row scales `1e-8`, `1`, and `1e8`.
- Correction thresholds `eta = 1 ... 1e-6` produced identical integration step counts and final errors in all 288 nonlinear tolerance-coupling cases.
- Mean initialization iterations increased only from `4.25` at `eta=1` to `5.19` at `eta=1e-6`; tighter initialization did not improve integration accuracy.
- `eta=0.1` is retained as a conservative default without material over-solving.

## Fixed-step convergence

The quadratic manifold in the original plan is solved to roundoff by the four-stage method, so the order experiment also uses

```text
x' = -z
z = S sin(x / S)
x(t) = 2 S atan(tan(1/2) exp(-t)).
```

At `row_scale=1`, `state_scale=1`, and `H=2^-8`:

| Method | Expected | Observed `x` order | Observed `z` order |
|---|---:|---:|---:|
| Four-stage DAE | 3 | 3.003 | 3.001 |
| Six-stage DAE | 4 | 3.884 | 3.943 |

The six-stage errors reach the floating-point floor at finer steps. Across row scales `1e-12`, `1`, and `1e12`, the worst pre-floor differential-error ratio was `1.00005`.

## Nonlinear basin and conditioning

- Backtracking converges from `z=+/-1e-6` for `z^2-1=0`; undamped Newton exhausted ten iterations from those guesses.
- Initial `z=0` is outside the supported Newton regime because the algebraic Jacobian is exactly singular.
- For the coupled linear system at `epsilon=1e-14`, the computed residual is exactly zero in working precision while the algebraic forward error is `2.22e-2`.
- Residual or correction convergence alone cannot certify forward error once the algebraic block is numerically singular. A condition estimate or pivot-quality diagnostic is the appropriate future safeguard.

## Controller, cadence, and clamping

- Adding `0`, `1`, `9`, or `99` algebraic copies did not change accepted steps or differential error for the slaved-variable test.
- Splitting the same interval into `1`, `10`, `100`, or `1000` public `Solve()` calls changed total accepted steps from `9` to `90`, `900`, and `9000`. Accuracy remained at roundoff. This is the documented cost of restarting `h_start` on every public call.
- With one active cell and all other cells exactly stationary, the current global WRMS controller produced:

  | Cells | Accepted steps | Active-cell relative error |
  |---:|---:|---:|
  | 1 | 44 | `3.73e-7` |
  | 10 | 32 | `1.20e-6` |
  | 100 | 24 | `3.81e-6` |
  | 1000 | 18 | `1.22e-5` |

  This is the expected square-root dilution of a global RMS norm, but it makes per-cell accuracy depend on batch composition. The current implementation is preserved pending a separate decision between global WRMS and maximum cellwise WRMS.

- The clamp-activation experiment previously returned `X=0`, `Z=-0.5`. Conditional post-clamp projection now returns `X=Z=0`, also covers partial solves with accepted steps, and adds no work when clamping is inactive.

## Chemistry validation

- Focused DAE/constraint tests pass: initialization, overshoot, algebraic-error insensitivity, algebraic step economy, and Robertson QSSA.
- Robertson QSSA converged across `k2=3e5 ... 3e9`; accepted-step counts remained `480, 461, 417, 350, 241` at `rtol=1e-6`.
- Tropospheric O3-NOx-HOx QSSA converged across photolysis scales `0.1 ... 2.0`; accepted-step counts remained `376 ... 566`.
- The equilibrium-efficiency benchmark converged for `1 ... 256` independent triples.

## Reproduction

```bash
cmake -S . -B build
cmake --build build --target dae_convergence_experiments -j2
./build/dae_convergence_experiments build/dae_experiments
ctest --test-dir build --output-on-failure \
  -R 'constraint_initialization|dae_constraint_overshoot|dae_algebraic_error_(insensitivity|step_economy)|robertson_qssa'
```

Primary output files:

- `projection_scaling.csv`
- `fixed_step_order.csv`
- `tolerance_coupling.csv`
- `nonlinear_initialization.csv`
- `controller_cadence.csv`
- `summary.md`

## Remaining decisions

1. Preserve global WRMS across cells or adopt maximum cellwise WRMS for batch-invariant per-cell accuracy.
2. Add an algebraic-block condition or pivot-quality diagnostic before treating residual convergence as a forward-error guarantee near singularity.
3. Carry an internal step-size estimate across public `Solve()` boundaries if dense-output cadence is expected to be work-invariant.
