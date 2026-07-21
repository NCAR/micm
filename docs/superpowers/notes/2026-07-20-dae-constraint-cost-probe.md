# DAE Constraint Cost Probe — Phase 3 Profiling

- **Date:** 2026-07-20
- **Branch:** `dae-constraint-cost`
- **Probe:** `benchmark/constraint_cost_probe.cpp` (Robertson, rtol 1e-6,
  k2 = 3e7, medians of 21 interleaved repetitions)

## Measurements

Dense-output run (19 log-spaced `Solve()` segments to t = 1e6):

| Variant | µs | accepted | decompositions | µs/step |
|---|---:|---:|---:|---:|
| full ODE | 137.6 | 430 | 430 | 0.320 |
| full ODE + null external model | 138.2 | 430 | 430 | 0.321 |
| QSSA-DAE | 156.8 | 417 | 436 | 0.376 |

Single-call run (one `Solve(1e6)`, initialization amortized to once):

| Variant | µs | accepted | µs/step |
|---|---:|---:|---:|
| full ODE | 99.0 | 308 | 0.321 |
| QSSA-DAE | 109.4 | 292 | 0.375 |

## Findings

1. **External-model dispatch is negligible.** Routing the ProcessSet sweeps
   through a null external model's `std::function` callbacks costs
   0.6–0.9 µs of a 138 µs run (< 0.7%). The plan's "devirtualize the
   callback layer" task would buy essentially nothing and is dropped.
2. **Constraint initialization is effectively free.** The dense-output run
   performs 19 extra factorizations (one single-iteration Newton per public
   `Solve()` call, `init_iters = 19`), yet the implied per-init cost from the
   dense-vs-single-call comparison is −0.19 µs/call — zero within
   measurement noise. No "skip redundant reinitialization" optimization is
   warranted on cost grounds.
3. **The overhead is per-step machinery: ~53 ns/step (+16.6%) on this
   3-species system**, identical in the dense-output and single-call runs.
   Candidates, in likely order of importance on small systems:
   - the mass-coupling path in the stage loop (six invocations per 4-stage
     step, each through a `DenseMatrixPolicy::Function` closure with
     column-view/ForEachRow machinery — heavy relative to a 3×1 column);
   - the constraint residual (3×/step) and Jacobian (1×/step) sweeps;
   - algebraic-row handling in `AlphaMinusJacobian` and the error norm.

## Revised phase-3 plan

- Replace the `mass_coupling` closure with a direct loop over a precomputed
  list of differential-variable indices (no per-call closure dispatch, no
  column-view construction), then re-run this probe — expected to recover
  the bulk of the 53 ns/step.
- Re-profile on the 16-species tropospheric system (work-precision data
  shows ~10% overhead at rtol 1e-6 there) — the balance between coupling
  machinery and constraint math shifts with matrix size.
- Keep the ±5% wall-clock acceptance criterion from the plan; measure with
  the phase-2 interleaved work-precision rig.

## Reproduction

```bash
cmake -S . -B build -D MICM_ENABLE_BENCHMARKS=ON
cmake --build build --target constraint_cost_probe -j
./build/constraint_cost_probe
```
