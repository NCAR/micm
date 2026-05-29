# DAE Rosenbrock Benchmark — Algebraic Constraints vs. Full ODE on a Stiff System

**Date:** 2026-05-28
**Branch:** `dae-rosenbrock-benchmark`
**Status:** Design — awaiting review

## Purpose

Produce publication- and documentation-quality evidence that replacing a stiff
fast mode with an algebraic constraint lets the Rosenbrock DAE solver take far
fewer, larger steps than a full ODE solve of the same system — at matched
accuracy after the initial transient.

The deliverable is reproducible numbers (CSV) and figures, not a CI guard. A
thin GoogleTest wrapper that asserts the headline inequality may follow later
but is out of scope here.

## Scope

**Phase 1 (this spec): QSSA reduction of the Robertson problem.** The canonical
stiff ODE benchmark, reduced by quasi-steady-state approximation on the fast
intermediate species B.

**Phase 2 (noted, separate spec later): reversible-equilibrium case.** A system
whose fast process is a genuine reversible reaction, where the algebraic
`EquilibriumConstraint` is *exact* in the stiff limit (not an approximation).
This is deliberately deferred; the structure below is built so the sweep/metrics/
plotting harness is reusable for it.

## Background and honest framing

The Robertson system (matching `test/tutorial/configs/robertson`):

```
r1:  A      -> B            rate = k1 * [A]            (k1 = 0.04)
r2:  2B     -> B + C        rate = k2 * [B]^2          (k2 = 3e7,  the stiff knob)
r3:  B + C  -> A + C        rate = k3 * [B] * [C]      (k3 = 1e4)
```

Forcing:

```
dA/dt = -k1*A           + k3*B*C
dB/dt =  k1*A - k2*B^2  - k3*B*C
dC/dt =         k2*B^2
```

Conservation: `A + B + C = 1` (with A(0)=1, B(0)=C(0)=0).

Stiffness arises because B is a fast transient: it rises sharply, then sits in
quasi-steady-state at a tiny value (~1e-5). A full ODE solver must resolve the
fast B mode with small steps; its step count grows with k2.

**The DAE reduction is a QSSA on B**, replacing B's ODE row with `dB/dt = 0`:

```
G(A,B,C) = k1*A - k2*B^2 - k3*B*C = 0        (quadratic in B)
```

This is an *approximation*: it discards B's initial fast transient. Therefore
the comparison is **not** "identical answer, faster." It is:

> The QSSA-reduced DAE matches the full ODE for A and C to within a small
> tolerance *after the transient*, while costing a fraction of the solver steps,
> and that cost gap widens as the system gets stiffer.

That accuracy-vs-cost tradeoff is the publishable result and must be reported as
such — both the agreement and the savings, side by side.

## Why a custom constraint (not a built-in type)

`G` is quadratic in B, so neither `EquilibriumConstraint`
(`K_eq*[reactants]^nu - [products]`) nor `LinearConstraint` expresses it. The
external/custom-constraint path — the pattern in
`test/integration/stub_aerosol_with_constraints.hpp` and exercised by
`test/integration/test_external_model_constraints.cpp` — supplies an arbitrary
residual and Jacobian, so **no new core constraint types are required.**

Constraint Jacobian entries (state convention; solver subtracts dG/dy):

```
dG/dA = k1
dG/dB = -2*k2*B - k3*C
dG/dC = -k3*B
```

Non-zero Jacobian elements for B's row: (B,A), (B,B), (B,C).

## Architecture

A standalone benchmark executable. Generation of numbers is fully separated from
any correctness assertion. Five components, each independently understandable
and testable:

1. **`RobertsonReactions(k1, k2, k3)`** — single source of truth for the species
   and the three reactions, used by *both* solver builds so they are provably the
   same chemistry. Returns the gas phase + `std::vector<Process>`. k2 is settable
   per case for the stiffness sweep (rates are user-defined parameters in the
   existing config, so they are set on the state, not baked into the solver).

2. **`RobertsonQssaConstraint`** — custom external-model constraint object
   following the `stub_aerosol_with_constraints.hpp` interface. Declares B as the
   algebraic variable, depends on {A, B, C}, and provides the residual `G` and
   Jacobian above. Holds k1, k2, k3 (k2 set per case).

3. **`RunCase(method, k2, rtol, atol)`** — builds the requested solver (full ODE
   or DAE), runs Robertson from the standard initial condition over a fixed set
   of logarithmically spaced output times, and returns:
   - aggregated `SolverStats` summed across all internal `Solve` calls
     (`number_of_steps_`, `accepted_`, `rejected_`, `function_calls_`,
     `jacobian_updates_`),
   - the sampled A and C trajectories,
   - optional median wall-clock over N repetitions (informational only).

4. **`ReferenceSolution(k2)`** — a full-ODE solve at very tight tolerances
   (e.g. rtol=1e-10), giving the trajectory the DAE accuracy is measured against.

5. **`main`** — sweeps k2 over orders of magnitude (3e5 … 3e9), runs both methods
   at matched tolerances, computes post-transient max relative error of the DAE A
   and C trajectories against the reference, and writes one CSV.

### Identical-configuration discipline

Both solvers use the same Rosenbrock order
(`FourStageDifferentialAlgebraicRosenbrockParameters`), the same rtol/atol, the
same output time grid, and the same initial condition. The *only* differences
are: (a) B is marked algebraic via the constraint in the DAE build, and (b) the
DAE replaces B's ODE row with `G=0`. This is what makes the step-count comparison
fair. r2/r3 remain in the process set for the DAE build because they feed the A
and C forcing; only B's row is replaced.

## Data flow

```
main
 ├─ for k2 in sweep:
 │    ├─ ref      = ReferenceSolution(k2)                      // tight-tol full ODE
 │    ├─ ode      = RunCase(FULL_ODE, k2, rtol, atol)
 │    ├─ dae      = RunCase(DAE_QSSA, k2, rtol, atol)
 │    ├─ err_dae  = max post-transient rel. error of dae vs ref (A, C)
 │    └─ write CSV rows: {method, k2, rtol, steps, accepted, rejected,
 │                        function_calls, jacobian_updates, wallclock_median_us,
 │                        max_rel_err}
 └─ done -> robertson_dae_benchmark.csv
```

`plot_robertson_dae.py` reads the CSV and emits two figures:
- **Fig 1 — cost vs stiffness:** `number_of_steps` (and `jacobian_updates`) vs k2
  for both methods on log–log axes. Expected: ODE rising, DAE ~flat.
- **Fig 2 — accuracy vs stiffness:** DAE max post-transient relative error vs k2,
  showing the reduction stays accurate as the savings grow.

## Metrics

**Primary (deterministic, hardware-independent):** the `SolverStats` counters.
These reproduce exactly on any machine and are the core evidence.

**Secondary (informational):** median wall-clock over N repetitions, reported but
explicitly labeled hardware-dependent. Never the sole basis of a claim.

**Accuracy:** max relative error of DAE A and C against the tight-tolerance
reference, sampled at output times after a transient cutoff `t_skip` (chosen past
B's initial rise, e.g. t > 1e-3 s). The transient cutoff is a reported parameter.

## Error handling

- Any `Solve` returning a state other than `Converged` is recorded in the CSV
  (the case is not silently dropped) and flagged; the sweep continues. A full ODE
  solve hitting `max_number_of_steps_` at extreme stiffness is itself a data
  point, not a failure to hide.
- No silent truncation of the sweep: every requested k2 produces a CSV row, even
  if it reports non-convergence.

## Build and layout

- `benchmark/robertson_dae.cpp` — the executable.
- `benchmark/plot_robertson_dae.py` — committed plotting script (matplotlib).
- `benchmark/CMakeLists.txt` — target behind a new `MICM_ENABLE_BENCHMARKS`
  option (OFF by default), linking the existing MICM target. Wired into the top
  `CMakeLists.txt`.
- Reuses `test/tutorial/configs/robertson` reaction/species forms (constructed in
  code rather than parsed, to keep k2 a free parameter).

## Testing

The benchmark is a measurement tool, not a unit under test, but it must be
*correct* to be trustworthy:

- **Constraint Jacobian verification:** use the existing
  `micm/util/jacobian_verification.hpp` to finite-difference-check
  `RobertsonQssaConstraint`'s analytic Jacobian (the same guard
  `test_external_model_constraints.cpp` relies on).
- **Sanity assertions inside the benchmark** (not gtest): post-transient DAE vs
  reference error below a loose ceiling at the nominal k2; conservation
  `A+B+C ≈ 1` for the full ODE; constraint residual `|G| ≈ 0` for the DAE.
- A dedicated GoogleTest guard of the headline inequality is **out of scope**
  (candidate Phase 1.5 follow-on).

## Out of scope

- GPU/CUDA path (CPU only for the paper numbers).
- Multi-grid-cell / vectorized comparisons.
- The reversible-equilibrium case (Phase 2, separate spec) — but the sweep,
  metrics, CSV schema, and plotting are designed to be reused by it.

## Open parameters to finalize during implementation

- Exact k2 sweep range and point count.
- rtol/atol values for the matched comparison and for the reference.
- Transient cutoff `t_skip` and the output time grid.
- Wall-clock repetition count N.
