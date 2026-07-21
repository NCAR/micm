# Schur Reduction of Algebraic Rows — Ceiling Measurement and Design

- **Date:** 2026-07-20
- **Branch:** `dae-schur-reduction`
- **Prototype:** `benchmark/schur_ceiling.cpp`
- **Status:** ceiling measured (go); core implementation designed, not yet built.

## Why

MICM keeps algebraic variables in the factored linear system: the equilibrium
family solves 3N rows with two-thirds algebraic and loses 1.5–2.9× wall-clock
to the full ODE despite taking 26% fewer steps. Index-1 guarantees the
algebraic block \(J_{zz}\) is nonsingular, so the stage solves can be reduced
to the differential block.

## Measured ceiling

`schur_ceiling` compares the full 3N ODE, the 3N DAE, and the exactly reduced
N-dimensional problem (equilibrium + conservation eliminated analytically —
the best any Schur implementation could approach) at rtol 1e-6, medians of 7
interleaved repetitions:

| N | ODE µs | DAE µs | reduced µs | DAE/ODE | ceiling reduced/ODE |
|---:|---:|---:|---:|---:|---:|
| 1 | 129 | 196 | 50 | 1.52 | 0.39 |
| 4 | 204 | 485 | 70 | 2.38 | 0.35 |
| 16 | 653 | 1899 | 360 | 2.91 | 0.55 |
| 64 | 2146 | 5346 | 418 | 2.49 | 0.19 |
| 256 | 18560 | 49220 | 3684 | 2.65 | 0.20 |

The ceiling is ~5× faster than the full ODE and ~13× faster than the current
DAE at N = 256. Even capturing half the ceiling clears the phase-4 acceptance
criterion (DAE/ODE < 1 for N ≥ 16) with a wide margin. **Recommendation: go.**

## Core design

Partition the stage matrix by `state.upper_left_identity_diagonal_`:

```
[alpha*M - J] = [ alpha*I - J_xx    -J_xz ]
                [      -J_zx        -J_zz ]
```

Per step (J fixed, alpha per stage-factorization):

1. Factor `J_zz` (n_z × n_z per block; for QSSA radical families it is small
   and often nearly diagonal; for the equilibrium family it is 2×2
   block-diagonal per triple).
2. Form the Schur complement `S = (alpha*I - J_xx) - J_xz * J_zz^{-1} * J_zx`
   (n_x × n_x). Sparsity: `pattern(J_xx) ∪ pattern(J_xz·J_zz^{-1}·J_zx)`;
   with block-diagonal `J_zz` the fill is the product pattern of the
   constraint couplings only — for the equilibrium family S is diagonal.
3. Factor `S`; each stage solve is then: eliminate the z part of the RHS,
   solve S for the x part, back-substitute z.

Implementation shape:

- A new linear-solver policy (`SchurLinearSolver`) selected by a builder
  option, so the default path is untouched. Symbolic phase at build time
  (patterns of `J_zz` factors and `S`), numeric phase in `LinearFactor`.
- The constraint-initialization solve (`-dG/dy` with identity differential
  rows) reduces trivially — the differential block is the identity — and can
  reuse the same machinery.
- The alpha shift only touches the differential diagonal, so `J_zz`'s factors
  are reusable across stage alphas within a step when the method re-factors
  (the non-inplace path's alpha-delta trick still applies to `S`).

## Risks

- Fill-in when constraints couple many differential species through `J_zz`
  (dense radical families); the symbolic phase must bound it and fall back to
  the unreduced path when `nnz(S)` exceeds a threshold.
- Vectorized (grouped-cell) orderings need the same treatment per group.
- CUDA parity: out of scope for the first CPU implementation; the policy
  switch keeps GPU builds on the existing path.

## Follow-on

Implement `SchurLinearSolver` (CPU standard ordering first), validate with the
fixed-step order harness (orders must be unchanged — the reduction is exact
linear algebra), and re-run `equilibrium_efficiency`, `robertson_dae`,
`tropospheric_dae`, and `schur_ceiling` to measure captured fraction of the
ceiling.
