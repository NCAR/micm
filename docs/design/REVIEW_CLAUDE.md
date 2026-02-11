# Code Review: MICM Rosenbrock DAE Solver

**Reviewer:** Claude Opus 4.6
**Date:** 2026-02-11 (updated after Codex Pass 3)
**Branch:** `dae-constraint-enforcement`
**Scope:** Full review of DAE constraint system integrated into the Rosenbrock solver

---

## Executive Summary

The DAE (Differential-Algebraic Equation) constraint system is architecturally sound. The "replace-state-rows" design correctly replaces algebraic species rows in-place rather than appending extra rows, keeping matrix dimensions consistent throughout the solver. The mass-matrix diagonal (`upper_left_identity_diagonal_`) properly distinguishes ODE rows (M=1) from algebraic rows (M=0), and is correctly applied in `AlphaMinusJacobian` and the c/H stage combination terms. The `ProcessSet` correctly skips kinetic contributions for algebraic variable rows.

Three review passes have been completed (two Claude, three Codex). Most bugs and code quality issues are now fixed. This document tracks the current status of all findings from both reviewers.

---

## 1. Bugs - FIXED

### 1.1 `Solver::Solve()` clamps algebraic variables to zero - FIXED

**File:** `include/micm/solver/solver.hpp`

`Max(0.0)` was clamping every variable to non-negative. Now selectively clamps only ODE variables (where `upper_left_identity_diagonal_[i_var] > 0.0`), leaving algebraic variables unclamped.

### 1.2 `NonZeroJacobianElements()` missing self-diagonal for replace mode - FIXED

**File:** `include/micm/constraint/constraint_set.hpp`

Added explicit insertion of `(constraint_row, constraint_row)` when `replace_state_rows_` is true, ensuring the diagonal exists for AlphaMinusJacobian and LU factorization.

### 1.3 Stray semicolon - FIXED

**File:** `include/micm/solver/rosenbrock.inl` (was line 256)

### 1.4 LOW: `AlgebraicSpecies()` returns a copy vs. const reference inconsistency - OPEN

**File:** `include/micm/constraint/constraint.hpp:77`

`Constraint::GetAlgebraicSpecies()` returns `std::string` by value while `EquilibriumConstraint::AlgebraicSpecies()` returns `const std::string&`. The `std::visit` lambda also returns by value. Not a bug, but creates an unnecessary copy.

---

## 2. Robustness Issues

### 2.1 Unguarded `std::pow()` with potentially negative concentrations in Residual - FIXED

**File:** `include/micm/constraint/equilibrium_constraint.hpp`

Added `std::max(0.0, ...)` guard before `std::pow()` calls in `Residual()` to prevent NaN from transient negative concentrations.

### 2.2 Residual and Jacobian use inconsistent concentration guards - FIXED (Codex Pass 3, Finding 2)

**File:** `include/micm/constraint/equilibrium_constraint.hpp`

Applied `std::max(0.0, ...)` guard consistently in the Jacobian method: aggregate product terms, per-species derivative terms, and `conc == 0` special-case fallbacks now all clamp negative concentrations to zero, matching the Residual policy.

### 2.3 MEDIUM: Equilibrium constant is compile-time static - OPEN (Enhancement)

**File:** `include/micm/constraint/equilibrium_constraint.hpp`

Real atmospheric chemistry equilibrium constants are temperature-dependent (Arrhenius/van't Hoff). The current implementation stores a single constant. This limits the system to isothermal problems.

**Recommendation:** Add support for a rate-constant-like callable that takes temperature as input, similar to `ArrheniusRateConstant`.

### 2.4 MEDIUM: `FourStageDifferentialAlgebraicRosenbrockParameters` needs citation - OPEN

**File:** `include/micm/solver/rosenbrock_solver_parameters.hpp`

The 4-stage DAE parameters have no literature citation. The `SixStageDifferentialAlgebraicRosenbrockParameters` correctly cites Hairer & Wanner (1996).

**Recommendation:** Add a citation or derive order conditions verification for the 4-stage method.

---

## 3. Design Issues

### 3.1 Template parameter typo `Dervied` - FIXED

**File:** `include/micm/solver/rosenbrock.hpp`

### 3.2 Unused constructor parameters - FIXED

**File:** `include/micm/solver/rosenbrock.hpp`, `include/micm/solver/backward_euler.hpp`

Removed `jacobian`, `number_of_species`, and `number_of_constraints` from both Rosenbrock and BackwardEuler constructors.

### 3.3 Stale comments referencing old architecture - FIXED

**File:** `include/micm/solver/rosenbrock.inl`

### 3.4 `std::cout` diagnostic output in tests - FIXED

**File:** `test/integration/test_equilibrium.cpp`

### 3.5 LOW: Dead append-rows code paths - OPEN

**Files:** `constraint_set.hpp`, `state.hpp/inl`

The builder unconditionally uses replace mode when constraints are present. The append-rows `else` branches are untested by integration tests and should eventually be removed.

### 3.6 Include of `<iostream>` in solver header - FIXED

**File:** `include/micm/solver/rosenbrock.hpp`

---

## 4. Findings from Codex Reviews

### 4.1 `State` copy operations slice solver-specific temporary buffers - FIXED

**Files:** `include/micm/solver/state.hpp`, `temporary_variables.hpp`, `rosenbrock_temporary_variables.hpp`, `backward_euler_temporary_variables.hpp`

Added virtual `Clone()` to `TemporaryVariables` base class, overridden in both derived types. State copy ctor/assignment now use `Clone()` to preserve the derived type.

### 4.2 Row-replacement Jacobian overwrite with duplicate dependency columns - FIXED

**File:** `include/micm/constraint/constraint_set.hpp`

Changed replace mode from `=` to `-=` so duplicate dependency columns accumulate correctly. The Jacobian row is already zeroed by `Fill(0)` and ProcessSet skips algebraic rows, so `-=` is safe.

### 4.3 Unused-species validation ignores constraint-only species - FIXED

**File:** `include/micm/solver/solver_builder.inl`

Added constraint dependencies and algebraic targets to the used-species set in `UnusedSpeciesCheck()`.

### 4.4 Overloaded `Solve(time_step, state, params)` bypasses clamping - FIXED (Codex Pass 3, Finding 1)

**File:** `include/micm/solver/solver.hpp`

Extracted clamping into a private `PostSolveClamp()` helper, called from both `Solve` overloads. Both paths now apply identical post-solve clamping (ODE-only for DAE, all variables otherwise).

### 4.5 LOW: No targeted regression tests for recently fixed issues - OPEN (Codex Pass 3, Finding 3)

The fixes for Clone() (4.1), duplicate dependency accumulation (4.2), and constraint-aware unused-species (4.3) lack dedicated regression tests.

**Needed tests:**
- Copy a solver-generated State and solve both original and copy (Rosenbrock + BackwardEuler)
- Constraint with repeated species (e.g., `A + B <-> A + C`) verifying Jacobian accumulation
- `.SetIgnoreUnusedSpecies(false)` with constraint-only species builds successfully

### 4.6 LOW: Jacobian sparsity includes kinetic algebraic-row entries - OPEN (Enhancement)

**File:** `include/micm/solver/solver_builder.inl`

Kinetic nonzero elements are collected before algebraic row IDs are applied, so rows that will be entirely replaced by constraints still carry their kinetic sparsity pattern. This causes extra fill in the LU factorization for no benefit.

**Proposed enhancement:** Filter kinetic nonzero entries whose row is algebraic before building the Jacobian structure.

---

## 5. Tests

### Tests Added (from Review Recommendations)

| Test | Description | Status |
|------|-------------|--------|
| `DAESolveWithFourStageDAEParameters` | 4-stage DAE Rosenbrock parameters | PASSING |
| `DAESolveWithSixStageDAEParameters` | 6-stage RODAS parameters | PASSING |
| `DAESolveWithTwoCoupledConstraints` | Two constraints sharing species B | PASSING |
| `DAEConservationLaw` | Verifies A+B conservation over time | PASSING |
| `DAESolveStiffCoupling` | Large K_eq (1000), no NaN/Inf | PASSING |
| `DAESolveWithNonUnitStoichiometry` | K_eq*[A]^2 - [B] = 0 | PASSING |
| `DAEClampingDoesNotBreakAlgebraicVariables` | Selective clamping verification | PASSING |
| `DAESolveMultiGridCell` | 3 grid cells with different ICs, per-cell checks | PASSING |

### Tests Still Needed

| Priority | Test | Reason |
|----------|------|--------|
| HIGH | State copy + solve (Rosenbrock) | Regression test for 4.1 fix |
| HIGH | State copy + solve (BackwardEuler) | Regression test for 4.1 fix |
| MEDIUM | Overlapping-species constraint Jacobian | Regression test for 4.2 fix |
| MEDIUM | Constraint-only species + strict unused check | Regression test for 4.3 fix |
| LOW | Reorder state with constraints | Markowitz reordering + DAE interaction |
| LOW | Error estimation for algebraic variables | Step size controller behavior |

---

## 6. Architecture Assessment

### What's Done Well

1. **Replace-state-rows design**: Correct DAE formulation with consistent matrix dimensions.
2. **Mass matrix diagonal**: Clean {1.0, 0.0} distinction applied in AlphaMinusJacobian and stage c/H terms.
3. **ProcessSet algebraic guards**: Efficiently prevents kinetic contributions on constraint rows.
4. **DAE-specific Rosenbrock parameters**: FourStage and SixStage (RODAS) parameter sets available.
5. **Validation**: Duplicate constraint row detection, unknown species errors, equilibrium constant validation.
6. **Polymorphic Clone**: State copy now preserves solver-specific temporary variables correctly.

### Architecture Risks

1. **NormalizedError weighting**: Treats algebraic variables identically to ODE variables. May need different scaling for stiffer systems.
2. **No index detection**: Assumes all constraints are index-1. Higher-index DAEs would require index reduction (not implemented; should be documented as a limitation).

---

## 7. Summary of Open Items

| Priority | Issue | Type | Section |
|----------|-------|------|---------|
| MEDIUM | Temperature-dependent equilibrium constant | Enhancement | 2.3 |
| MEDIUM | 4-stage DAE parameters citation | Documentation | 2.4 |
| LOW | No regression tests for Clone/dup-dep/unused-species fixes | Testing | 4.5 |
| LOW | Remove dead append-rows code paths | Simplification | 3.5 |
| LOW | `AlgebraicSpecies()` return type inconsistency | Code quality | 1.4 |
| LOW | Filter kinetic sparsity from algebraic rows | Enhancement | 4.6 |
