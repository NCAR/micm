# Code Review: MICM Rosenbrock DAE Solver

**Reviewer:** Claude Opus 4.6
**Date:** 2026-02-11
**Branch:** `dae-constraint-enforcement`
**Scope:** Full review of DAE constraint system integrated into the Rosenbrock solver

---

## Executive Summary

The DAE (Differential-Algebraic Equation) constraint system is architecturally sound. The "replace-state-rows" design correctly replaces algebraic species rows in-place rather than appending extra rows, keeping matrix dimensions consistent throughout the solver. The mass-matrix diagonal (`upper_left_identity_diagonal_`) properly distinguishes ODE rows (M=1) from algebraic rows (M=0), and is correctly applied in `AlphaMinusJacobian` and the c/H stage combination terms. The `ProcessSet` correctly skips kinetic contributions for algebraic variable rows.

Two passes of review have been completed. Most bugs and code quality issues from the first pass are now fixed. Remaining open issues are documented below alongside findings from the Codex review (`REVIEW_CODEX.md`).

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

**File:** `include/micm/constraint/constraint.hpp`

`Constraint::GetAlgebraicSpecies()` returns `std::string` by value while `EquilibriumConstraint::AlgebraicSpecies()` returns `const std::string&`. The `std::visit` lambda also returns by value. Not a bug, but creates an unnecessary copy.

---

## 2. Robustness Issues

### 2.1 Unguarded `std::pow()` with potentially negative concentrations - FIXED

**File:** `include/micm/constraint/equilibrium_constraint.hpp`

Added `std::max(0.0, ...)` guard before `std::pow()` calls in `Residual()` to prevent NaN from transient negative concentrations.

### 2.2 MEDIUM: Equilibrium constant is compile-time static - OPEN (Enhancement)

**File:** `include/micm/constraint/equilibrium_constraint.hpp`

Real atmospheric chemistry equilibrium constants are temperature-dependent (Arrhenius/van't Hoff). The current implementation stores a single constant. This limits the system to isothermal problems.

**Recommendation:** Add support for a rate-constant-like callable that takes temperature as input, similar to `ArrheniusRateConstant`.

### 2.3 MEDIUM: `FourStageDifferentialAlgebraicRosenbrockParameters` needs citation - OPEN

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

## 4. Findings from Codex Review (REVIEW_CODEX.md)

These findings were identified by the Codex reviewer and confirmed as still valid.

### 4.1 HIGH: `State` copy operations slice solver-specific temporary buffers - OPEN

**Files:** `include/micm/solver/state.hpp:110-111, 139-140`

```cpp
temporary_variables_ =
    other.temporary_variables_ ? std::make_unique<TemporaryVariables>(*other.temporary_variables_) : nullptr;
```

Both copy constructor and copy assignment recreate `temporary_variables_` as a base `TemporaryVariables` object. The Rosenbrock solver `static_cast`s to `RosenbrockTemporaryVariables` in `rosenbrock.inl:19`. Copying a state and then solving the copy triggers undefined behavior.

**Proposed fix:** Add `virtual std::unique_ptr<TemporaryVariables> Clone() const = 0;` to the base class and override in derived types. Use `Clone()` in State copy operations instead of `std::make_unique<TemporaryVariables>(...)`.

**Needed tests:**
- Copy a Rosenbrock state and solve both original and copy
- Equivalent test for BackwardEuler state copy

### 4.2 MEDIUM: Row-replacement Jacobian can overwrite when dependency column repeats - OPEN

**Files:** `include/micm/constraint/constraint_set.hpp:354-377`

Replace mode writes with `=` per dependency entry. If a species appears on both sides of an equilibrium (e.g., `A + B <-> A + C`), the dependency list has duplicate columns and only the last partial derivative for that `(row, col)` is retained.

Currently safe for typical equilibrium constraints (reactants and products are disjoint), but a latent bug for constraints with overlapping species.

**Proposed fix:** In replacement mode, clear target row entries once, then accumulate (`+=`) all dependency contributions. Or deduplicate dependencies and pre-aggregate derivatives before writeback.

**Needed tests:**
- Constraint with overlapping species (e.g., `A + B <-> A + C`) verifying Jacobian `dG/dA` includes both partial derivative terms

### 4.3 MEDIUM: Unused-species validation ignores constraint-only species - OPEN

**File:** `include/micm/solver/solver_builder.inl:174`

`UnusedSpeciesCheck()` only considers `RatesPolicy::SpeciesUsed(reactions_)`. A species that appears only in a constraint (not in any reaction) would be flagged as unused when `SetIgnoreUnusedSpecies(false)`.

**Proposed fix:** Union reaction-used species with constraint dependency/algebraic target species before the difference check.

**Needed tests:**
- Builder test with strict unused-species checking + constraint-only species should build successfully

### 4.4 LOW: Jacobian sparsity includes kinetic algebraic-row entries - OPEN (Enhancement)

**File:** `include/micm/solver/solver_builder.inl:341-354`

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

### Tests Still Needed

| Priority | Test | Reason |
|----------|------|--------|
| HIGH | State copy + solve (Rosenbrock) | Validates fix for 4.1 (slicing bug) |
| HIGH | State copy + solve (BackwardEuler) | Same slicing issue |
| HIGH | Multi-grid-cell DAE solve | Verifies vectorized paths |
| MEDIUM | Overlapping-species constraint Jacobian | Validates fix for 4.2 |
| MEDIUM | Constraint-only species + strict unused check | Validates fix for 4.3 |
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

### Architecture Risks

1. **State copy slicing** (4.1): Most critical risk — any code path that copies a state and solves the copy will crash.
2. **NormalizedError weighting**: Treats algebraic variables identically to ODE variables. May need different scaling for stiffer systems.
3. **No index detection**: Assumes all constraints are index-1. Higher-index DAEs would require index reduction (not implemented; should be documented as a limitation).

---

## 7. Summary of Open Items

| Priority | Issue | Type | Section |
|----------|-------|------|---------|
| HIGH | State copy slicing of temporary variables | Bug | 4.1 |
| MEDIUM | Jacobian overwrite with duplicate dependency columns | Bug | 4.2 |
| MEDIUM | Unused-species check ignores constraint species | Bug | 4.3 |
| MEDIUM | Temperature-dependent equilibrium constant | Enhancement | 2.2 |
| MEDIUM | 4-stage DAE parameters citation | Documentation | 2.3 |
| LOW | Remove dead append-rows code paths | Simplification | 3.5 |
| LOW | AlgebraicSpecies() return type inconsistency | Code quality | 1.4 |
| LOW | Filter kinetic sparsity from algebraic rows | Enhancement | 4.4 |
