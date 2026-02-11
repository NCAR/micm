# Rosenbrock DAE Solver Code Review (Codex) - Pass 2

## Scope
Re-reviewed the full Rosenbrock DAE path after the mass-matrix row-replacement changes:

- solver assembly and stage math
- constraint mapping/evaluation/Jacobian
- builder/state wiring
- unit/regression/integration coverage

Primary files reviewed:
- `include/micm/solver/rosenbrock.inl`
- `include/micm/constraint/constraint_set.hpp`
- `include/micm/constraint/equilibrium_constraint.hpp`
- `include/micm/process/process_set.hpp`
- `include/micm/solver/solver_builder.inl`
- `include/micm/solver/state.hpp`
- `include/micm/solver/state.inl`
- `test/integration/test_equilibrium.cpp`
- `test/unit/constraint/test_constraint_set.cpp`
- `test/unit/solver/test_rosenbrock.cpp`
- `test/regression/RosenbrockChapman/regression_test_solve.cpp`

## Validation
- Full suite: `ctest --output-on-failure` in `build` -> 56/56 passed.
- Disabled test scan: `rg -n "DISABLED_" test` -> none.

## Findings (ordered by severity)
1. High (Open): `State` copy operations still slice solver-specific temporary buffers.
   Evidence:
   - `State` copy ctor/assignment recreate `temporary_variables_` as `TemporaryVariables` base object.
   - Rosenbrock solve path unconditionally `static_cast`s back to `RosenbrockTemporaryVariables`.
   References:
   - `include/micm/solver/state.hpp:110`
   - `include/micm/solver/state.hpp:139`
   - `include/micm/solver/rosenbrock.inl:18`
   Risk:
   - Copying a state returned by `GetState(...)` and solving the copy can trigger undefined behavior/crash.
   Proposed fix:
   - Add `virtual std::unique_ptr<TemporaryVariables> Clone() const = 0;` and clone derived buffers in copy operations.
   - Alternative: make `State` non-copyable and force move semantics.
   Needed tests:
   - New test that copies a Rosenbrock state and solves both original and copy.
   - Equivalent test for Backward Euler state copy behavior.

2. Medium (Open): Row-replacement Jacobian path can overwrite contributions when a dependency column repeats.
   Evidence:
   - Replacement mode writes with `=` per dependency entry.
   - Equilibrium dependencies are built by concatenating reactants/products and are not deduplicated.
   References:
   - `include/micm/constraint/constraint_set.hpp:351`
   - `include/micm/constraint/constraint_set.hpp:366`
   - `include/micm/constraint/equilibrium_constraint.hpp:97`
   Risk:
   - If one species appears multiple times in dependencies (e.g., species on both sides), only the last partial derivative is retained for that `(row, col)`.
   Proposed fix:
   - In replacement mode, clear target row entries once, then accumulate (`+=`) all dependency contributions.
   - Or deduplicate dependencies and pre-aggregate derivatives before writeback.
   Needed tests:
   - Constraint with overlapping species across sides (e.g., `A + B <-> A + C`) verifying Jacobian `dG/dA` includes both terms.
   - Run in standard and vectorized sparse policies.

3. Medium (Open/UX): strict unused-species validation does not account for constraint-only species.
   Evidence:
   - Unused-species check is based only on reaction usage, not constraints.
   References:
   - `include/micm/solver/solver_builder.inl:174`
   - `include/micm/solver/solver_builder.inl:338`
   Risk:
   - With `.SetIgnoreUnusedSpecies(false)`, species used only by algebraic constraints are flagged as unused.
   Proposed fix:
   - Union reaction-used species with constraint dependencies/algebraic targets before check.
   Needed tests:
   - Builder test with strict unused-species checking + constraint-only species should build successfully.

4. Low (Enhancement): Jacobian sparsity still includes kinetic algebraic-row structure even when those rows are masked.
   Evidence:
   - Kinetic nonzeros are collected before algebraic row IDs are applied.
   References:
   - `include/micm/solver/solver_builder.inl:341`
   - `include/micm/solver/solver_builder.inl:354`
   Impact:
   - Extra fill pattern/factorization work on rows that are later overwritten by constraints.
   Proposed enhancement:
   - Filter kinetic nonzero entries whose row is algebraic before building Jacobian structure.

## Confirmed Improvements From Prior Pass
- Full mass-matrix row-replacement enforcement is active and integration-tested (`DAESolveWithConstraint` enabled).
- `alpha * M - J` and stage `c/H` coupling are now mass-matrix scaled.
- Constraint/vectorized layout correctness and normalized-error indexing issues remain fixed.
- Six-stage DA regression test now uses six-stage DA parameters.

## Additional Needed Tests (beyond current suite)
1. `State` copy + solve regression (Rosenbrock and Backward Euler).
2. Duplicate-dependency constraint Jacobian correctness in replacement mode.
3. Strict unused-species + constraints integration case.
4. Jacobian-structure regression checking algebraic-row pruning if enhancement #4 is implemented.
