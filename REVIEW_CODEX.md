# Rosenbrock DAE Solver Code Review (Codex)

## Scope
Reviewed the Rosenbrock DAE path across solver core logic, constraint handling, state wiring, and regression/unit/integration tests.

Primary files reviewed:
- `include/micm/solver/rosenbrock.inl`
- `include/micm/constraint/constraint_set.hpp`
- `include/micm/solver/solver_builder.inl`
- `include/micm/solver/state.hpp`
- `test/unit/solver/test_rosenbrock.cpp`
- `test/unit/constraint/test_constraint_set.cpp`
- `test/regression/RosenbrockChapman/regression_test_solve.cpp`
- `test/integration/test_equilibrium.cpp`

Validation run:
- `ctest --output-on-failure -R "^(constraint_set|rosenbrock|equilibrium)$"` (all passed)

## Findings (ordered by severity)
1. Critical (Fixed): Constraint evaluation and Jacobian assembly were incorrect for vectorized layouts.
   Evidence: `ConstraintSet` assumed row-contiguous dense storage and contiguous per-cell sparse blocks, which is invalid for `VectorMatrix` and vector-ordered sparse matrices. This affects both forcing and Jacobian terms in constrained/vectorized runs.
   References: `include/micm/constraint/constraint_set.hpp:204`, `include/micm/constraint/constraint_set.hpp:260`, `include/micm/constraint/constraint_set.hpp:291`.
   Fix applied:
   - Added contiguous row buffering for vectorized dense matrices before calling constraint `Residual`/`Jacobian`.
   - Added vectorized sparse path that writes via `jacobian[i_cell][row][col]` instead of flat-offset arithmetic.

2. High (Fixed): Normalized error used misaligned flattened indexing when constraints add extra columns.
   Evidence: Error normalization iterated flattened `Y` length against flattened `Ynew`/`errors` with larger column counts, causing cross-row mixing in multi-cell constrained states and distorted step control.
   References: fixed logic now at `include/micm/solver/rosenbrock.inl:337` and `include/micm/solver/rosenbrock.inl:374`.
   Fix applied:
   - Rewrote both normalized-error paths to row/column loops over species dimensions (`Y.NumRows() x Y.NumColumns()`), independent of matrix storage layout.

3. Medium (Fixed): The six-stage DAE regression test instantiated the four-stage DAE parameters.
   Evidence: The `SixStageDASolve` test called `FourStageDifferentialAlgebraicRosenbrockParameters`.
   Reference: `test/regression/RosenbrockChapman/regression_test_solve.cpp:42`.
   Fix applied:
   - Updated test to use `SixStageDifferentialAlgebraicRosenbrockParameters`.

4. High (Open): Purely algebraic variable enforcement is still incomplete by design.
   Evidence: Integration test is explicitly disabled and documents that algebraic variables with no kinetic equation are not directly enforced by current row-augmentation approach.
   Reference: `test/integration/test_equilibrium.cpp:235`.
   Recommended enhancement:
   - Replace/augment ODE rows for algebraic variables with constraint rows, or add a projection step after each accepted solve.

5. Medium (Open): `State` copy operations slice `temporary_variables_` and can invalidate Rosenbrock casts.
   Evidence: Copy constructor/assignment copy `TemporaryVariables` base only; Rosenbrock later `static_cast`s to derived temporary variable type.
   References: `include/micm/solver/state.hpp:106`, `include/micm/solver/state.hpp:134`, `include/micm/solver/rosenbrock.inl:18`.
   Recommended enhancement:
   - Add polymorphic cloning (`virtual Clone`) on `TemporaryVariables` and derived implementations, or make `State` non-copyable.

## Fixes Implemented
Code updates completed:
- `include/micm/constraint/constraint_set.hpp`
- `include/micm/solver/rosenbrock.inl`
- `test/regression/RosenbrockChapman/regression_test_solve.cpp`

Test additions completed:
- Added vectorized constraint-layout correctness coverage in `test/unit/constraint/test_constraint_set.cpp:620`.
- Added constrained normalized-error coverage (standard + vectorized) in `test/unit/solver/test_rosenbrock.cpp:61` and `test/unit/solver/test_rosenbrock.cpp:209`.

## Additional Needed Tests
1. Add an enabled integration test for algebraic-variable enforcement once the DAE row/projection strategy is implemented.
   Reference baseline: `test/integration/test_equilibrium.cpp:253`.

2. Add a regression test that copies a `State` and then calls Rosenbrock `Solve`, to catch temporary-variable slicing issues.
   Suggested location: `test/unit/solver/test_state.cpp` or `test/unit/solver/test_rosenbrock.cpp`.

3. Add mixed-policy constrained tests (standard dense + vector sparse, and vector dense + standard sparse) to ensure no hidden layout assumptions remain.
   Suggested locations: `test/unit/constraint/test_constraint_set.cpp` and `test/unit/solver/test_rosenbrock.cpp`.

## Notes
This review prioritized correctness and regression risk in the Rosenbrock DAE path. The implemented fixes address concrete correctness bugs in constrained multi-cell/vectorized operation and improve regression coverage for those paths.
