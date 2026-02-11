# Rosenbrock DAE Solver Code Review (Codex) - Pass 6

## Scope
Re-reviewed `dae-constraint-enforcement` after applying the high-priority fix for algebraic-row sparsity pruning in `ProcessSet::SetJacobianFlatIds`.

Primary files in this pass:
- `include/micm/process/process_set.hpp`
- `test/unit/solver/test_solver_builder.cpp`
- `include/micm/solver/solver_builder.inl` (interaction point)

## Validation
- Targeted unit test: `ctest --output-on-failure -R solver_builder` -> passed.
- Full suite: `ctest --output-on-failure` -> 56/56 passed.
- Disabled tests scan: `rg -n "DISABLED_" test` -> none.

## Findings (ordered by severity)
1. High (Fixed): algebraic-row sparsity pruning could throw `Zero element access` in solver build.
   Root cause:
   - Builder pruned kinetic Jacobian rows for algebraic variables.
   - `ProcessSet::SetJacobianFlatIds` still attempted `VectorIndex(...)` lookups for those pruned entries.
   Fix:
   - `SetJacobianFlatIds` is now algebraic-row aware and skips sparse lookup for algebraic rows while preserving loop alignment with placeholder flat IDs.
   References:
   - `include/micm/process/process_set.hpp:245`
   - `include/micm/process/process_set.hpp:252`
   - `include/micm/process/process_set.hpp:259`
   - `include/micm/process/process_set.hpp:269`

2. High (Fixed): missing regression coverage for the exact build-time throw topology.
   Added test:
   - `SolverBuilder.DAEKineticAlgebraicRowPruningDoesNotThrow` reproduces the prior failure pattern (`A -> C` kinetics plus `K_eq*B - C = 0` constraint) and verifies build succeeds.
   Reference:
   - `test/unit/solver/test_solver_builder.cpp:159`

## Remaining Non-Blocking Items
- Medium (deferred): equilibrium constants are static and not temperature-dependent.
- Medium (deferred): citation/doc provenance for the 4-stage DAE Rosenbrock parameter set.

## Overall Assessment
No open high-severity defects remain from this review cycle. The previously identified build-time regression is fixed and covered by a targeted unit test, and the full test suite is passing.
