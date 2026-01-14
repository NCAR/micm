# Claude Code Project Context - MICM DAE Constraint System

## Project Overview

MICM (Model Independent Chemistry Module) is a C++ chemistry solver library. We've implemented a DAE (Differential-Algebraic Equation) constraint system to allow algebraic constraints to be solved alongside the existing ODE kinetics solver.

## Current State

### What's Complete
- **Constraint infrastructure**: Base `Constraint` class, `EquilibriumConstraint`, and `ConstraintSet` manager
- **SolverBuilder API**: `SetConstraints()` method to inject constraints into the solver
- **State integration**: Constraint variables included in state with proper identity diagonal (1s for ODE, 0s for algebraic)
- **Rosenbrock integration**: Constraint forcing and Jacobian terms added to solver loop
- **AlphaMinusJacobian fix**: Modified to add alpha to ALL diagonals for DAE support
- **Test coverage**: 55 tests passing, including 8 constraint unit tests and 6 integration tests

### DAE Usage Notes
- Use smaller time steps (dt ~0.001) for stiff constraint coupling
- Apply projection step after each solve to enforce constraints exactly
- See `DAESolveWithConstraint` test for working example

## Key Files

### Constraint Implementation
- `include/micm/constraint/constraint.hpp` - Base class
- `include/micm/constraint/equilibrium_constraint.hpp` - Equilibrium constraint
- `include/micm/constraint/constraint_set.hpp` - Collection manager

### Solver Integration
- `include/micm/solver/solver_builder.hpp` - `SetConstraints()` declaration
- `include/micm/solver/solver_builder.inl` - Build logic with constraints
- `include/micm/solver/rosenbrock.inl` - Constraint terms in solve loop, `AlphaMinusJacobian()`
- `include/micm/solver/rosenbrock.hpp` - Solver class declaration
- `include/micm/solver/rosenbrock_temporary_variables.hpp` - Matrix sizing

### Tests
- `test/unit/constraint/test_constraint_set.cpp` - Unit tests
- `test/integration/test_equilibrium.cpp` - Integration tests

### Documentation
- `ARCHITECTURE.md` - Implementation details
- `TODO.md` - Current status and remaining work
- `TESTS.md` - Complex test case specifications

---

## Resume Prompt

Copy and paste this to continue working on the DAE constraint system:

```
I'm continuing work on the MICM DAE constraint system. Please read:
- TODO.md - for current priorities and next steps
- TESTS.md - for test case specifications

The next task is implementing the Chapman mechanism with QSSA constraint.
Build and run tests first to verify the current state, then proceed with the TODO.
```

---

## Build Commands

```bash
cd /Users/fillmore/EarthSystem/MICM/build
cmake --build . -j$(sysctl -n hw.ncpu)
ctest --output-on-failure
```

## Available Skills

The following custom skills are available for this project:

- `/micm-test` - Build and run tests (all, specific pattern, or build-only)
- `/debug-test` - Debug a failing test with tracing support
- `/solver-debug` - Investigate numerical issues in solvers

---

## Technical Context

### Sign Convention
- `SubtractJacobianTerms()` stores `-J_true` in the jacobian matrix
- `AlphaMinusJacobian()` adds alpha to ALL diagonals: `jacobian[i][i] += alpha`
- Result: `alpha - J[i][i]` for all rows (both ODE and algebraic)

### Matrix Structure
```
For N species + M constraints:

Jacobian (N+M) x (N+M):
[  J_kinetics    |      0         ]  <- N rows (ODE)
[  dG/d[species] | dG/d[constraint]]  <- M rows (algebraic)

After AlphaMinusJacobian: all diagonals have alpha added
```

### Test Commands
```bash
# Quick constraint tests
ctest --output-on-failure -R "constraint|equilibrium"

# All tests
ctest --output-on-failure

# Single test with verbose output
ctest --output-on-failure -R "test_name" -V
```

---

## Numerical Debugging Strategies

### When solver produces NaN/Inf or exploding values:
1. Add tracing to `rosenbrock.inl` Solve() loop to print:
   - Matrix diagonal values after AlphaMinusJacobian
   - K values before/after linear solve
   - Forcing vector values
2. Check matrix conditioning: compare magnitude of diagonal elements
3. Verify sign conventions: SubtractJacobianTerms stores -J, not J

### Key diagnostic points in Rosenbrock:
- `initial_forcing` after AddForcingTerms - should be O(rate × concentration)
- `jacobian` diagonals after AlphaMinusJacobian - should all be similar magnitude
- `K[stage]` after linear solve - should scale with H

### Common DAE issues:
- Constraint rows not regularized → ill-conditioned matrix
- c/H terms amplify unscaled K values → explosion
- Fix: ensure all diagonals get alpha added

### Debug tracing template for rosenbrock.inl:
```cpp
// Add after AlphaMinusJacobian call:
std::cerr << "Diagonals after AlphaMinusJacobian:" << std::endl;
for (std::size_t i = 0; i < state.jacobian_.NumColumns(); ++i)
  std::cerr << "  [" << i << "]: " << state.jacobian_.DiagonalElement(0, i) << std::endl;

// Add after linear solve:
std::cerr << "K[" << stage << "] after solve:" << std::endl;
for (std::size_t i = 0; i < K[stage].NumColumns(); ++i)
  std::cerr << "  [" << i << "]: " << K[stage][0][i] << std::endl;
```

---

## Session History

The plan file at `~/.claude/plans/abundant-snacking-puddle.md` contains the original implementation plan if needed for reference.
