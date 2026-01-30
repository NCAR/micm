# PR Plan: DAE Constraint System

This document outlines the strategy for merging the DAE constraint system into main via 3 incremental, reviewable PRs.

## PR 1: Constraint Base Classes

**Branch**: `dae-1-constraint-classes`
**Target**: `main`
**Estimated Size**: ~1,300 lines (all new files)

### Description
Introduces the constraint class hierarchy with full unit test coverage. Pure addition with zero changes to existing solver code.

### Files
```
include/micm/constraint/constraint.hpp          # Base class
include/micm/constraint/equilibrium_constraint.hpp  # Equilibrium implementation
include/micm/constraint/constraint_set.hpp      # Collection manager
include/micm/Constraint.hpp                     # Umbrella header

test/unit/constraint/CMakeLists.txt
test/unit/constraint/test_constraint.cpp
test/unit/constraint/test_constraint_set.cpp
test/unit/CMakeLists.txt                        # Add add_subdirectory(constraint)
```

### Review Focus
- API design for Constraint base class
- EquilibriumConstraint formulation (K_eq = [products]/[reactants])
- ConstraintSet management of indices and Jacobian elements
- Unit test coverage

### Acceptance Criteria
- [ ] All existing tests pass
- [ ] New constraint unit tests pass
- [ ] No changes to existing solver behavior

---

## PR 2: State/Builder Infrastructure

**Branch**: `dae-2-state-infrastructure`
**Target**: `main` (after dae-1-constraint-classes merged)
**Estimated Size**: ~200 lines (modifications to existing files)

### Description
Adds state and builder infrastructure to support constraints without changing solver behavior. When no constraints are configured, the solver works exactly as before.

### Files
```
include/micm/solver/state.hpp                   # Add constraint_size_, upper_left_identity_diagonal_
include/micm/solver/state.inl                   # Sizing logic for constraints
include/micm/solver/solver_builder.hpp          # SetConstraintCount(), SetConstraintNames()
include/micm/solver/solver_builder.inl          # Implement above (store values only)
include/micm/solver/rosenbrock_temporary_variables.hpp  # Sizing updates if needed
```

### Review Focus
- State struct additions for tracking constraint count
- Builder API design (fluent interface consistency)
- Backward compatibility (default constraint_count_=0)

### Acceptance Criteria
- [ ] All existing tests pass unchanged
- [ ] New API methods compile and store values
- [ ] No solver behavior changes when constraints not used

---

## PR 3: Rosenbrock DAE Integration

**Branch**: `dae-3-rosenbrock-integration`
**Target**: `main` (after dae-2-state-infrastructure merged)
**Estimated Size**: ~900 lines

### Description
Wires constraints into the Rosenbrock solver loop and adds integration tests demonstrating DAE solving. Also includes all documentation.

### Files

**Solver Integration**:
```
include/micm/solver/solver_builder.hpp          # SetConstraints() declaration
include/micm/solver/solver_builder.inl          # Full Build() with ConstraintSet creation
include/micm/solver/rosenbrock.hpp              # ConstraintSet member, constructor update
include/micm/solver/rosenbrock.inl              # AddForcingTerms(), SubtractJacobianTerms() in loop
```

**Integration Tests**:
```
test/integration/test_equilibrium.cpp
test/integration/CMakeLists.txt
```

**Documentation**:
```
ARCHITECTURE.md
TODO.md
TESTS.md
CLAUDE.md
PLAN.md
PR-PLAN.md
.claude/commands/debug-test.md
.claude/commands/micm-test.md
.claude/commands/solver-debug.md
```

**Minor Test Updates** (if any):
```
test/unit/solver/test_linear_solver_in_place_policy.hpp
test/unit/solver/test_linear_solver_policy.hpp
test/unit/solver/test_lu_decomposition_in_place_policy.hpp
test/unit/solver/test_lu_decomposition_policy.hpp
```

### Review Focus
- DAE formulation in Rosenbrock solve loop
- AlphaMinusJacobian behavior (adds alpha to ALL diagonals)
- ConstraintSet integration in Build()
- Integration test correctness (equilibrium actually achieved)
- Sign conventions documented in ARCHITECTURE.md

### Acceptance Criteria
- [ ] All existing tests pass
- [ ] Integration tests demonstrate correct equilibrium behavior
- [ ] Documentation accurately describes implementation

---

## Merge Order

```
main
  ↑
  dae-1-constraint-classes
  ↑
  dae-2-state-infrastructure
  ↑
  dae-3-rosenbrock-integration
  ↑
rosenbrock_dae (current branch)
```

## Creating the Branches

```bash
# From rosenbrock_dae branch:

# PR 1: Cherry-pick or checkout only constraint files
git checkout main
git checkout -b dae-1-constraint-classes
# Add constraint files...

# PR 2: Based on PR 1
git checkout dae-1-constraint-classes
git checkout -b dae-2-state-infrastructure
# Add state/builder infrastructure...

# PR 3: Based on PR 2 (or just rebase rosenbrock_dae)
git checkout dae-2-state-infrastructure
git checkout -b dae-3-rosenbrock-integration
# Add remaining files...
```

## Notes

- Each PR should be independently reviewable and testable
- PRs can be reviewed in parallel but must merge sequentially
- If reviews identify issues, fixes go into the appropriate PR in the chain
