# MICM DAE Constraint System - TODO

## Current Status

### Completed Features

- [x] **Constraint Base Class** (`include/micm/constraint/constraint.hpp`)
  - Abstract base with `Residual()` and `Jacobian()` methods
  - Species dependencies tracking

- [x] **EquilibriumConstraint** (`include/micm/constraint/equilibrium_constraint.hpp`)
  - Implements G = K_eq * [reactants] - [products] = 0
  - Proper Jacobian calculation with stoichiometry handling
  - Edge case handling for zero concentrations

- [x] **ConstraintSet** (`include/micm/constraint/constraint_set.hpp`)
  - Follows ProcessSet pattern
  - `AddForcingTerms()` - adds constraint residuals to forcing vector
  - `SubtractJacobianTerms()` - subtracts constraint Jacobian
  - `NonZeroJacobianElements()` - returns sparsity pattern
  - `SetJacobianFlatIds()` - maps to sparse matrix indices

- [x] **SolverBuilder SetConstraints API** (`include/micm/solver/solver_builder.hpp/inl`)
  - `SetConstraints(std::vector<std::unique_ptr<Constraint>>&&)` method
  - Extended variable map including constraint variables
  - Deep copying for builder reusability
  - Jacobian sparsity pattern merging

- [x] **State Infrastructure**
  - Identity diagonal with 1s for ODE, 0s for algebraic variables
  - Constraint variables in state and variable_map
  - Proper matrix sizing including constraints

- [x] **Rosenbrock Integration**
  - Constraint forcing added after ODE forcing
  - Constraint Jacobian subtracted after ODE Jacobian
  - Temporary variables sized for species + constraints

### Test Coverage

#### Unit Tests (test/unit/constraint/test_constraint_set.cpp)
- [x] `ConstraintSetBasicAPI` - Basic construction and setup
- [x] `ConstraintSetForcingTerms` - Residual calculation
- [x] `ConstraintSetJacobianTerms` - Jacobian structure
- [x] `ConstraintSetJacobianValues` - Jacobian numerical values
- [x] `MultipleConstraints` - Multiple constraints in one set
- [x] `ThreeDStateOneConstraint` - 3D state with 1 constraint
- [x] `FourDStateTwoConstraints` - 4D state with 2 constraints
- [x] `CoupledConstraintsSharedSpecies` - Constraints sharing species

#### Integration Tests (test/integration/test_equilibrium.cpp)
- [x] `ReversibleReactionReachesEquilibrium` - ODE equilibrium test
- [x] `SimpleIsomerization` - ODE isomerization test
- [x] `ConstraintSetAPITest` - Direct ConstraintSet API test
- [x] `SetConstraintsAPIWorks` - SolverBuilder with 1 constraint
- [x] `SetConstraintsAPIMultipleConstraints` - SolverBuilder with 2 constraints

**All 55 tests passing**

---

## Remaining Work

### High Priority: Solver Modifications for DAE Support

The current Rosenbrock solver needs modifications to properly solve DAE systems. The `AlphaMinusJacobian()` function in `rosenbrock.hpp` has issues:

#### Issue 1: Off-diagonal Scaling

**Current behavior**: Only diagonal elements are modified
```cpp
// Current implementation
for (std::size_t i = 0; i < jacobian.NumColumns(); ++i)
{
  double val = upper_left_identity_diagonal[i] - alpha * jacobian.DiagonalElement(cell, i);
  jacobian.DiagonalElement(cell, i) = val;
}
```

**Required behavior**: All elements should be scaled by alpha for proper (M - hγJ) formation
```cpp
// Required for DAE support
// Diagonal: M[i][i] - alpha * J[i][i]
// Off-diagonal: 0 - alpha * J[i][j] = -alpha * J[i][j]
```

**Impact**: For ODE systems (M=I everywhere), the diagonal dominates stability. For DAE systems with M[i][i]=0 for algebraic equations, the off-diagonal scaling becomes critical.

#### Issue 2: Sign Convention for Constraint Rows

**Problem**: After `SubtractJacobianTerms`, jacobian stores `-J_true`. For constraint rows where M[i][i]=0:
- Current: `0 - alpha * jacobian[i][i]` = `0 - alpha * (-J[i][i])` = `+alpha * J[i][i]`
- For G with dG/d[C] = -1: result is `-alpha` (wrong sign)
- Needed: `+alpha` for proper matrix conditioning

#### Proposed Fix

Modify `AlphaMinusJacobian()` to:
1. Scale ALL Jacobian elements by alpha, not just diagonal
2. Correctly handle sign for constraint rows

**Warning**: Changes must preserve behavior for existing ODE tests.

### Medium Priority: Additional Constraint Types

- [ ] **ConservationConstraint**: Mass conservation (sum of species = constant)
- [ ] **BoundConstraint**: Keep species within bounds
- [ ] **CustomConstraint**: User-defined constraint functions

### Low Priority: Enhancements

- [ ] Constraint-aware error estimation
- [ ] Specialized DAE solver parameters (RODAS variants)
- [ ] Index-2 DAE support (if needed)
- [ ] GPU constraint support (CUDA)

---

## Test Results Summary

```
Test Suite                              Tests    Status
---------------------------------------- ------- -------
constraint_set (unit)                    8       PASS
equilibrium (integration)                5       PASS
All other existing tests                 42      PASS
---------------------------------------- ------- -------
TOTAL                                    55      PASS
```

Last test run: All 55 tests passing

---

## Notes

### Why Full DAE Solving Isn't Working Yet

The constraint infrastructure is complete, but actual DAE time integration fails because:

1. **Matrix Formation**: `AlphaMinusJacobian()` doesn't properly form (M - hγJ) for constraint rows
2. **Numerical Stability**: Large K_eq values cause stiff coupling that amplifies numerical errors
3. **Step Size Control**: The existing error estimator doesn't account for algebraic variables

### Workaround

For now, fast equilibria can be handled using:
1. **Reversible reactions**: Model A + B <-> C with forward/backward rate constants
2. **Operator splitting**: Solve equilibria separately after kinetics
3. **QSSA**: Quasi-steady-state approximation for fast intermediates

### References

- Hairer & Wanner, "Solving Ordinary Differential Equations II: Stiff and DAE Problems"
- RODAS: Rosenbrock methods for DAEs
- CAM-Chem equilibrium handling
