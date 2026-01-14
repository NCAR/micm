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

- [x] **AlphaMinusJacobian DAE Fix** (`include/micm/solver/rosenbrock.inl`)
  - Modified to add alpha to ALL diagonal elements (not just ODE rows)
  - Treats algebraic constraints as stiff ODEs for numerical stability
  - Prevents K values from exploding due to c/H terms in stage computation
  - Both vectorized and non-vectorized versions updated

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
- [x] `DAESolveWithConstraint` - Full DAE integration with projection step

**All 55 tests passing**

---

## Remaining Work

### ~~High Priority: Solver Modifications for DAE Support~~ ✅ COMPLETED

The `AlphaMinusJacobian()` function has been fixed. See "How the Fix Works" section below.

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
equilibrium (integration)                6       PASS
All other existing tests                 41      PASS
---------------------------------------- ------- -------
TOTAL                                    55      PASS
```

Last test run: All 55 tests passing

---

## How the Fix Works

### The Problem

The original `AlphaMinusJacobian()` function only added alpha (= 1/(h*gamma)) to ODE rows based on `upper_left_identity_diagonal_`. For constraint rows where M[i][i]=0, the diagonal wasn't regularized. This caused:

1. **Massive ill-conditioning**: ODE diagonals were ~10^10 while constraint diagonals were ~1
2. **K value explosion**: The c/H terms in Rosenbrock stage computation amplified K values for constraints
3. **Numerical instability**: Constraint variables would explode to values like 10^16

### The Solution

Modified `AlphaMinusJacobian()` to add alpha to ALL diagonal elements:

```cpp
// Add alpha to ALL diagonals for proper DAE support.
// For ODE variables (M[i][i]=1): forms αM - J = α - J
// For algebraic variables (M[i][i]=0): also adds α to regularize the constraint.
// This treats algebraic constraints as stiff ODEs (ε*z' = g(y,z) with ε=hγ),
// which ensures K values scale with H and prevents numerical instability
// from the c/H terms in Rosenbrock stage computation.
for (const auto& i_elem : state.jacobian_diagonal_elements_)
  jacobian_vector[i_elem] += alpha;
```

This treats algebraic constraints as stiff ODEs with ε = hγ, ensuring K values scale properly with the step size H.

### Using DAE Constraints

For best results when using algebraic constraints:

1. **Use smaller time steps**: Start with dt = 0.001 or smaller
2. **Apply projection step**: After each solve, enforce the constraint exactly:
   ```cpp
   // For constraint C = K_eq * B
   state.variables_[0][C_idx] = K_eq * state.variables_[0][B_idx];
   ```
3. **Monitor convergence**: Check `result.state_ == SolverState::Converged`

### Alternative Approaches

For some use cases, these alternatives may work better:

1. **Reversible reactions**: Model A + B <-> C with forward/backward rate constants
2. **Operator splitting**: Solve equilibria separately after kinetics
3. **QSSA**: Quasi-steady-state approximation for fast intermediates

### References

- Hairer & Wanner, "Solving Ordinary Differential Equations II: Stiff and DAE Problems"
- RODAS: Rosenbrock methods for DAEs
- CAM-Chem equilibrium handling
