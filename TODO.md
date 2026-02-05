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

### Next Up: Chapman Mechanism with QSSA Constraint

Implement the Chapman stratospheric ozone mechanism with atomic oxygen as a QSSA algebraic constraint. This is a real-world atmospheric chemistry test case that exercises the core DAE capability.

**Reactions**:
```
O2 + hv -> 2O           j1 (slow)
O + O2 + M -> O3        k2 (fast)
O3 + hv -> O2 + O       j3 (fast)
O + O3 -> 2O2           k4 (slow)
```

**QSSA Constraint**: `d[O]/dt ≈ 0`
```
2*j1*[O2] + j3*[O3] - k2*[O][O2][M] - k4*[O][O3] = 0
```

**Requirements**:
- [ ] Implement `QSSAConstraint` or `CustomConstraint` class
- [ ] Chapman mechanism test with O as algebraic variable
- [ ] Validate against analytical solution and pure ODE reference
- [ ] Diurnal cycle test with time-varying photolysis rates

**CSV Output for Plotting**:
Numerical integration tests should write CSV data to stdout (one line per timestep):
```
time,O2,O,O3,O_analytical,O3_analytical
0.000000,4.0e+17,0.0e+00,5.0e+12,1.2e+07,5.0e+12
0.001000,4.0e+17,1.2e+07,5.0e+12,1.2e+07,5.0e+12
...
```
- First line: header with column names
- Subsequent lines: time, all state variables, analytical values (if available)
- Use scientific notation for consistency
- Redirect to file when running: `./test_chapman > chapman_results.csv`
- Enables plotting with Python/gnuplot for visual validation

See `TESTS.md` for full specification and additional test cases.

### Medium Priority: Additional Constraint Types

- [ ] **ConservationConstraint**: Mass conservation (sum of species = constant)
- [ ] **QSSAConstraint**: Quasi-steady-state approximation (d[X]/dt = 0)
- [ ] **CustomConstraint**: User-defined constraint functions via lambdas
- [ ] **BoundConstraint**: Keep species within bounds

### Low Priority: Enhancements

- [ ] Constraint-aware error estimation
- [ ] Automatic projection (built into solver)
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
