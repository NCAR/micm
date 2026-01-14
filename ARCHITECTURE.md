# MICM DAE Constraint System Architecture

## Overview

This document describes the implementation of the Differential-Algebraic Equation (DAE) constraint system added to MICM. The system allows algebraic constraints to be solved alongside ordinary differential equations (ODEs) for chemical kinetics.

## Component Architecture

### 1. Constraint Classes (`include/micm/constraint/`)

#### Base Class: `Constraint`
- **File**: `constraint.hpp`
- **Purpose**: Abstract base class for all constraint types
- **Key Members**:
  - `name_`: Constraint identifier
  - `species_dependencies_`: List of species this constraint depends on
- **Virtual Methods**:
  - `Residual()`: Evaluate constraint residual G(x)
  - `Jacobian()`: Compute partial derivatives dG/dx

#### Equilibrium Constraint: `EquilibriumConstraint`
- **File**: `equilibrium_constraint.hpp`
- **Purpose**: Enforce chemical equilibrium relationships
- **Formula**: `G = K_eq * prod([reactants]^stoich) - prod([products]^stoich) = 0`
- **Key Members**:
  - `reactants_`: Vector of (species_name, stoichiometry) pairs
  - `products_`: Vector of (species_name, stoichiometry) pairs
  - `equilibrium_constant_`: K_eq value
- **Jacobian Calculation**:
  - For reactant R with stoichiometry n: `dG/d[R] = K_eq * n * [R]^(n-1) * prod(others)`
  - For product P with stoichiometry m: `dG/d[P] = -m * [P]^(m-1) * prod(others)`

#### Constraint Set: `ConstraintSet`
- **File**: `constraint_set.hpp`
- **Purpose**: Manage collection of constraints, following `ProcessSet` pattern
- **Key Members**:
  - `constraints_`: Vector of unique_ptr to Constraint objects
  - `constraint_row_offset_`: Starting row index for constraints (= number of species)
  - `dependency_ids_`: Flattened vector of species indices for all constraints
  - `jacobian_flat_ids_`: Sparse matrix indices for Jacobian entries
- **Key Methods**:
  - `AddForcingTerms()`: Add constraint residuals to forcing vector
  - `SubtractJacobianTerms()`: Subtract constraint Jacobian from sparse matrix
  - `NonZeroJacobianElements()`: Return set of (row, col) pairs for sparsity pattern
  - `SetJacobianFlatIds()`: Map constraint Jacobian to sparse matrix flat indices

### 2. Solver Builder Integration (`include/micm/solver/solver_builder.hpp/inl`)

#### New Members Added to `SolverBuilder`:
```cpp
std::shared_ptr<std::vector<std::unique_ptr<Constraint>>> constraints_;
std::size_t constraint_count_ = 0;
std::vector<std::string> constraint_names_{};
```

#### New Method: `SetConstraints()`
```cpp
SolverBuilder& SetConstraints(std::vector<std::unique_ptr<Constraint>>&& constraints);
```
- Takes ownership of constraint vector via move semantics
- Stores constraints in shared_ptr for builder copyability
- Sets `constraint_count_` automatically

#### Build() Modifications:
1. **Extended Variable Map**: Creates map including constraint variables:
   ```cpp
   extended_variable_map[names[i]] = number_of_species + i;
   ```
   This allows constraints to reference their own variables (e.g., "constraint_0")

2. **Constraint Deep Copy**: Clones constraints for builder reusability:
   ```cpp
   if (auto* eq = dynamic_cast<const EquilibriumConstraint*>(c.get()))
     constraints_copy.push_back(std::make_unique<EquilibriumConstraint>(*eq));
   ```

3. **Jacobian Merging**: Combines ODE and constraint Jacobian sparsity patterns:
   ```cpp
   auto constraint_jac_elements = constraint_set.NonZeroJacobianElements();
   nonzero_elements.insert(constraint_jac_elements.begin(), constraint_jac_elements.end());
   ```

4. **State Parameters**: Includes constraint variables in state:
   ```cpp
   StateParameters state_parameters = {
     .number_of_species_ = number_of_species,
     .number_of_constraints_ = number_of_constraints,
     // ...
   };
   ```

### 3. State Modifications (`include/micm/solver/state.hpp/inl`)

#### Identity Diagonal for DAE:
```cpp
for (std::size_t i = 0; i < state_size_; i++)
  upper_left_identity_diagonal_.push_back(1.0);  // ODE variables
for (std::size_t i = 0; i < constraint_size_; i++)
  upper_left_identity_diagonal_.push_back(0.0);  // Algebraic variables
```

#### Variable Names:
State includes both species and constraint variable names in `variable_names_` and `variable_map_`.

### 4. Rosenbrock Solver Integration (`include/micm/solver/rosenbrock.inl`)

#### Forcing Terms:
```cpp
initial_forcing.Fill(0);
rates_.AddForcingTerms(state.rate_constants_, Y, initial_forcing);
if (constraints_.Size() > 0)
{
  constraints_.AddForcingTerms(Y, initial_forcing);
}
```

#### Jacobian Terms:
```cpp
state.jacobian_.Fill(0);
rates_.SubtractJacobianTerms(state.rate_constants_, Y, state.jacobian_);
if (constraints_.Size() > 0)
{
  constraints_.SubtractJacobianTerms(Y, state.jacobian_);
}
```

#### Matrix Formation:
```cpp
AlphaMinusJacobian(
    state.jacobian_,
    state.upper_left_identity_diagonal_,  // M: 1 for ODE, 0 for algebraic
    singular,
    H * parameters.gamma_[0]);
```

### 5. Temporary Variables (`include/micm/solver/rosenbrock_temporary_variables.hpp`)

Modified to include constraint dimensions:
```cpp
Ynew_(number_of_grid_cells, state_parameters.number_of_species_ + state_parameters.number_of_constraints_),
initial_forcing_(number_of_grid_cells, state_parameters.number_of_species_ + state_parameters.number_of_constraints_),
// ...
```

## Data Flow

```
User Code                    SolverBuilder                    Solver
    |                             |                              |
    |--SetConstraints()---------->|                              |
    |                             |--Build()                     |
    |                             |   |--Create extended_variable_map
    |                             |   |--Deep copy constraints
    |                             |   |--Create ConstraintSet
    |                             |   |--Merge Jacobian elements
    |                             |   |--Create StateParameters
    |                             |   |--Return Solver----------->|
    |                             |                              |
    |--solver.Solve()------------------------------------------------>|
    |                             |                              |--rates_.AddForcingTerms()
    |                             |                              |--constraints_.AddForcingTerms()
    |                             |                              |--rates_.SubtractJacobianTerms()
    |                             |                              |--constraints_.SubtractJacobianTerms()
    |                             |                              |--AlphaMinusJacobian()
    |                             |                              |--LinearSolver.Solve()
    |<--result---------------------------------------------------------|
```

## Matrix Structure

For a system with N species and M constraints, the augmented system is:

```
State vector: [x_1, x_2, ..., x_N, c_1, c_2, ..., c_M]
              |---- species ----|  |-- constraints --|

Jacobian (N+M) x (N+M):
              species cols     constraint cols
            [  J_kinetics    |      0         ]  species rows
            [----------------|----------------]
            [  dG/d[species] | dG/d[constraint]]  constraint rows

Identity diagonal (mass matrix M):
            [1, 1, ..., 1, 0, 0, ..., 0]
             |-- N 1s --|  |-- M 0s --|
```

## Sign Conventions

Following the `ProcessSet` convention:
- `SubtractJacobianTerms()` subtracts the true Jacobian: `jacobian -= J_true`
- After subtraction: `jacobian = -J_true`
- `AlphaMinusJacobian()` computes: `M[i][i] - alpha * jacobian[i][i]`

## Known Limitations

1. **Off-diagonal scaling**: The current `AlphaMinusJacobian()` only modifies diagonal elements. For proper DAE support, all elements should be scaled by alpha (= h*gamma).

2. **Diagonal sign**: For constraint rows with M[i][i]=0, the diagonal computation produces `-alpha` instead of `+alpha`, which can cause numerical instability.

3. **Stiff coupling**: Large equilibrium constants (K_eq >> 1) create stiff algebraic-differential coupling that the current solver struggles with.

## File Listing

```
include/micm/constraint/
  constraint.hpp              - Base Constraint class
  constraint_set.hpp          - ConstraintSet manager
  equilibrium_constraint.hpp  - EquilibriumConstraint implementation

include/micm/solver/
  solver_builder.hpp          - SetConstraints() declaration
  solver_builder.inl          - SetConstraints() and Build() implementation
  rosenbrock.inl              - Constraint integration in solve loop
  rosenbrock_temporary_variables.hpp - Matrix sizing with constraints
  state.hpp/inl               - Identity diagonal and state structure

test/unit/constraint/
  test_constraint_set.cpp     - Unit tests for ConstraintSet

test/integration/
  test_equilibrium.cpp        - API integration tests
```
