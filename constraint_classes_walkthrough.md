# Constraint Classes for DAE Support in MICM

This document provides an educational walkthrough of the new `Constraint`, `EquilibriumConstraint`, and `ConstraintSet` classes introduced in the `dae-1-constraint-classes` branch. These classes enable Differential-Algebraic Equation (DAE) support in MICM's Rosenbrock solver.

## Table of Contents

1. [Overview](#overview)
2. [Class Hierarchy](#class-hierarchy)
3. [The Constraint Base Class](#the-constraint-base-class)
4. [The EquilibriumConstraint Class](#the-equilibriumconstraint-class)
5. [The ConstraintSet Class](#the-constraintset-class)
6. [How ConstraintSet Mirrors ProcessSet](#how-constraintset-mirrors-processset)
7. [Integration with the DAE System](#integration-with-the-dae-system)

---

## Overview

### What are DAEs?

A Differential-Algebraic Equation (DAE) system combines:
- **Differential equations**: `dy/dt = F(y)` — rates of change for species concentrations
- **Algebraic equations**: `G(y) = 0` — constraints that must always be satisfied

In atmospheric chemistry, algebraic constraints are useful for enforcing:
- **Chemical equilibrium** (fast reversible reactions)
- **Mass conservation**
- **Charge balance**

### Why Constraints?

Standard ODE solvers only handle differential equations. When chemical reactions are much faster than the timescale of interest, they effectively reach equilibrium instantaneously. Rather than using extremely small timesteps to resolve these fast reactions, we can enforce the equilibrium condition as an algebraic constraint.

---

## Class Hierarchy

```
Constraint (Abstract Base Class)
    │
    └── EquilibriumConstraint (Concrete Implementation)

ConstraintSet (Container/Manager)
```

**File Locations:**
- `include/micm/constraint/constraint.hpp` — Abstract base class
- `include/micm/constraint/equilibrium_constraint.hpp` — Equilibrium implementation
- `include/micm/constraint/constraint_set.hpp` — Set manager

---

## The Constraint Base Class

The `Constraint` class defines the interface that all constraints must implement.

### Purpose

Constraints define algebraic relations `G(y) = 0` that must be satisfied by the species concentrations. Each constraint provides:
- A **residual function** `G(y)` that equals zero when satisfied
- **Jacobian entries** `dG/dy` for each dependent species

### Class Definition

```cpp
class Constraint
{
 public:
  /// Name of the constraint (for identification/debugging)
  std::string name_;

  /// Names of species this constraint depends on
  std::vector<std::string> species_dependencies_;

  /// Evaluate the constraint residual G(y)
  /// @return Residual value (should be 0 when constraint is satisfied)
  virtual double Residual(const std::vector<double>& concentrations,
                          const std::vector<std::size_t>& indices) const = 0;

  /// Compute partial derivatives dG/d[species] for each dependent species
  /// @return Vector of partial derivatives in same order as species_dependencies_
  virtual std::vector<double> Jacobian(const std::vector<double>& concentrations,
                                       const std::vector<std::size_t>& indices) const = 0;

  /// Get the number of species this constraint depends on
  std::size_t NumberOfDependencies() const
  {
    return species_dependencies_.size();
  }
};
```

### Key Design Points

1. **Pure virtual methods** — `Residual()` and `Jacobian()` must be implemented by derived classes
2. **Species dependencies** — Each constraint declares which species it depends on
3. **Index-based access** — Methods receive indices mapping species names to positions in the concentration vector

---

## The EquilibriumConstraint Class

The `EquilibriumConstraint` class implements equilibrium constraints for reversible chemical reactions.

### Mathematical Background

For a reversible reaction:

```
aA + bB ⇌ cC + dD
```

At equilibrium, the concentrations satisfy:

```
K_eq = [C]^c × [D]^d / ([A]^a × [B]^b)
```

The constraint equation is formulated as:

```
G = K_eq × [A]^a × [B]^b - [C]^c × [D]^d = 0
```

### Class Definition

```cpp
class EquilibriumConstraint : public Constraint
{
 public:
  /// Reactant species names and their stoichiometric coefficients
  std::vector<std::pair<std::string, double>> reactants_;

  /// Product species names and their stoichiometric coefficients
  std::vector<std::pair<std::string, double>> products_;

  /// Equilibrium constant K_eq = k_forward / k_backward
  double equilibrium_constant_;

 private:
  /// Indices for efficient Jacobian computation
  std::vector<std::size_t> reactant_dependency_indices_;
  std::vector<std::size_t> product_dependency_indices_;
};
```

### Constructor

```cpp
EquilibriumConstraint(
    const std::string& name,
    const std::vector<std::pair<std::string, double>>& reactants,
    const std::vector<std::pair<std::string, double>>& products,
    double equilibrium_constant);
```

The constructor:
1. Validates that `K_eq > 0`
2. Builds `species_dependencies_` by concatenating reactants then products
3. Stores index mappings for efficient Jacobian computation

### Residual Calculation

```cpp
double Residual(...) const override
{
  // Compute product of reactant concentrations raised to stoichiometric powers
  double reactant_product = 1.0;
  for (each reactant i)
    reactant_product *= pow([reactant_i], stoich_i);

  // Compute product of product concentrations raised to stoichiometric powers
  double product_product = 1.0;
  for (each product j)
    product_product *= pow([product_j], stoich_j);

  // G = K_eq × [reactants] - [products]
  return equilibrium_constant_ * reactant_product - product_product;
}
```

### Jacobian Calculation

The Jacobian entries are partial derivatives of `G` with respect to each species:

**For reactant R with stoichiometry n:**
```
dG/d[R] = K_eq × n × [R]^(n-1) × prod([other_reactants])
```

**For product P with stoichiometry m:**
```
dG/d[P] = -m × [P]^(m-1) × prod([other_products])
```

The implementation handles special cases when concentration is zero.

### Usage Example

```cpp
// Create an equilibrium constraint: A + B ⇌ AB with K_eq = 1000
EquilibriumConstraint constraint(
    "A_B_equilibrium",
    {{"A", 1.0}, {"B", 1.0}},  // Reactants with stoichiometry
    {{"AB", 1.0}},              // Products
    1000.0);                    // K_eq

// At equilibrium with [A] = 0.001, [B] = 0.001, [AB] = 0.001:
// K_eq × [A] × [B] = 1000 × 0.001 × 0.001 = 0.001 = [AB]
// Residual = 0 ✓
```

---

## The ConstraintSet Class

The `ConstraintSet` class manages a collection of constraints and integrates them with the Rosenbrock solver.

### Purpose

- Owns and manages multiple `Constraint` objects
- Maps species names to variable indices
- Pre-computes Jacobian structure for efficient sparse matrix operations
- Provides batch operations over grid cells

### Class Definition

```cpp
class ConstraintSet
{
 protected:
  struct ConstraintInfo
  {
    std::size_t constraint_index_;       // Index in constraints_ vector
    std::size_t constraint_row_;         // Row in forcing/Jacobian (state_size + i)
    std::size_t number_of_dependencies_; // Number of species dependencies
  };

  std::vector<std::unique_ptr<Constraint>> constraints_;
  std::vector<std::size_t> dependency_ids_;           // Flat list of species indices
  std::vector<ConstraintInfo> constraint_info_;       // Info for each constraint
  std::vector<std::size_t> jacobian_flat_ids_;        // Pre-computed sparse indices
  std::size_t constraint_row_offset_{ 0 };            // = number of species
};
```

### Key Methods

#### 1. Constructor — Maps Species to Indices

```cpp
ConstraintSet(
    std::vector<std::unique_ptr<Constraint>>&& constraints,
    const std::map<std::string, std::size_t>& variable_map,
    std::size_t constraint_row_offset);
```

- Takes ownership of constraints via move semantics
- Maps each species dependency to its index in the state vector
- Throws if a constraint depends on an unknown species

#### 2. NonZeroJacobianElements — Sparse Structure Discovery

```cpp
std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements() const;
```

Returns all (row, column) pairs where constraints contribute to the Jacobian:
- **Row**: `constraint_row_offset + constraint_index`
- **Column**: index of each dependent species

#### 3. SetJacobianFlatIds — Pre-compute Sparse Indices

```cpp
template<typename OrderingPolicy>
void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix);
```

Converts 2D (row, col) indices to 1D flat indices for efficient sparse matrix access during solving.

#### 4. AddForcingTerms — Compute Residuals

```cpp
template<typename DenseMatrixPolicy>
void AddForcingTerms(const DenseMatrixPolicy& state_variables,
                     DenseMatrixPolicy& forcing) const;
```

For each constraint `G_i`, adds `G_i(x)` to `forcing[constraint_row_offset + i]`.

#### 5. SubtractJacobianTerms — Compute Jacobian Contributions

```cpp
template<class DenseMatrixPolicy, class SparseMatrixPolicy>
void SubtractJacobianTerms(const DenseMatrixPolicy& state_variables,
                           SparseMatrixPolicy& jacobian) const;
```

**Subtracts** constraint Jacobian entries (matching ProcessSet convention for Rosenbrock solver).

---

## How ConstraintSet Mirrors ProcessSet

The `ConstraintSet` class is intentionally designed to mirror the `ProcessSet` class, ensuring consistent integration with the solver.

### Architectural Comparison

| Aspect | ProcessSet | ConstraintSet |
|--------|-----------|---------------|
| **Purpose** | Collection of chemical reactions | Collection of algebraic constraints |
| **Base Unit** | `Process` (ChemicalReaction) | `Constraint` |
| **Container** | `vector<Process>` | `vector<unique_ptr<Constraint>>` |
| **Dependency Tracking** | `reactant_ids_`, `product_ids_` | `dependency_ids_` |
| **Sparse Index Caching** | `jacobian_flat_ids_` | `jacobian_flat_ids_` |
| **Info Structure** | `ProcessInfo` | `ConstraintInfo` |

### Method-Level Comparison

Both classes implement the same key methods with matching signatures:

#### NonZeroJacobianElements

```cpp
// ProcessSet
std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements() const;

// ConstraintSet
std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements() const;
```

Both return the sparsity pattern for their contributions to the Jacobian.

#### SetJacobianFlatIds

```cpp
// ProcessSet
template<typename OrderingPolicy>
void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix);

// ConstraintSet
template<typename OrderingPolicy>
void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix);
```

Both pre-compute flat indices for efficient sparse matrix access.

#### AddForcingTerms

```cpp
// ProcessSet — adds reaction rates to forcing (dy/dt contributions)
template<typename DenseMatrixPolicy>
void AddForcingTerms(const DenseMatrixPolicy& rate_constants,
                     const DenseMatrixPolicy& state_variables,
                     DenseMatrixPolicy& forcing) const;

// ConstraintSet — adds constraint residuals to forcing (G(y) contributions)
template<typename DenseMatrixPolicy>
void AddForcingTerms(const DenseMatrixPolicy& state_variables,
                     DenseMatrixPolicy& forcing) const;
```

Note: ConstraintSet doesn't need `rate_constants` — constraints are evaluated directly from concentrations.

#### SubtractJacobianTerms

```cpp
// ProcessSet — subtracts dF/dy from Jacobian
template<class DenseMatrixPolicy, class SparseMatrixPolicy>
void SubtractJacobianTerms(const DenseMatrixPolicy& rate_constants,
                           const DenseMatrixPolicy& state_variables,
                           SparseMatrixPolicy& jacobian) const;

// ConstraintSet — subtracts dG/dy from Jacobian
template<class DenseMatrixPolicy, class SparseMatrixPolicy>
void SubtractJacobianTerms(const DenseMatrixPolicy& state_variables,
                           SparseMatrixPolicy& jacobian) const;
```

### Grid Cell Loop Pattern

Both classes iterate over grid cells in the same way:

```cpp
// ProcessSet pattern
for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
{
  auto cell_rate_constants = rate_constants[i_cell];
  auto cell_state = state_variables[i_cell];
  auto cell_forcing = forcing[i_cell];
  // ... process reactions ...
}

// ConstraintSet pattern
for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
{
  auto cell_state = state_variables[i_cell];
  auto cell_forcing = forcing[i_cell];
  // ... evaluate constraints ...
}
```

### Row Offset Distinction

The key structural difference:
- **ProcessSet** operates on rows `0` to `N_species - 1`
- **ConstraintSet** operates on rows `N_species` to `N_species + N_constraints - 1`

This is managed by `constraint_row_offset_`:

```cpp
// In ConstraintInfo
info.constraint_row_ = constraint_row_offset_ + i;  // e.g., if 3 species, first constraint is row 3
```

---

## Integration with the DAE System

### Extended State Vector

The DAE solver extends the state to include constraint variables:

```
State = [y_0, y_1, ..., y_{N-1}, g_0, g_1, ..., g_{M-1}]
         ├── Species (N) ──┤    ├── Constraints (M) ─┤
```

### Extended Jacobian Structure

```
           y_0  y_1  ...  y_N  g_0  g_1  ...
    y_0  [  ×    ×         ×                  ]  ← df_0/dy
    y_1  [  ×    ×         ×                  ]
    ...
    g_0  [  ×    ×         ×                  ]  ← dG_0/dy (ConstraintSet rows)
    g_1  [  ×    ×         ×                  ]
```

### Forcing Vector Structure

```
Forcing = [dy_0/dt, dy_1/dt, ..., G_0(y), G_1(y), ...]
           ├─ ProcessSet ─┤    ├─ ConstraintSet ─┤
```

### Solver Integration Flow

```
1. Build Jacobian sparsity pattern:
   - process_set.NonZeroJacobianElements()
   - constraint_set.NonZeroJacobianElements()
   - Combine and create sparse matrix

2. Pre-compute flat indices:
   - process_set.SetJacobianFlatIds(jacobian)
   - constraint_set.SetJacobianFlatIds(jacobian)

3. Each solver iteration:
   a. Zero forcing vector
   b. process_set.AddForcingTerms(...)      // F(y) rows
   c. constraint_set.AddForcingTerms(...)   // G(y) rows

   d. Zero Jacobian
   e. process_set.SubtractJacobianTerms(...)      // dF/dy entries
   f. constraint_set.SubtractJacobianTerms(...)   // dG/dy entries

   g. Solve linear system
```

---

## Design Patterns Summary

1. **Abstract Factory** — `Constraint` base class allows different constraint types
2. **Move Semantics** — `ConstraintSet` takes ownership via `unique_ptr` move
3. **Flat Array Optimization** — Pre-computed indices for cache-friendly iteration
4. **Template Policy** — Supports various dense/sparse matrix implementations
5. **Parallel Structure** — Mirrors `ProcessSet` for consistent solver integration

This architecture enables clean extension of MICM's Rosenbrock solver from ODEs to DAEs while maintaining performance and code consistency.

---

## Note on DAE Notation Conventions

In DAE literature, you may encounter a distinction between **differential variables** (y) and **algebraic variables** (z), written as a semi-explicit index-1 DAE:

```
dy/dt = F(y, z)
    0 = G(y, z)
```

Here, y evolves according to differential equations while z is determined implicitly by the algebraic constraints. The full system Jacobian then has a block structure:

```
         ∂/∂y      ∂/∂z
       ┌─────────┬─────────┐
 F     │  ∂F/∂y  │  ∂F/∂z  │   ← Differential equations
       ├─────────┼─────────┤
 G     │  ∂G/∂y  │  ∂G/∂z  │   ← Algebraic constraints
       └─────────┴─────────┘
```

MICM uses a **unified state vector y** for simplicity, where all species concentrations (whether governed by ODEs or determined by constraints) are stored together. This reflects how the code actually works — `ProcessSet` and `ConstraintSet` both operate on the same state vector. The distinction between "differential" and "algebraic" species is implicit in which equations govern them, rather than explicit in the notation.
