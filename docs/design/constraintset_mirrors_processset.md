# Parallel Structure: ProcessSet vs ConstraintSet Jacobian Computation

This document shows how `ConstraintSet::SubtractJacobianTerms` mirrors `ProcessSet::SubtractJacobianTerms`, following the same structural pattern for consistent solver integration.

## Method Signatures

```cpp
// ProcessSet: Jacobian contributions from chemical reactions (dF/dy)
template<class DenseMatrixPolicy, class SparseMatrixPolicy>
void ProcessSet::SubtractJacobianTerms(
    const DenseMatrixPolicy& rate_constants,
    const DenseMatrixPolicy& state_variables,
    SparseMatrixPolicy& jacobian) const;

// ConstraintSet: Jacobian contributions from constraints (dG/dy)
template<class DenseMatrixPolicy, class SparseMatrixPolicy>
void ConstraintSet::SubtractJacobianTerms(
    const DenseMatrixPolicy& state_variables,
    SparseMatrixPolicy& jacobian) const;
```

Note: ConstraintSet doesn't need `rate_constants` — constraints are evaluated directly from concentrations.

---

## Side-by-Side Code Comparison

### ProcessSet::SubtractJacobianTerms

```cpp
auto cell_jacobian = jacobian.AsVector().begin();

// loop over grid cells
for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
{
  auto cell_rate_constants = rate_constants[i_cell];
  auto cell_state = state_variables[i_cell];

  auto react_id = jacobian_reactant_ids_.begin();
  auto yield = jacobian_yields_.begin();
  auto flat_id = jacobian_flat_ids_.begin();

  // loop over process-dependent variable pairs
  for (const auto& process_info : jacobian_process_info_)
  {
    double d_rate_d_ind = cell_rate_constants[process_info.process_id_];
    for (std::size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
      d_rate_d_ind *= cell_state[*(react_id++)];
    for (std::size_t i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
      cell_jacobian[*(flat_id++)] += d_rate_d_ind;
    for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
      cell_jacobian[*(flat_id++)] -= *(yield++) * d_rate_d_ind;
  }
  // increment cell_jacobian after each grid cell
  cell_jacobian += jacobian.FlatBlockSize();
}
```

### ConstraintSet::SubtractJacobianTerms

```cpp
// Allocate reusable buffer for constraint Jacobian values
std::vector<double> jac_buffer(max_dependencies_);

auto cell_jacobian = jacobian.AsVector().begin();

// Loop over grid cells
for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
{
  auto cell_state = state_variables[i_cell];

  // Get pointer to concentration data for this cell (no copy!)
  const double* concentrations = &cell_state[0];

  for (const auto& info : constraint_info_)
  {
    // Get pointer to pre-computed indices for this constraint
    const std::size_t* indices = dependency_ids_.data() + info.dependency_offset_;

    // Compute constraint Jacobian into reusable buffer
    constraints_[info.constraint_index_]->Jacobian(concentrations, indices, jac_buffer.data());

    // Get pointer to pre-computed flat indices for this constraint
    const std::size_t* flat_ids = jacobian_flat_ids_.data() + info.jacobian_flat_offset_;

    // Subtract Jacobian entries (matching ProcessSet convention)
    for (std::size_t i = 0; i < info.number_of_dependencies_; ++i)
    {
      cell_jacobian[flat_ids[i]] -= jac_buffer[i];
    }
  }

  // Advance to next grid cell's Jacobian block
  cell_jacobian += jacobian.FlatBlockSize();
}
```

---

## Structural Parallel

| Step | ProcessSet | ConstraintSet |
|------|-----------|---------------|
| 1. Get Jacobian iterator | `auto cell_jacobian = jacobian.AsVector().begin();` | `auto cell_jacobian = jacobian.AsVector().begin();` |
| 2. Grid cell loop | `for (i_cell = 0; i_cell < NumRows(); ++i_cell)` | `for (i_cell = 0; i_cell < NumRows(); ++i_cell)` |
| 3. Get cell state | `auto cell_state = state_variables[i_cell];` | `auto cell_state = state_variables[i_cell];` |
| 4. Init flat index iterator | `auto flat_id = jacobian_flat_ids_.begin();` | `auto flat_id = jacobian_flat_ids_.begin();` |
| 5. Info loop | `for (const auto& process_info : jacobian_process_info_)` | `for (const auto& info : constraint_info_)` |
| 6. Compute Jacobian values | Inline calculation of `d_rate_d_ind` | Call `Jacobian(concentrations, indices)` |
| 7. Store to sparse matrix | `cell_jacobian[*(flat_id++)] += / -=` | `cell_jacobian[*(flat_id++)] -=` |
| 8. Advance to next cell | `cell_jacobian += jacobian.FlatBlockSize();` | `cell_jacobian += jacobian.FlatBlockSize();` |

---

## Key Parallels

### 1. Same Sparse Matrix Access Pattern

Both use pre-computed flat indices (`jacobian_flat_ids_`) to write directly into the sparse Jacobian:

```cpp
// ProcessSet
cell_jacobian[*(flat_id++)] += d_rate_d_ind;

// ConstraintSet
cell_jacobian[*(flat_id++)] -= jac[i];
```

### 2. Same Grid Cell Iteration

Both loop over grid cells and advance the Jacobian pointer by `FlatBlockSize()`:

```cpp
// Both methods
for (std::size_t i_cell = 0; i_cell < state_variables.NumRows(); ++i_cell)
{
  // ... process cell ...
  cell_jacobian += jacobian.FlatBlockSize();
}
```

### 3. Same Info Structure Pattern

Both use a nested struct to track per-item metadata:

```cpp
// ProcessSet
struct ProcessInfo {
  std::size_t process_id_;
  std::size_t independent_id_;
  std::size_t number_of_dependent_reactants_;
  std::size_t number_of_products_;
};

// ConstraintSet
struct ConstraintInfo {
  std::size_t constraint_index_;
  std::size_t constraint_row_;
  std::size_t number_of_dependencies_;
};
```

### 4. Same Flat Index Pre-computation

Both use `SetJacobianFlatIds()` to pre-compute sparse matrix indices:

```cpp
// ProcessSet
template<typename OrderingPolicy>
void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix);

// ConstraintSet
template<typename OrderingPolicy>
void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix);
```

---

## Key Differences

| Aspect | ProcessSet | ConstraintSet |
|--------|-----------|---------------|
| **Input data** | Needs `rate_constants` + `state_variables` | Only needs `state_variables` |
| **Jacobian rows** | Species rows (0 to N-1) | Constraint rows (N to N+M-1) |
| **Jacobian calculation** | Inline: `dF_d_ind = k * [reactants]` | Delegated to `Constraint::Jacobian()` |
| **Entry signs** | `+=` for reactants, `-=` for products | `-=` for all (matches solver convention) |
| **Complexity** | Handles reactants and products separately | Uniform handling of all dependencies |

---

## Visual: Where Each Writes in the Jacobian

```
Extended Jacobian Matrix:

         species columns
       ┌─────────────────┐
       │  ×   ×   ×   ×  │  ← species row 0    ┐
       │  ×   ×   ×   ×  │  ← species row 1    │ ProcessSet writes here
       │  ×   ×   ×   ×  │  ← species row 2    │ (dF/dy entries)
       │  ×   ×   ×   ×  │  ← species row 3    ┘
       ├─────────────────┤
       │  ×   ×   ×   ×  │  ← constraint row 0 ┐ ConstraintSet writes here
       │  ×   ×   ×   ×  │  ← constraint row 1 ┘ (dG/dy entries)
       └─────────────────┘
```

Both contribute to the same sparse Jacobian matrix, but to different row regions.

---

## Why Mirror the Pattern?

1. **Solver compatibility** — The Rosenbrock solver calls both methods the same way
2. **Sparse matrix efficiency** — Both benefit from pre-computed flat indices
3. **Grid parallelism** — Both support multi-cell computation with the same loop structure
4. **Maintainability** — Developers familiar with ProcessSet immediately understand ConstraintSet
5. **Future vectorization** — ConstraintSet can follow ProcessSet's vectorized implementation pattern

---

## Mathematical Context

**ProcessSet** computes derivatives of the forcing function F(y):
```
For reaction: A + B → C
Rate = k[A][B]

∂(d[A]/dt)/∂[A] = -k[B]      (reactant loses)
∂(d[A]/dt)/∂[B] = -k[A]      (reactant loses)
∂(d[C]/dt)/∂[A] = +k[B]      (product gains)
∂(d[C]/dt)/∂[B] = +k[A]      (product gains)
```

**ConstraintSet** computes derivatives of constraint residuals:
```
For constraint: G = K_eq[A][B] - [C] = 0

∂G/∂[A] = K_eq[B]
∂G/∂[B] = K_eq[A]
∂G/∂[C] = -1
```

Both contribute partial derivatives to the system Jacobian, just for different equation types (differential vs algebraic).
