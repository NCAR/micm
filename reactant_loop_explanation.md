# Line-by-Line Explanation: Reactant Product Loop

This document explains the indexing logic in the `EquilibriumConstraint::Residual()` method.

## The Code

```cpp
// Compute product of reactant concentrations raised to stoichiometric powers
double reactant_product = 1.0;
for (std::size_t i = 0; i < reactants_.size(); ++i)
{
  double conc = concentrations[indices[reactant_dependency_indices_[i]]];
  reactant_product *= std::pow(conc, reactants_[i].second);
}
```

## Context

For an equilibrium constraint like `A + 2B ⇌ C`, we need to compute:

```
K_eq × [A]^1 × [B]^2
```

This loop computes the `[A]^1 × [B]^2` part (the product of reactant concentrations raised to their stoichiometric powers).

---

## Line-by-Line Breakdown

### Line 1: Comment
```cpp
// Compute product of reactant concentrations raised to stoichiometric powers
```
Describes the purpose: multiply together each reactant's concentration raised to its stoichiometry.

---

### Line 2: Initialize accumulator
```cpp
double reactant_product = 1.0;
```
Start with 1.0 because we're multiplying. Each iteration multiplies in another term.

---

### Line 3: Loop over reactants
```cpp
for (std::size_t i = 0; i < reactants_.size(); ++i)
```
Iterate over each reactant in the reaction. For `A + 2B ⇌ C`:
- `i = 0` → reactant A
- `i = 1` → reactant B

---

### Line 5: The triple indirection to get concentration
```cpp
double conc = concentrations[indices[reactant_dependency_indices_[i]]];
```

This is the key line. Let's unpack it from inside out:

#### Step 1: `reactant_dependency_indices_[i]`

This maps from "which reactant" to "which position in `species_dependencies_`".

**Why is this needed?** The `species_dependencies_` vector contains all species (reactants first, then products). For the constraint `A + 2B ⇌ C`:

```
species_dependencies_ = ["A", "B", "C"]
                          0    1    2
```

The `reactant_dependency_indices_` tells us where each reactant lives in that list:
```
reactant_dependency_indices_ = [0, 1]  // A is at index 0, B is at index 1
```

So `reactant_dependency_indices_[0]` = 0 (A's position), `reactant_dependency_indices_[1]` = 1 (B's position).

#### Step 2: `indices[...]`

The `indices` parameter maps from "position in `species_dependencies_`" to "position in the global `concentrations` vector".

**Why is this needed?** The solver has a global state vector with ALL species in the system, not just the ones this constraint cares about. For example, if the system has species `[X, A, Y, B, C, Z]`:

```
concentrations = [X_conc, A_conc, Y_conc, B_conc, C_conc, Z_conc]
                    0       1       2       3       4       5
```

The `indices` vector tells us where our dependent species are in this global vector:
```
indices = [1, 3, 4]  // A is at global index 1, B at 3, C at 4
```

So `indices[0]` = 1 (A's global position), `indices[1]` = 3 (B's global position).

#### Step 3: `concentrations[...]`

Finally, we index into the global concentrations vector to get the actual concentration value.

#### Full example trace for reactant A (i=0):

```
reactant_dependency_indices_[0] = 0     // A is the 0th species dependency
indices[0] = 1                          // A is at global index 1
concentrations[1] = A_conc              // Get A's concentration
```

#### Full example trace for reactant B (i=1):

```
reactant_dependency_indices_[1] = 1     // B is the 1st species dependency
indices[1] = 3                          // B is at global index 3
concentrations[3] = B_conc              // Get B's concentration
```

---

### Line 6: Multiply with stoichiometry
```cpp
reactant_product *= std::pow(conc, reactants_[i].second);
```

- `reactants_[i]` is a `std::pair<std::string, double>` where:
  - `.first` = species name (e.g., "A")
  - `.second` = stoichiometric coefficient (e.g., 1.0 for A, 2.0 for B)

- `std::pow(conc, reactants_[i].second)` raises the concentration to the stoichiometric power

- `reactant_product *= ...` accumulates the product

For `A + 2B`:
```
Iteration 0: reactant_product = 1.0 × [A]^1 = [A]
Iteration 1: reactant_product = [A] × [B]^2 = [A] × [B]²
```

---

## Visual Summary

```
                    reactants_
                   ┌─────────────────┐
                   │ ("A", 1.0)      │  i=0
                   │ ("B", 2.0)      │  i=1
                   └─────────────────┘
                           │
                           ▼
            reactant_dependency_indices_
                   ┌─────────────────┐
                   │  0              │  → position of A in species_dependencies_
                   │  1              │  → position of B in species_dependencies_
                   └─────────────────┘
                           │
                           ▼
               species_dependencies_
                   ┌─────────────────┐
                   │ "A"             │  index 0
                   │ "B"             │  index 1
                   │ "C"             │  index 2 (product)
                   └─────────────────┘
                           │
                           ▼
                      indices
                   ┌─────────────────┐
                   │  1              │  → A is at global index 1
                   │  3              │  → B is at global index 3
                   │  4              │  → C is at global index 4
                   └─────────────────┘
                           │
                           ▼
                  concentrations
          ┌──────────────────────────────┐
          │ [X]  [A]  [Y]  [B]  [C]  [Z] │
          │  0    1    2    3    4    5  │
          └──────────────────────────────┘
```

## Complete State vs Constraint-Specific Vectors

Understanding which vectors hold the complete system state versus only the species in this constraint:

**Complete state (all species in the system):**
- `concentrations` — the parameter passed into `Residual()`/`Jacobian()`, contains every species in the solver

**Constraint-specific (only species in this constraint):**
- `reactants_` — just the reactants (name + stoichiometry pairs)
- `products_` — just the products (name + stoichiometry pairs)
- `species_dependencies_` — reactants + products concatenated (names only)
- `reactant_dependency_indices_` — positions of reactants within `species_dependencies_`
- `product_dependency_indices_` — positions of products within `species_dependencies_`

**The bridge between them:**
- `indices` — the parameter that maps constraint-specific positions → global positions

```
Constraint-specific           Bridge              Complete state
─────────────────────        ────────            ──────────────
reactant_dependency_indices_    │
         │                      │
         ▼                      │
species_dependencies_  ──────► indices ──────► concentrations
         ▲                      │
         │                      │
product_dependency_indices_     │
```

The `indices` vector has the same length as `species_dependencies_` (constraint-specific), but its values are positions into `concentrations` (complete state).

---

## Why This Indirection?

The triple indirection exists because:

1. **Constraints are self-contained** — They store their reactants/products by name
2. **The solver uses indices** — For performance, the solver works with integer indices, not string lookups
3. **Species ordering varies** — Different systems have different species in different orders
4. **Separation of concerns** — The constraint doesn't need to know about the global system layout; the `indices` parameter handles the mapping at runtime
