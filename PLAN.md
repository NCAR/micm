# Rosenbrock DAE Solver Implementation Plan

## Executive Summary

This plan details the implementation of algebraic constraint support for the Rosenbrock solver in MICM. The DAE formulation allows solving systems with both differential equations (dx/dt = F(x,y,t)) and algebraic constraints (G(x,y) = 0), enabling applications like equilibrium chemistry and conservation laws.

## 1. Mathematical Foundation

### 1.1 DAE System Structure

From `ODE.wiki/Rosenbrock-Methods.md`:

```
State vector: [x; y] where
  - x: differential variables (ODEs): dx/dt = F(x, y, t)
  - y: algebraic variables (constraints): G(x, y) = 0

Block Jacobian structure:
[F_x  F_y]   where F_x = dF/dx, F_y = dF/dy
[G_x  G_y]         G_x = dG/dx, G_y = dG/dy

Rosenbrock stage system:
[I - hγF_x    -hγF_y ] [k]     [F]
[  -hγG_x     -hγG_y ] [l]  =  [G] + correction terms
```

**Key insight**: The upper-left identity block (I) applies only to ODE rows, not constraint rows. This is already implemented via `upper_left_identity_diagonal_` (1.0 for ODEs, 0.0 for constraints).

### 1.2 Stiffly-Accurate Property

The RODAS methods (`FourStageDifferentialAlgebraicRosenbrockParameters()` and `SixStageDifferentialAlgebraicRosenbrockParameters()`) are "stiffly accurate" meaning the last stage equals the solution, ensuring consistent algebraic constraints at step completion.

### 1.3 KPP Reference

KPP provides ROS-2, ROS-3, ROS-4, RODAS-3, RODAS-4 variants targeting relative errors ~10^-2 to 10^-5 for atmospheric chemistry.

## 2. Design Decision: Species-Only Constraints

Constraints operate **only on existing species concentrations** without adding new algebraic variables to the state vector. This means:
- No Lagrange multipliers or auxiliary constraint variables
- Constraints enforce relations like equilibrium ratios directly on species
- State size remains equal to number of species
- Simpler implementation and clearer semantics

## 3. Current State (rosenbrock_dae branch)

The branch already has infrastructure in place:
- `StateParameters::number_of_constraints_` field
- `State::constraint_size_` and `upper_left_identity_diagonal_`
- `AlphaMinusJacobian` uses the identity mask correctly
- DAE Rosenbrock parameter sets exist
- `SolverBuilder::SetConstraintCount()` and `SetConstraintNames()` methods

**What's missing**:
- Constraint equation specification API
- Constraint forcing term computation (G(x,y))
- Constraint Jacobian term computation (G_x, G_y)
- Integration of constraint contributions into the solver loop
- Test cases for DAE systems

## 3. API Design

### 3.1 Constraint Base Class (New)

```cpp
// include/micm/constraint/constraint.hpp
namespace micm
{
  class Constraint
  {
   public:
    std::string name_;
    std::vector<std::string> variable_dependencies_;

    virtual ~Constraint() = default;

    /// @brief Evaluate constraint residual G(x,y)
    virtual double Residual(const std::map<std::string, double>& concentrations) const = 0;

    /// @brief Compute partial derivatives dG/d[species]
    virtual std::map<std::string, double> Jacobian(
        const std::map<std::string, double>& concentrations) const = 0;
  };
}
```

### 3.2 EquilibriumConstraint (Built-in)

For reversible reaction equilibrium: K_eq = [products]/[reactants]

```cpp
// include/micm/constraint/equilibrium_constraint.hpp
class EquilibriumConstraint : public Constraint
{
 public:
    std::vector<std::pair<std::string, double>> reactants_;   // species, stoichiometry
    std::vector<std::pair<std::string, double>> products_;    // species, stoichiometry
    double equilibrium_constant_;

    double Residual(const std::map<std::string, double>& conc) const override
    {
        // G = k_f * prod([reactants]) - k_b * prod([products]) = 0
    }
};
```

### 3.3 ConservationConstraint (Built-in)

For mass/element conservation: sum(c_i * [X_i]) = constant

```cpp
// include/micm/constraint/conservation_constraint.hpp
class ConservationConstraint : public Constraint
{
 public:
    std::vector<std::pair<std::string, double>> species_coefficients_;
    double conserved_total_;
};
```

### 3.4 ConstraintSet (Analogous to ProcessSet)

```cpp
// include/micm/constraint/constraint_set.hpp
class ConstraintSet
{
 public:
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements() const;

    template<typename OrderingPolicy>
    void SetJacobianFlatIds(const SparseMatrix<double, OrderingPolicy>& matrix,
                            std::size_t constraint_row_offset);

    template<typename DenseMatrixPolicy>
    void AddForcingTerms(const DenseMatrixPolicy& state_variables,
                         DenseMatrixPolicy& forcing,
                         std::size_t constraint_row_offset) const;

    template<class DenseMatrixPolicy, class SparseMatrixPolicy>
    void SubtractJacobianTerms(const DenseMatrixPolicy& state_variables,
                               SparseMatrixPolicy& jacobian,
                               std::size_t constraint_row_offset) const;
};
```

### 3.5 Updated SolverBuilder API

```cpp
template<...>
class SolverBuilder
{
 public:
    SolverBuilder& SetConstraints(std::vector<std::unique_ptr<Constraint>>&& constraints);
    SolverBuilder& AddConstraint(std::unique_ptr<Constraint> constraint);
};
```

## 4. Implementation Details

### 4.1 New Files to Create

| File | Purpose |
|------|---------|
| `include/micm/constraint/constraint.hpp` | Abstract base class |
| `include/micm/constraint/equilibrium_constraint.hpp` | Equilibrium K_eq constraint |
| `include/micm/constraint/conservation_constraint.hpp` | Mass conservation constraint |
| `include/micm/constraint/constraint_set.hpp` | Collection with forcing/Jacobian |
| `include/micm/constraint/constraint_set.inl` | Template implementations |
| `test/unit/constraint/test_equilibrium_constraint.cpp` | Unit tests |
| `test/unit/constraint/test_constraint_set.cpp` | Unit tests |
| `test/integration/test_analytical_dae_rosenbrock.cpp` | Integration tests |

### 4.2 Files to Modify

| File | Modifications |
|------|---------------|
| `include/micm/solver/rosenbrock.hpp` | Add ConstraintSet member |
| `include/micm/solver/rosenbrock.inl` | Call constraint forcing/Jacobian in Solve() |
| `include/micm/solver/solver_builder.hpp` | Add constraint storage |
| `include/micm/solver/solver_builder.inl` | Build ConstraintSet, merge Jacobian |

### 4.3 Solver Loop Modifications

```cpp
// rosenbrock.inl - In Solve() method

// After ODE forcing:
initial_forcing.Fill(0);
rates_.AddForcingTerms(state.rate_constants_, Y, initial_forcing);

// NEW: Add constraint residuals (constraint rows)
if (state.constraint_size_ > 0)
{
    constraints_.AddForcingTerms(Y, initial_forcing, state.state_size_);
}

// After ODE Jacobian:
state.jacobian_.Fill(0);
rates_.SubtractJacobianTerms(state.rate_constants_, Y, state.jacobian_);

// NEW: Add constraint Jacobian terms
if (state.constraint_size_ > 0)
{
    constraints_.SubtractJacobianTerms(Y, state.jacobian_, state.state_size_);
}
```

### 4.4 Jacobian Structure

```
Species:    s0  s1  s2  c0  c1    (si = species, ci = constraint vars)
Row 0 (s0): [x   x   .   x   .]   <- dF0/d(all)
Row 1 (s1): [x   x   x   .   .]   <- dF1/d(all)
Row 2 (s2): [.   x   x   .   x]   <- dF2/d(all)
Row 3 (c0): [x   x   .   x   .]   <- dG0/d(all)  CONSTRAINT ROW
Row 4 (c1): [.   .   x   .   x]   <- dG1/d(all)  CONSTRAINT ROW

upper_left_identity_diagonal_: [1, 1, 1, 0, 0]
```

## 5. Test Strategy

### 5.1 Preserve Existing Tests

All 52 existing tests must pass unchanged. The constraint system is additive.

### 5.2 Key Test: Reversible Reaction Equilibrium

**System**: A + B <-> AB

- Forward rate: k_f = 1e6 M^-1 s^-1
- Backward rate: k_b = 1e3 s^-1
- K_eq = k_f/k_b = 1000 M^-1

**Initial conditions**: [A]_0 = [B]_0 = 1.0 M, [AB]_0 = 0

**Analytical equilibrium solution**:
For equal initial concentrations, solve K_eq = x/(1-x)^2 where x = [AB]_eq:
- [AB]_eq ≈ 0.96875 M
- [A]_eq = [B]_eq ≈ 0.03125 M

**Verification**:
1. K_eq = [AB]/([A][B]) = 1000 at equilibrium
2. Mass conservation: [A] + [AB] = 1.0, [B] + [AB] = 1.0

```cpp
TEST(DAERosenbrock, ReversibleReactionEquilibrium)
{
    // A + B <-> AB with K_eq = 1000
    auto equilibrium = std::make_unique<EquilibriumConstraint>(
        {{"A", 1.0}, {"B", 1.0}},  // reactants
        {{"AB", 1.0}},             // products
        1000.0                      // K_eq
    );

    auto solver = builder
        .SetSystem(system)
        .SetReactions({forward, backward})
        .AddConstraint(std::move(equilibrium))
        .Build();

    // Solve and verify K_eq = [AB]/([A][B]) = 1000
}
```

### 5.3 Additional Tests

1. **Conservation constraint**: Verify sum([A] + [B] + [C]) = constant
2. **ODE-only vs DAE**: Same results when no constraints added
3. **AlphaMinusJacobian**: Verify alpha NOT added to constraint diagonals
4. **Vectorized solvers**: Test with VectorMatrix types

## 6. Implementation Sequence

### Phase 1: Core Constraint Infrastructure
1. Create `Constraint` base class
2. Create `EquilibriumConstraint`
3. Create `ConstraintSet` with forcing/Jacobian methods
4. Add unit tests

### Phase 2: Solver Integration
5. Modify `SolverBuilder` to accept constraints
6. Modify solver constructor to include `ConstraintSet`
7. Update `Solve()` loop
8. Verify `AlphaMinusJacobian` works correctly

### Phase 3: Testing
9. Run all 52 existing tests (must pass)
10. Add constraint unit tests
11. Implement reversible reaction test
12. Add vectorized solver tests

### Phase 4: Polish
13. Add `ConservationConstraint`
14. Documentation/examples
15. Performance optimization if needed

## 7. Critical Files

- `include/micm/solver/rosenbrock.inl` - Core solver loop
- `include/micm/solver/solver_builder.inl` - Builder logic
- `include/micm/solver/state.hpp` - Reference for identity mask usage
- `include/micm/process/process_set.hpp` - Pattern for ConstraintSet
- `test/integration/analytical_policy.hpp` - Pattern for tests

## 8. Verification Checklist

- [ ] All 52 existing tests pass
- [ ] Unit tests for Constraint classes
- [ ] Unit tests for ConstraintSet
- [ ] Reversible reaction reaches correct equilibrium
- [ ] Mass conservation verified
- [ ] AlphaMinusJacobian uses identity mask correctly
- [ ] Vectorized solver variants work
- [ ] No performance regression for ODE-only mode
