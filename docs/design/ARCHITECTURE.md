# MICM Rosenbrock + DAE Architecture

## 1. Purpose
This document explains the mathematical model and implementation architecture for MICM's Rosenbrock solver with algebraic constraints.

It covers:
- the ODE and DAE equations being solved,
- how constraints are represented,
- how Jacobian/forcing terms are assembled,
- how full mass-matrix DAE enforcement is implemented,
- the solver flow at runtime,
- vectorized/non-vectorized implementation details,
- current limitations and extension points.

The focus is the CPU Rosenbrock path and constraint integration used by `SolverBuilder`.

## 2. Mathematical Model

### 2.1 Kinetic ODE model
For species state vector `y` (length `n_species`), kinetic chemistry is modeled as

`dy/dt = f_chem(y)`

where `f_chem` is assembled from reaction rates and stoichiometric yields.

### 2.2 Algebraic constraints
For constrained systems, one or more algebraic equations are added:

`g_k(y) = 0`, for `k = 1..m`

Current implemented constraint type:
- `EquilibriumConstraint`, with residual form
  `G = K_eq * product(reactants^stoich) - product(products^stoich)`.

### 2.3 Mass-matrix DAE form
The enforced system is represented as:

`M * dy/dt = F(y)`

where:
- `M` is diagonal,
- `M_ii = 1` for ODE species rows,
- `M_ii = 0` for algebraic rows,
- `F_i(y) = f_chem,i(y)` for ODE rows,
- `F_i(y) = g_i(y)` for algebraic rows (row-replacement mode).

Important: constrained variables are still part of the species state. Constraints replace selected species rows; they are not separate physical unknowns.

## 3. Rosenbrock Discretization (as implemented)

### 3.1 Linear system form
At each accepted/rejected internal step, the factored matrix is:

`A = alpha * M - J`

where:
- `alpha = 1 / (H * gamma_0)`,
- `H` is internal time step,
- `J = dF/dy` (with sign convention in assembly matching existing ProcessSet/ConstraintSet behavior).

### 3.2 Stage equation structure
Stages use right-hand-side assembly with coefficients `a_ij`, `c_ij`, `m_i`, `e_i` from selected Rosenbrock parameters.

A key DAE detail: the `c_ij / H` coupling is multiplied by `M` row-wise in implementation, i.e. algebraic rows (where `M_ii=0`) do not receive stiff `1/H` coupling terms.

### 3.3 Error control
Adaptive step size uses embedded error estimate:
- `Y_error = sum(e_i * K_i)`
- `error = normalized_norm(Y, Y_new, Y_error)`
- `H_new = H * clamp(...)`

Error normalization is computed over species-state dimensions (row/column indexing), avoiding layout-dependent flattening errors.

## 4. Data Model and Dimensions

### 4.1 State dimensions
`State` stores:
- `state_size_` = number of species,
- `constraint_size_` = number of constraints (metadata),
- `variables_` shape = `(n_cells, state_size_)`.

In full row-replacement mode:
- Jacobian size is `state_size_ x state_size_` (not `state_size_ + constraint_size_`),
- temporary Rosenbrock matrices (`Ynew_`, `initial_forcing_`, `K_`, `Yerror_`) are species-sized,
- mass-matrix diagonal length is `state_size_`.

### 4.2 Mass matrix storage
Mass diagonal is stored in `state.upper_left_identity_diagonal_` and populated from `StateParameters.mass_matrix_diagonal_`.

This vector is used in:
- `AlphaMinusJacobian`: diagonal shift by `alpha * M_ii`,
- stage coupling: `c_ij/H * M_ii * K_j`.

## 5. Constraint Representation

### 5.1 Constraint variant
`Constraint` is a variant wrapper currently over `EquilibriumConstraint`.

Each constraint provides:
- residual evaluation,
- Jacobian entries with respect to dependent species,
- target algebraic species via `GetAlgebraicSpecies()`.

### 5.2 Algebraic row selection
For `EquilibriumConstraint`, the current mapping policy is:
- algebraic row = first product species.

For example:
- constraint `K_eq * [B] - [C] = 0` maps to algebraic species row `C`.

## 6. Build-Time Wiring (`SolverBuilder`)

`SolverBuilder::Build` does the following:

1. Build kinetic `RatesPolicy` (`ProcessSet`) and kinetic Jacobian sparsity.
2. Build `ConstraintSet` in row-replacement mode when constraints are present.
3. Collect algebraic variable IDs from constraints.
4. Pass algebraic IDs to `ProcessSet` so kinetic forcing/Jacobian does not write those rows.
5. Merge constraint Jacobian sparsity into overall Jacobian sparsity.
6. Build Jacobian with species-only dimension (`n_species`).
7. Build mass-matrix diagonal with zeros at algebraic row IDs.
8. Store state parameters including:
   - `mass_matrix_diagonal_`,
   - `constraints_replace_state_rows_ = true` when constraints exist.

## 7. Forcing and Jacobian Assembly

### 7.1 Kinetic (`ProcessSet`)
Kinetic contributions are unchanged mathematically for ODE rows.

New behavior:
- if row `i` is algebraic (`is_algebraic_variable_[i] == true`), kinetic forcing and Jacobian writes for that row are skipped.

This is implemented in both:
- non-vectorized dense/sparse paths,
- vectorized dense/sparse paths.

### 7.2 Constraints (`ConstraintSet`)
ConstraintSet supports two semantics:
- augmentation mode (`replace_state_rows_ = false`): add/subtract on constraint rows,
- row-replacement mode (`replace_state_rows_ = true`): set row values directly.

In row-replacement mode:
- forcing row is assigned `g(y)` (not incremented),
- Jacobian row entries are assigned `-dg/dy` with current sign convention (not incremented).

This guarantees constraint equations define the algebraic rows.

### 7.3 Vectorized layout correctness
For vectorized dense matrices (`VectorMatrix`), rows are not contiguous in memory.

ConstraintSet handles this by:
- building a contiguous per-cell concentration buffer before evaluating residual/Jacobian.

For vectorized sparse matrices, flat block arithmetic is not used for per-cell updates; cell updates use indexed access (`jacobian[cell][row][col]`) to avoid layout assumptions.

## 8. Rosenbrock Runtime Flow (`Solve`)

Per internal step:

1. Compute `initial_forcing`:
   - kinetic `AddForcingTerms`,
   - constraint `AddForcingTerms`.

2. Compute Jacobian:
   - kinetic `SubtractJacobianTerms`,
   - constraint `SubtractJacobianTerms`.

3. Form/factor matrix:
   - `A = alpha * M - J` via `AlphaMinusJacobian`.

4. Stage loop:
   - build stage RHS,
   - add `c_ij/H` coupling scaled by `M_ii`,
   - linear solve for `K_i`.

5. Candidate update:
   - `Ynew = Y + sum(m_i * K_i)`.

6. Embedded error estimate and adapt step size.

7. Accept/reject handling and retry if needed.

## 9. Why previous augmentation was insufficient
Historically, constraints were appended as extra rows (`n_species + i`) while constrained species remained ODE rows.

That architecture could evaluate constraint residuals but did not strictly enforce purely algebraic species dynamics in the persistent species state.

Full enforcement requires row replacement in the species system with `M_ii=0` for algebraic rows, which is what this architecture now does.

## 10. Validation and Tests

### 10.1 Key integration behavior
`EquilibriumIntegration.DAESolveWithConstraint` is enabled and now passes under full row-replacement mass-matrix behavior.

### 10.2 Additional coverage
Added/updated coverage includes:
- vectorized constraint forcing/Jacobian layout behavior,
- constrained normalized-error behavior,
- constraint row replacement mapping,
- six-stage DA regression parameter selection.

### 10.3 Full-suite status (current branch)
`ctest --output-on-failure` passes across the full suite.

## 11. Constraints and Limitations

1. One constraint per algebraic species row.
   - `ConstraintSet` throws if multiple constraints map to the same algebraic row.

2. Algebraic row target policy for equilibrium constraints is currently heuristic:
   - first product species.

3. Constraint type support is currently variant-limited (`EquilibriumConstraint`).

4. `constraint_size_` remains useful metadata but no longer implies augmented solve dimension in row-replacement mode.

## 12. Extension Guidance

### 12.1 Adding new constraint types
A new constraint type should provide:
- `Residual(...)`,
- `Jacobian(...)`,
- `AlgebraicSpecies()`.

It must be added to `Constraint::ConstraintVariant`.

### 12.2 Optional direct substitution mode
A future optimization mode may eliminate some algebraic variables analytically (direct substitution). That mode should remain separate from full mass-matrix enforcement and should not be mixed on the same constraint row.

## 13. File-Level Map

Core files for this architecture:
- `include/micm/solver/solver_builder.inl`
- `include/micm/solver/state.hpp`
- `include/micm/solver/state.inl`
- `include/micm/solver/rosenbrock.inl`
- `include/micm/solver/rosenbrock_temporary_variables.hpp`
- `include/micm/process/process_set.hpp`
- `include/micm/constraint/constraint.hpp`
- `include/micm/constraint/equilibrium_constraint.hpp`
- `include/micm/constraint/constraint_set.hpp`

Primary tests exercising the path:
- `test/integration/test_equilibrium.cpp`
- `test/unit/constraint/test_constraint_set.cpp`
- `test/unit/solver/test_rosenbrock.cpp`
- `test/regression/RosenbrockChapman/regression_test_solve.cpp`

