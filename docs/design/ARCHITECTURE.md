# MICM Rosenbrock DAE Architecture

## 1. Purpose
This document captures the current architecture for MICM's CPU Rosenbrock solver with algebraic constraints, including recent fixes for full mass-matrix row replacement and Jacobian sparsity pruning robustness.

Scope:
- kinetic ODE + algebraic constraint model,
- species-row replacement DAE formulation,
- builder-time sparsity/mass-matrix wiring,
- runtime forcing/Jacobian assembly,
- coverage and known limitations.

## 2. Mathematical Form

### 2.1 Kinetics
For species state vector `y`:

`dy/dt = f_chem(y)`

### 2.2 Constraints
For `m` constraints:

`g_k(y) = 0`, `k=1..m`

Implemented constraint type today:
- `EquilibriumConstraint`:
  `G = K_eq * product(reactants^stoich) - product(products^stoich)`

### 2.3 DAE formulation
MICM enforces constraints with a diagonal mass matrix:

`M * dy/dt = F(y)`

with:
- `M_ii = 1` for ODE rows,
- `M_ii = 0` for algebraic rows,
- `F_i = f_chem,i` on ODE rows,
- `F_i = g_i` on algebraic rows.

Important: constrained species remain part of the species state vector; there is no extra appended unknown block.

## 3. State and Dimensions

- `state_size_ = n_species`
- `constraint_size_ = n_constraints` (metadata only)
- `variables_` shape: `(n_cells, n_species)`
- Jacobian shape: `n_species x n_species`
- Rosenbrock temporary matrices are species-sized (not species+constraints sized)

Mass diagonal is stored in `state.upper_left_identity_diagonal_` and derived from `StateParameters.mass_matrix_diagonal_`.

## 4. Build-Time Wiring (`SolverBuilder::Build`)

Current build sequence:

1. Build kinetic `RatesPolicy` (`ProcessSet`) and kinetic Jacobian sparsity.
2. Build `ConstraintSet` (row-replacement mode only) when constraints exist.
3. Collect `algebraic_variable_ids` from constraints.
4. Pass algebraic IDs to `ProcessSet` so kinetic forcing/Jacobian skips algebraic rows.
5. Prune kinetic sparsity entries whose dependent row is algebraic.
6. Merge constraint Jacobian sparsity into the global sparsity set.
7. Build species-sized sparse Jacobian and linear-solver structures.
8. Set Jacobian flat IDs for rates and constraints.
9. Build mass-matrix diagonal (`1` ODE rows, `0` algebraic rows).

Key robustness detail (recent fix):
- After step 5, some algebraic-row kinetic entries no longer exist in the sparse structure.
- `ProcessSet::SetJacobianFlatIds` now skips `VectorIndex(...)` lookups for algebraic dependent rows and stores placeholder IDs to keep update-loop indexing aligned.
- This avoids `Zero element access` throws for valid constrained systems.

## 5. Forcing and Jacobian Assembly

### 5.1 Kinetic contributions (`ProcessSet`)
- Forcing and Jacobian are computed as before for ODE rows.
- Writes into algebraic rows are skipped in both scalar and vectorized paths.

### 5.2 Constraint contributions (`ConstraintSet`)
ConstraintSet is row-replacement only:
- forcing row is assigned to constraint residual `g(y)`,
- Jacobian row contributions are applied with solver sign convention (`-dg/dy`),
- duplicate dependency columns accumulate correctly via subtraction on pre-zeroed rows.

Constraint row mapping:
- each constraint targets one algebraic species row (`GetAlgebraicSpecies()`).
- duplicate mapping to the same row is rejected.

## 6. Rosenbrock Runtime

Per internal step:

1. Build forcing:
   - `ProcessSet::AddForcingTerms`
   - `ConstraintSet::AddForcingTerms`
2. Build Jacobian:
   - `ProcessSet::SubtractJacobianTerms`
   - `ConstraintSet::SubtractJacobianTerms`
3. Form/factor:
   - `A = alpha * M - J`
4. Stage loop:
   - includes `c_ij / H` coupling scaled row-wise by `M_ii`
   - algebraic rows (`M_ii = 0`) do not receive the stiff `1/H` coupling term
5. Candidate update + embedded error + adaptive step control.

## 7. Post-Solve Clamping

Both `Solver::Solve` overloads share `PostSolveClamp`:
- DAE mode (`constraints_replace_state_rows_ == true`): clamp only ODE rows (`M_ii > 0`)
- non-DAE mode: clamp all rows

This keeps overload behavior consistent and avoids clamping algebraic rows.

## 8. Reordering and Constraints

State reordering (`SetReorderState(true)`) is supported with constraints:
- species map can be reordered by sparsity heuristics,
- algebraic row IDs are resolved using reordered indices,
- mass-matrix diagonal and row replacement follow reordered mapping.

## 9. Coverage and Validation

Current branch validation:
- `ctest --output-on-failure` passes full suite.

Notable coverage for DAE architecture:
- `EquilibriumIntegration.DAESolveWithConstraint`
- `EquilibriumIntegration.DAESolveWithConstraintAndReorderState`
- `SolverBuilder.DAEKineticAlgebraicRowPruningDoesNotThrow`
- `SolverBuilder.DAEKineticAlgebraicReactantRowPruningDoesNotThrow`
- `EquilibriumIntegration.DAEStateCopyAndSolve`
- `EquilibriumIntegration.DAEOverlappingSpeciesJacobian`
- `EquilibriumIntegration.DAEConstraintOnlySpeciesNotUnused`

## 10. Known Limitations

1. Constraint variant support is currently limited to `EquilibriumConstraint`.
2. Equilibrium constants are static in the current model (no built-in temperature dependence).
3. The implementation assumes index-1 style algebraic constraints and does not perform automatic index reduction.
4. Algebraic-row target policy for equilibrium constraints is currently based on constraint metadata (`GetAlgebraicSpecies()` mapping; for equilibrium constraints this is the first product species).

## 11. Extension Guidance

To add a new constraint type:

1. Implement:
   - `Residual(...)`
   - `Jacobian(...)`
   - `AlgebraicSpecies()`
   - dependency/species-name accessors
2. Add it to `Constraint::ConstraintVariant`.
3. Ensure dependency indexing and Jacobian sparsity rows are represented in `ConstraintSet`.
4. Add integration coverage for:
   - mixed kinetic/constraint coupling,
   - reordered-state build/solve,
   - vectorized and scalar paths where applicable.

## 12. Core File Map

- `include/micm/solver/solver_builder.inl`
- `include/micm/solver/solver.hpp`
- `include/micm/solver/rosenbrock.inl`
- `include/micm/solver/state.hpp`
- `include/micm/solver/state.inl`
- `include/micm/solver/rosenbrock_temporary_variables.hpp`
- `include/micm/process/process_set.hpp`
- `include/micm/constraint/constraint.hpp`
- `include/micm/constraint/constraint_set.hpp`
- `include/micm/constraint/equilibrium_constraint.hpp`
