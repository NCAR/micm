# Fix DAE Algebraic Variable Error Estimation

## Problem

The Rosenbrock DAE solver's error estimator produces near-zero `Yerror` entries
for algebraic variables, regardless of tolerance settings. The embedded error
formula `Yerror = sum(e_i * K_i)` inherits this from the stage computation: for
algebraic rows the mass matrix diagonal `M_ii = 0`, which zeroes out the
inter-stage coupling terms `(c/H) * M_ii * K[j]`. The resulting `K[stage]` values
for algebraic rows are tiny, and so is `Yerror`.

**Consequence:** The solver cannot detect when a conservation-constraint "balance"
variable (e.g., SO2(g) = C_total - sum(other S species)) overshoots through zero.
Step acceptance is controlled entirely by differential variable tolerances, which
may be far too loose to prevent algebraic overshoot.

**Evidence:** Changing `atol` for an algebraic variable by 9 orders of magnitude
(1e-3 → 1e-12) produces identical solver behavior (same step count, same final
values, 0 rejections). The error entry is effectively zero regardless.

See `/home/user/git-repos/musica/.github/prompts/dae_algebraic_overshoot.md` and
`/home/user/git-repos/musica/.github/issues/rosenbrock_dae_algebraic_error.md`
for the full downstream analysis.

## Approach: Constraint Jacobian Error Propagation (Option B)

Use the **implicit function theorem** to propagate known differential variable
errors through the constraint Jacobian to obtain algebraic variable errors.

### Mathematical basis

For each algebraic constraint `g(y_d, y_a) = 0`, the implicit function theorem
gives us the sensitivity of the algebraic variable to differential variables:

    dy_a/dy_d = -(∂g/∂y_a)^{-1} · (∂g/∂y_d)

The algebraic error is therefore:

    δy_a = -(∂g/∂y_a)^{-1} · Σ_d (∂g/∂y_d · Yerror[d])

where `Yerror[d]` is the embedded error for each differential variable (from
the standard `sum(e_i * K_i)` formula, which works correctly for ODE rows).

### Why this works

- **Mathematically correct for index-1 DAEs.** The constraint Jacobian ∂g/∂y_a
  is nonsingular (that's what makes it index-1), so the inversion is well-defined.
- **Exact for linear constraints.** For `c1*x1 + c2*x2 + ... + xN = C`, the
  partials are just the coefficients: `∂g/∂x_i = c_i`. The formula gives
  `δxN = -(1/cN) · Σ(c_i · Yerror[i])`, which is the exact propagated error.
- **First-order accurate for nonlinear constraints** (e.g., equilibrium). The
  linearization is the same order of accuracy as the embedded error estimate.
- **Each constraint has exactly one algebraic variable**, so `∂g/∂y_a` is a
  scalar — "inversion" is just division. No linear system solve needed.
- **No extra residual evaluation needed.** The Jacobian is already computed and
  stored in `state.jacobian_`. We just read the existing entries.

### Why not Option A (constraint residual injection)?

The Rosenbrock stages implicitly maintain the constraint through Jacobian
coupling, so `g(Ynew)` is near-zero even when the algebraic variable has
overshot badly. The residual measures constraint satisfaction, not error.

### Why not step-change `Ynew - Y`?

This measures total change, not error. During legitimate fast transients (e.g.,
rapid aqueous uptake), the algebraic variable *should* change a lot. Using total
change as the error metric forces step rejection during physically correct
behavior.

## Implementation Plan

### Task 1: C++ Test — Prove `Yerror` Is Insensitive to Algebraic `atol` ✅

Already implemented in `test/integration/test_dae_algebraic_error_insensitivity.cpp`.
Tests `ErrorSensitiveToAlgebraicAtol` and `AlgebraicVariableDoesNotOvershootDeeply`
demonstrate the defect with a 5-species cascade system.

### Task 2: Add `ErrorPropagationFunction` to Constraint Types

Add a new method to `LinearConstraint` and `EquilibriumConstraint` that returns
a pre-compiled function to propagate differential errors to the algebraic variable.

**Signature** (same for both constraint types):

```cpp
template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
std::function<void(const DenseMatrixPolicy& Yerror, const SparseMatrixPolicy& jacobian, DenseMatrixPolicy& Yerror_out)>
ErrorPropagationFunction(
    const ConstraintInfo& info,
    const auto& state_variable_indices,
    const SparseMatrixPolicy& jacobian) const;
```

The function reads Jacobian entries for the constraint row and computes:

    Yerror[algebraic_col] = -(1 / J[row][algebraic_col]) · Σ_{d ∈ differential} J[row][d] · Yerror[d]

where:
- `row` = the algebraic variable's row in the Jacobian
- `J[row][algebraic_col]` = `∂g/∂y_a` (the diagonal entry for the algebraic var)
- `J[row][d]` = `∂g/∂y_d` for each differential dependency
- The Jacobian is stored with the solver's sign convention: `SubtractJacobianTerms`
  subtracts `∂g/∂y`, so the stored values are `-∂g/∂y`. The formula must account
  for this: stored `J_stored = -∂g/∂y`, so `∂g/∂y = -J_stored`, and the formula
  becomes `δy_a = -(1 / (-J_stored[row][a])) · Σ(-J_stored[row][d] · Yerror[d])`
  = `-(1 / J_stored[row][a]) · Σ(J_stored[row][d] · Yerror[d])`.

**Important:** The sign convention simplifies nicely. Since both numerator and
denominator use the same stored (negated) values, the negations cancel:

    Yerror[a] = -(1 / J_stored[row][a]) · Σ_d (J_stored[row][d] · Yerror[d])

**LinearConstraint implementation:**

For `c1*x1 + c2*x2 + ... + cN*xN = C` where `xN` is algebraic:
- `J_stored[row][i] = -c_i` (from `SubtractJacobianTerms`)
- `Yerror[xN] = -(1/(-cN)) · Σ_{i≠N} (-c_i · Yerror[i]) = (1/cN) · Σ(c_i · Yerror[i])`
- Since coefficients are constant, this can be pre-computed without reading the
  Jacobian at runtime. But for uniformity with equilibrium, we read from the
  Jacobian anyway (the values are there).

Actually, for linear constraints the Jacobian entries are constants, so we can
use the pre-captured coefficients directly. For equilibrium constraints, we must
read from the Jacobian because the partials are state-dependent.

**EquilibriumConstraint implementation:**

For `K_eq · ∏R^s - ∏P^s = 0` where `P[0]` is algebraic:
- Partials are concentration-dependent, already computed by `SubtractJacobianTerms`
  and stored in `state.jacobian_`
- Read them at error-computation time via sparse matrix flat IDs

**Implementation approach:** Use `SparseMatrixPolicy::Function()` to read Jacobian
entries. The function needs:
- The flat ID for `J[row][algebraic_col]` (denominator)
- The flat IDs for `J[row][d]` for each differential dependency
- Access to `Yerror` to read differential errors and write the algebraic error

Since the Jacobian is sparse and we need specific entries, pre-capture the flat
IDs during compilation (same pattern as `JacobianFunction`).

### Task 3: Add `PropagateAlgebraicErrors` to `ConstraintSet`

**File:** `include/micm/constraint/constraint_set.hpp`

Add a new method and corresponding function storage:

```cpp
// New member variables
std::vector<std::function<void(const DenseMatrixPolicy&, const SparseMatrixPolicy&, DenseMatrixPolicy&)>>
    constraint_error_propagation_functions_;
std::vector<std::function<void(const DenseMatrixPolicy&, const SparseMatrixPolicy&, DenseMatrixPolicy&)>>
    external_constraint_error_propagation_functions_;

// New method
void PropagateAlgebraicErrors(
    const DenseMatrixPolicy& Yerror_in,
    const SparseMatrixPolicy& jacobian,
    DenseMatrixPolicy& Yerror_out) const
{
  for (const auto& fn : constraint_error_propagation_functions_)
    fn(Yerror_in, jacobian, Yerror_out);
  for (const auto& fn : external_constraint_error_propagation_functions_)
    fn(Yerror_in, jacobian, Yerror_out);
}
```

Wire the compilation in `SetConstraintFunctions()` and
`SetExternalModelConstraintFunctions()`, calling each constraint's
`ErrorPropagationFunction()`.

### Task 4: Fix the DAE Mass-Matrix Coupling to Use `MatrixPolicy::Function()` ✅

Already implemented. The raw `[]` indexing in the stage coupling loop has been
replaced with `DenseMatrixPolicy::Function()`.

### Task 5: Replace Step-Change Injection with Error Propagation in `Solve()`

**File:** `include/micm/solver/rosenbrock.inl`

**Replace the current block** (after `Yerror.Axpy(...)`, before `NormalizedError()`):

```cpp
// Current (WRONG — measures total change, not error):
if (has_constraints) { ... inject_step_change ... Yerror = Ynew - Y for algebraic vars ... }
```

**With:**

```cpp
if (has_constraints)
{
  // Propagate differential variable errors through constraint Jacobian
  // to obtain algebraic variable errors via implicit function theorem:
  //   δy_a = -(∂g/∂y_a)^{-1} · Σ_d(∂g/∂y_d · Yerror[d])
  constraints_.PropagateAlgebraicErrors(Yerror, state.jacobian_, Yerror);
}
```

Note: `Yerror` is used as both input and output. The propagation functions read
differential entries (which are correct from the embedded formula) and write
algebraic entries (which were near-zero). This is safe because each constraint
writes exactly one algebraic row and reads only differential rows.

### Task 6: External Model API — New `ErrorPropagationFunction`

External models that provide constraints must implement a new method for error
propagation. This is needed so downstream applications (e.g., musica) can
participate in the algebraic error estimation.

**New concept requirement** (added to `HasConstraints` or as a separate concept):

External models already provide:
- `ConstraintResidualFunction<DMP>(param_indices, var_indices)` → residual function
- `ConstraintJacobianFunction<DMP, SMP>(param_indices, var_indices, jacobian)` → Jacobian function

They must now also provide:

```cpp
template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
std::function<void(const DenseMatrixPolicy&, const SparseMatrixPolicy&, DenseMatrixPolicy&)>
ConstraintErrorPropagationFunction(
    const std::unordered_map<std::string, std::size_t>& state_variable_indices,
    const SparseMatrixPolicy& jacobian) const;
```

**What the function must do:**

For each algebraic constraint the model provides:
1. Read `Yerror[d]` for each differential species the constraint depends on
2. Read `J[constraint_row][d]` from the Jacobian for each differential dependency
3. Read `J[constraint_row][a]` for the algebraic variable (denominator)
4. Write `Yerror[a] = -(1/J[a]) · Σ(J[d] · Yerror[d])`

**Implementation guidance for external models:**

The function should use `SparseMatrixPolicy::Function()` for Jacobian access and
`DenseMatrixPolicy::Function()` for reading/writing `Yerror`. Pre-capture flat
IDs for the needed Jacobian entries during compilation (when `state_variable_indices`
and `jacobian` sparsity are available).

**ExternalModelConstraintSet changes:**

Add to `ExternalModelConstraintSet<DMP, SMP>`:

```cpp
// New member
std::function<std::function<void(const DenseMatrixPolicy&, const SparseMatrixPolicy&, DenseMatrixPolicy&)>(
    const std::unordered_map<std::string, std::size_t>& state_variable_indices,
    const SparseMatrixPolicy& jacobian)>
    get_error_propagation_function_;
```

Wire in the constructor:

```cpp
get_error_propagation_function_ = [shared_model](
    const std::unordered_map<std::string, std::size_t>& state_variable_indices,
    const SparseMatrixPolicy& jacobian)
    -> std::function<void(const DenseMatrixPolicy&, const SparseMatrixPolicy&, DenseMatrixPolicy&)>
{
  return shared_model->template ConstraintErrorPropagationFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
      state_variable_indices, jacobian);
};
```

Wire compilation in `SetExternalModelConstraintFunctions()`:

```cpp
external_constraint_error_propagation_functions_.push_back(
    model.get_error_propagation_function_(state_variable_indices, jacobian));
```

### Task 7: Update Tests

- **`test_dae_algebraic_error_insensitivity.cpp`**: After the fix, the
  `ErrorSensitiveToAlgebraicAtol` test should show different step counts for
  different `atol` values. Update assertions accordingly.
- **`test_dae_constraint_overshoot.cpp`**: Verify the algebraic variable stays
  non-negative with practical tolerances.
- **All existing Rosenbrock tests**: Must continue to pass (regression check).

### Task 8: Test Equilibrium Constraint Scaling

The equilibrium constraint `K_eq · ∏R^s - ∏P^s = 0` has state-dependent
partials. Verify that:
1. The Jacobian entries are read correctly at error-estimation time
2. The error propagation produces reasonable values (same order as differential
   errors, not wildly different due to scaling)
3. The `test_dae_constraint_overshoot.cpp:EquilibriumPlusConservation` test
   (if it exists) works correctly

## Implementation Order

1. **Task 2** — Add `ErrorPropagationFunction` to `LinearConstraint` and
   `EquilibriumConstraint`.
2. **Task 3** — Add `PropagateAlgebraicErrors` to `ConstraintSet`, wire compilation.
3. **Task 6** — Update `ExternalModelConstraintSet` and external model API.
4. **Task 5** — Replace step-change injection in `Solve()` with error propagation.
5. **Task 7** — Update tests and verify.
6. **Task 8** — Equilibrium scaling verification.

## Files Modified

- `include/micm/constraint/types/linear_constraint.hpp` — Task 2 (new method)
- `include/micm/constraint/types/equilibrium_constraint.hpp` — Task 2 (new method)
- `include/micm/constraint/constraint_set.hpp` — Task 3 (new storage + method)
- `include/micm/external_model.hpp` — Task 6 (new wrapper member + wiring)
- `include/micm/solver/rosenbrock.inl` — Task 5 (replace step-change block)
- `test/integration/test_dae_algebraic_error_insensitivity.cpp` — Task 7
- `test/integration/test_dae_constraint_overshoot.cpp` — Task 7 (verify)

## Coding Constraints

- **Never use square bracket `[]` indexing** on matrices in solver code. All
  matrix iteration must use `DenseMatrixPolicy::Function()` / `SparseMatrixPolicy::Function()`.
- **Never assume matrix ordering.** Use `ForEachRow`/`GetColumnView`/`ForEachBlock`/`GetBlockView`.
- Follow existing patterns in `rosenbrock.inl` and constraint types.
- Test code (gtest) may use `[]` indexing — this constraint applies only to
  the solver implementation.

## Risks and Open Questions

1. **Jacobian freshness**: The Jacobian in `state.jacobian_` is computed once per
   step (at the beginning of the step acceptance loop) and used for all stages.
   Using it for error propagation at the end of the step is an approximation —
   the Jacobian may have changed during the step. This is acceptable because the
   embedded error estimate itself is only first-order accurate, and the Jacobian
   change over one step is O(H²).

2. **Sign convention verification**: The `SubtractJacobianTerms` convention stores
   `-∂g/∂y` in the Jacobian. Must verify this carefully for both linear and
   equilibrium constraints. The negations cancel in the ratio, but getting the
   sign wrong would flip the error direction.

3. **Division by zero**: If `J[row][algebraic_col]` is zero, the constraint is
   degenerate (not index-1). Add a guard: if `|J[a]| < epsilon`, skip the error
   propagation for that constraint (fall back to near-zero error, which is the
   current behavior). Log a warning if feasible.

4. **Performance**: Reading sparse Jacobian entries is fast (pre-captured flat IDs).
   No matrix solves or extra evaluations needed. Should be cheaper than Option A's
   extra `AddForcingTerms()` call.
