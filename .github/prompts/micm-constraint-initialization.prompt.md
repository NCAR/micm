# Feature Request: DAE Constraint Initialization for the Rosenbrock Solver

## Summary

Add automatic constraint initialization to MICM's Rosenbrock DAE solver.
When the solver detects algebraic constraints, it should Newton-iterate
the algebraic variables to satisfy G(y) = 0 before the first time step.
This prevents a class of silent failures where the solver reports
`SolverState::Converged` but produces wildly wrong results because the
initial conditions were not on the constraint manifold.

## Motivation

MICM's Rosenbrock DAE solver uses a mass-matrix formulation where
algebraic variables (M_ii = 0) are determined by constraint equations
G(y) = 0 rather than by time integration. The error estimator in
`NormalizedError` correctly excludes algebraic variables from the error
norm — but this means there is **no mechanism to detect or correct
inconsistent initial conditions** for algebraic variables.

In practice, when a user provides initial conditions where G(y₀) ≠ 0,
the solver's first step attempts a large Newton-like correction that
overshoots catastrophically. The error estimator sees no problem (it
only checks ODE variables), reports convergence, and returns garbage.

This was discovered during MIAM's CAM cloud chemistry integration test:
a system with 5 algebraic constraints and 1 differential variable
(water, constant). The Kw equilibrium product was off by 10,000× in the
initial conditions. The solver reported `Converged` but concentrations
jumped by 25 orders of magnitude in one step.

Other DAE solvers (SUNDIALS IDA, DASSL) include constraint initialization
routines for exactly this reason. MICM already has all the infrastructure
needed — it just needs to use it before the first step.

## Design

### Approach: Automatic initialization inside `Solve()`

Before the main time-stepping loop in `AbstractRosenbrockSolver::Solve()`
(file: `include/micm/solver/rosenbrock.inl`), add a constraint
initialization phase that runs **once** on the first call (or whenever
requested). This Newton-iterates the algebraic variables to satisfy
G(y_alg) = 0, holding differential variables fixed.

### Algorithm

```
function InitializeConstraints(state, constraints, linear_solver):
    if constraints.Size() == 0:
        return success

    max_iterations = 10
    tolerance = 1e-10

    for iter = 1 to max_iterations:
        // 1. Evaluate constraint residuals
        residual.Fill(0)
        constraints.AddForcingTerms(state.variables_, state.custom_rate_parameters_, residual)

        // 2. Check convergence: ||G||_∞ over algebraic rows only
        max_residual = 0
        for each algebraic variable i (where upper_left_identity_diagonal_[i] == 0):
            max_residual = max(max_residual, |residual[i]|)
        if max_residual < tolerance:
            return success

        // 3. Compute constraint Jacobian
        state.jacobian_.Fill(0)
        constraints.SubtractJacobianTerms(state.variables_, state.custom_rate_parameters_, state.jacobian_)

        // 4. Form the system matrix for algebraic rows
        //    For algebraic rows: W = -J (since M_ii = 0, AlphaMinusJacobian gives -J)
        //    For ODE rows: we need identity to preserve them (no update)
        //    Use alpha = 0 in AlphaMinusJacobian, then manually set ODE diagonals to 1
        //    OR: just negate the Jacobian and set ODE rows to identity manually
        //
        //    Simpler: use the existing AlphaMinusJacobian with alpha=0.
        //    This gives (0*M - J) = -J for all rows.
        //    Then set residual[i] = 0 for all ODE variables (don't update them).
        AlphaMinusJacobian(state, 0.0)  // forms -J

        // 5. Factor and solve: -J * delta = residual  →  delta = -J^{-1} * G
        linear_solver.Factor(...state.jacobian_...)
        // Copy residual to a work vector, zero out ODE rows
        for each i where upper_left_identity_diagonal_[i] > 0:
            residual[i] = 0  // don't update differential variables
        linear_solver.Solve(residual, ...state.jacobian_...)

        // 6. Update algebraic variables: y_new = y + delta
        state.variables_ += residual  // residual now contains delta

    return failure  // did not converge in max_iterations
```

**Important subtlety**: The approach above of zeroing out residuals for
ODE rows and using the full Jacobian is simple but may not work cleanly
because the Jacobian may have cross-coupling between ODE and algebraic
rows. A cleaner approach: always compute the full Newton step but only
**apply** the update to algebraic variables:

```
// After linear_solver.Solve(delta, ...):
for each variable i:
    if upper_left_identity_diagonal_[i] == 0:  // algebraic
        state.variables_[cell][i] += delta[cell][i]
    // else: don't touch differential variables
```

This is simpler, robust, and matches how SUNDIALS IDA does it.

### Where to add this

**Primary change**: `include/micm/solver/rosenbrock.inl` in the `Solve()` method.

Add the initialization phase after the temporary variable setup (line ~22)
and before the main `while (present_time - time_step + round_off <= 0)`
loop (line ~32). The method should:

1. Check `if (constraints_.Size() > 0)` — skip entirely for pure ODE systems
2. Run the Newton iteration described above
3. If initialization fails after `max_iterations`, return a new
   `SolverState` indicating the failure (don't silently proceed)

**Tracking first-call**: Add a `bool needs_constraint_initialization_`
flag. Two options:
- (a) Set it to `true` in the constructor and `false` after first
  successful initialization. This means initialization runs once
  automatically.
- (b) Make `Solve()` always check `G(y₀)` at the start and only run
  Newton if the residual exceeds tolerance. This is stateless and
  handles the case where the user modifies state between calls (e.g.,
  after a discontinuous event like cloud activation). **This option
  is preferred.**

### API Changes

1. **New `SolverState` enum value**:
   ```cpp
   // In include/micm/solver/solver_result.hpp
   enum class SolverState
   {
     // ... existing values ...
     ConstraintInitializationFailed,  // NEW
   };
   ```
   Update `SolverStateToString()` accordingly.

2. **New parameters** (optional, with good defaults):
   ```cpp
   // In RosenbrockSolverParameters
   std::size_t constraint_init_max_iterations_ = 10;
   double constraint_init_tolerance_ = 1e-10;
   ```

3. **New stats fields** (optional):
   ```cpp
   // In SolverStats
   uint64_t constraint_init_iterations_{};
   ```

4. **No changes to the public `Solve()` signature.** The initialization
   is transparent to callers.

## Existing Infrastructure to Use

All of these already exist and should be reused — no new math code is
needed:

| What | Where | How to use |
|---|---|---|
| Algebraic variable IDs | `state.upper_left_identity_diagonal_` | `== 0.0` means algebraic |
| Constraint residuals | `constraints_.AddForcingTerms(Y, params, forcing)` | Fills forcing with G(y) for algebraic rows |
| Constraint Jacobian | `constraints_.SubtractJacobianTerms(Y, params, jac)` | Fills jac with -∂G/∂y |
| System matrix formation | `AlphaMinusJacobian(state, alpha)` | With alpha=0 gives -J |
| LU factorization | `linear_solver_.Factor(...)` | Already used in stage loop |
| Linear solve | `linear_solver_.Solve(delta, ...)` | Already used in stage loop |
| Temporary storage | `state.temporary_variables_` → `initial_forcing_` | Can reuse as residual/delta workspace |

The `ConstraintSet::AddForcingTerms` method computes constraint
residuals for algebraic rows using **direct assignment** (replaces
the row), so Fill(0) before calling it is sufficient.

## Codebase Reference

These are the files you will need to read and understand:

### Files to modify

- `include/micm/solver/rosenbrock.hpp` — Add initialization method declaration
- `include/micm/solver/rosenbrock.inl` — Implement initialization in Solve() flow
- `include/micm/solver/solver_result.hpp` — Add `ConstraintInitializationFailed` enum value, update `SolverStateToString()`
- `include/micm/solver/rosenbrock_solver_parameters.hpp` — Add `constraint_init_max_iterations_` and `constraint_init_tolerance_` parameters (optional)

### Files for reference (read but don't modify)

- `include/micm/constraint/constraint_set.hpp` — `AddForcingTerms()`, `SubtractJacobianTerms()`, `Size()`
- `include/micm/solver/solver_builder.inl` — How algebraic variables are registered, mass matrix diagonal set
- `include/micm/solver/state.hpp` — `upper_left_identity_diagonal_`, `variables_`, `jacobian_`, `temporary_variables_`

### Key types and signatures

```cpp
// The solver class (CRTP pattern)
template<class RatesPolicy, class LinearSolverPolicy, class ConstraintSetPolicy, class Derived>
class AbstractRosenbrockSolver

// Solve method
SolverResult Solve(double time_step, auto& state,
                   const RosenbrockSolverParameters& parameters) const noexcept;

// Constraint residuals: fills forcing[algebraic_row] = G(y)
constraints_.AddForcingTerms(
    const DenseMatrixPolicy& state_variables,
    const DenseMatrixPolicy& state_parameters,
    DenseMatrixPolicy& forcing) const;

// Constraint Jacobian: fills jacobian entries with -dG/dy
constraints_.SubtractJacobianTerms(
    const DenseMatrixPolicy& state_variables,
    const DenseMatrixPolicy& state_parameters,
    SparseMatrixPolicy& jacobian) const;

// Mass matrix diagonal: 1.0 = ODE, 0.0 = algebraic
state.upper_left_identity_diagonal_  // std::vector<double>

// System matrix: forms alpha*M - J in state.jacobian_
AlphaMinusJacobian(auto& state, const double& alpha) const;

// Linear solver (two variants, dispatched at compile time):
// Separate L/U:
linear_solver_.Factor(state.jacobian_, state.lower_matrix_, state.upper_matrix_);
linear_solver_.Solve(delta, state.lower_matrix_, state.upper_matrix_);
// In-place:
linear_solver_.Factor(state.jacobian_);
linear_solver_.Solve(delta, state.jacobian_);
```

## Testing

### 1. Unit test: `test_constraint_initialization.cpp`

Test the initialization in isolation with a simple 3-variable system:

```
Variables: A (differential), B (algebraic), C (algebraic)
Constraint 1: B = 2*A  (linear constraint, B algebraic)
Constraint 2: C = A*A  (quadratic constraint, C algebraic)

Test cases:
a) Consistent ICs: A=1, B=2, C=1 → no initialization needed, verify
   state unchanged
b) Inconsistent ICs: A=1, B=5, C=10 → after initialization, B=2, C=1,
   A unchanged
c) Mildly inconsistent: A=1, B=2.001, C=1.001 → converges in 1 iteration
d) Severely inconsistent: A=1, B=1000, C=-500 → converges in a few iterations
e) Pure ODE system (no constraints): verify initialization is skipped
```

### 2. Unit test: `test_constraint_initialization_failure.cpp`

Test that initialization failure is properly reported:

```
a) Impossible constraints (e.g., B = 1 AND B = 2 simultaneously) →
   returns ConstraintInitializationFailed
b) Singular constraint Jacobian → returns appropriately
   (RepeatedlySingularMatrix or ConstraintInitializationFailed)
```

### 3. Integration test: Inconsistent initial conditions recovery

Reproduce the exact failure scenario from MIAM's CAM cloud chemistry:

```
System: A(g) ⇌ A(aq) with Henry's Law + mass conservation
- Set A_g = 1e-8, A_aq = 0 (inconsistent with HLC equilibrium)
- Solve for 1 second
- Verify: mass conservation holds, HLC equilibrium reached,
  result matches a run started from consistent ICs
```

### 4. Integration test: Verify existing tests still pass

Run all existing Rosenbrock tests to verify no regressions:
- `test_equilibrium_constraint`
- `test_solver_configuration`
- All aerosol model tests
- All tutorial tests

The initialization should be a no-op (or at most 1 iteration) for
systems that already have consistent ICs.

### 5. Integration test: Multi-cell systems

Verify initialization works correctly with multiple grid cells
(DenseMatrix with NumRows() > 1), where each cell may have different
initial conditions and different levels of inconsistency.

### 6. Performance test (optional)

Verify that for systems with consistent ICs, the overhead of checking
the residual at the start of Solve() is negligible (one constraint
evaluation + one max-norm check, no Jacobian/factor/solve).

## Edge Cases to Handle

1. **No algebraic variables**: `constraints_.Size() == 0` → skip entirely
2. **All variables algebraic, no ODE variables**: Should still work —
   the Newton solve is the entire system
3. **Multiple grid cells**: Each cell is independent; iterate per-cell
   or check convergence per-cell
4. **Constraint Jacobian is singular**: Return
   `ConstraintInitializationFailed` (don't hang or crash)
5. **NaN/Inf during initialization**: Detect and return appropriate
   `SolverState`
6. **Subsequent Solve() calls**: If the user calls Solve() multiple
   times on the same state, the second call should detect that ICs are
   already consistent and skip the Newton iteration (just one residual
   evaluation, no Jacobian/factor/solve)
7. **User modifies state between Solve() calls**: The residual check at
   entry handles this naturally — if the state was perturbed off-manifold,
   re-initialization triggers automatically
