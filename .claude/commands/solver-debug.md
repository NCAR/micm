# Solver Debug Agent

Investigate numerical issues in MICM solvers.

## Usage

- `/solver-debug` - Start numerical debugging investigation

## Instructions

When the user invokes this skill, launch an investigation of solver numerical issues.

### Phase 1: Gather Context

1. **Ask the user** what symptom they're seeing:
   - Values exploding (which variable?)
   - NaN or Inf appearing
   - Wrong steady state
   - Solver not converging
   - Step size going to minimum

2. **Read the relevant solver file**:
   - For Rosenbrock: `include/micm/solver/rosenbrock.inl`
   - For Backward Euler: `include/micm/solver/backward_euler.inl`

### Phase 2: Trace Data Flow

For Rosenbrock solver issues, trace through:

1. **Forcing calculation** (`AddForcingTerms`):
   - ODE forcing: `rates_.AddForcingTerms()`
   - Constraint forcing: `constraints_.AddForcingTerms()`
   - Expected magnitude: O(rate × concentration)

2. **Jacobian calculation** (`SubtractJacobianTerms`):
   - Stores `-J_true` (negative of true Jacobian)
   - ODE Jacobian: `rates_.SubtractJacobianTerms()`
   - Constraint Jacobian: `constraints_.SubtractJacobianTerms()`

3. **Matrix formation** (`AlphaMinusJacobian`):
   - Adds `alpha = 1/(H * gamma)` to ALL diagonals
   - Result should have all diagonals of similar magnitude
   - If constraint diagonals are much smaller → ill-conditioning

4. **Linear solve**:
   - Solves `(αI - J) * K = forcing`
   - K values should scale with H

5. **Stage computation**:
   - Uses `c/H * K` terms
   - If K doesn't scale with H → c/H amplifies errors

### Phase 3: Identify Root Cause

Common patterns:

| Symptom | Root Cause | Solution |
|---------|------------|----------|
| Constraint var explodes | Missing alpha on constraint diagonal | Check AlphaMinusJacobian |
| All vars explode | Singular Jacobian | Check for zero rate constants |
| Wrong equilibrium | Sign error | Check SubtractJacobianTerms |
| Slow convergence | Step too large | Reduce h_max |

### Phase 4: Recommend Fix

Based on the root cause:

1. **Code fix**: If a bug is found, show the fix
2. **Parameter tuning**: Suggest solver parameter changes
3. **Usage fix**: Suggest projection steps or smaller time steps
4. **Test addition**: Suggest a test case to prevent regression

## Key Files Reference

- `rosenbrock.inl:52-60` - Forcing calculation
- `rosenbrock.inl:62-71` - Jacobian calculation
- `rosenbrock.inl:224-261` - AlphaMinusJacobian (both versions)
- `constraint_set.hpp` - Constraint forcing/Jacobian methods
