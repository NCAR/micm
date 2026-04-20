# DAE Algebraic Variable Tolerance Guide

## What changed

MICM's Rosenbrock DAE solver now uses **step-change error estimation** for
algebraic variables. Previously, the embedded error formula (`Yerror = Σ e_i·K_i`)
produced near-zero entries for algebraic rows because the mass matrix diagonal
`M_ii = 0` zeroes out the inter-stage coupling. The solver was completely
insensitive to algebraic tolerances — changing `atol` by 9 orders of magnitude
had no effect on step acceptance.

Now, after computing the embedded error, the solver replaces each algebraic
entry with the step change:

```
Yerror[a] = Ynew[a] - Y[a]
```

This means **algebraic variable tolerances now matter** and must be set
appropriately.

## How this affects downstream applications

Any code that uses `LinearConstraint` or `EquilibriumConstraint` (or external
model constraints) is affected. The solver will now reject steps where algebraic
variables change by more than their tolerance permits, which:

1. **Prevents overshoot** — balance variables can no longer go deeply negative
2. **May increase step counts** — if algebraic tolerances are set too tight
3. **May cause convergence failure** — if algebraic tolerances are unreasonably
   small relative to the legitimate per-step changes

## Required tolerance updates

### Key principle

Algebraic variables are **not integrated** — they are determined by constraint
equations at each step. Their `atol` controls how much the variable is allowed
to change per internal step before the solver rejects and retries with a smaller
step size.

**Rule of thumb:** Set `atol` for algebraic variables to the smallest physically
meaningful value, not to a numerical precision target. A value of `1e-8` to
`1e-10` is typical for atmospheric chemistry concentrations in mol/m³.

### Setting per-species tolerances

Instead of a uniform tolerance for all species:

```cpp
// OLD — same tight tolerance for everything (may cause convergence failure)
state.SetAbsoluteTolerances(std::vector<double>(N, 1.0e-12));
```

Set per-species tolerances using the variable map:

```cpp
auto state = solver.GetState(1);
state.SetRelativeTolerance(1.0e-6);

// Start with a tight default for differential variables
std::vector<double> atols(num_species, 1.0e-12);

// Use moderate tolerance for algebraic variables
// (equilibrium species, conservation balance species)
atols[state.variable_map_.at("A_gas")] = 1.0e-8;
atols[state.variable_map_.at("A_aq")]  = 1.0e-8;

state.SetAbsoluteTolerances(atols);
```

### Tolerance sizing guidance

| Variable type | Recommended `atol` | Rationale |
|---|---|---|
| Differential (ODE) species | `1e-12` to `1e-14` | Standard numerical precision |
| Equilibrium algebraic species | `1e-8` to `1e-10` | These change by O(concentration) per step during transients |
| Conservation balance species | `1e-8` to `1e-10` | Same reasoning — tracks the residual budget |

If the solver returns `SolverState::ConvergenceExceededMaxSteps` after this
change, the algebraic `atol` is likely too tight. Loosen it or increase
`max_number_of_steps_`.

## Testing your constraints

### 1. Verify convergence

After updating tolerances, confirm the solver still converges:

```cpp
auto result = solver.Solve(dt, state);
EXPECT_EQ(result.state_, SolverState::Converged);
```

### 2. Check conservation

For every `LinearConstraint` with constant `C_total`, verify the sum is
preserved after each solve:

```cpp
double sum = 0;
for (const auto& [species, coeff] : constraint_terms)
  sum += coeff * state.variables_[0][state.variable_map_.at(species.name_)];
EXPECT_NEAR(sum, C_total, 1.0e-12);
```

### 3. Check non-negativity of balance variables

The balance variable should never go deeply negative. A small negative value
(within `atol`) is acceptable:

```cpp
double balance_val = state.variables_[0][balance_idx];
EXPECT_GE(balance_val, -atol_for_balance)
    << "Balance variable went negative: " << balance_val;
```

### 4. Sensitivity test (optional)

Verify the solver responds to algebraic tolerance changes:

```cpp
auto r_loose = RunSystem(balance_atol = 1e-3);
auto r_tight = RunSystem(balance_atol = 1e-8);

// Tight tolerance should require more internal steps
EXPECT_GT(r_tight.accepted_steps, r_loose.accepted_steps);
```

## Reference tests

See these MICM integration tests for complete working examples:

- `test/integration/test_dae_constraint_overshoot.cpp` — Simple conservation
  (3 species) and equilibrium + conservation (3 species with 2 constraints).
  Demonstrates per-species tolerance setting and non-negativity checks.

- `test/integration/test_dae_algebraic_error_insensitivity.cpp` — Cascade
  equilibrium system (4 species, 2 equilibria, 1 conservation). Tests tolerance
  sensitivity and overshoot prevention with the step-change error estimate.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ConvergenceExceededMaxSteps` after update | Algebraic `atol` too tight for the step-change error estimate | Loosen algebraic `atol` to `1e-8` or increase `max_number_of_steps_` |
| Balance variable goes slightly negative (within `atol`) | Expected — the step-change error is O(H), not exact | Tighten `atol` if this is unacceptable |
| Many more steps than before | Step-change error is conservative (O(H) vs O(H^(p+1))) | Expected for algebraic variables; loosen `atol` if too slow |
| No change in behavior | Variable may not actually be algebraic | Check constraint ordering — the algebraic species must be last (linear) or first product (equilibrium) |
