# MICM DAE Constraint System - Test Cases

This document describes complex test cases for validating the DAE constraint system. These tests go beyond the basic API tests to exercise realistic atmospheric chemistry scenarios.

---

## Atmospheric Chemistry Equilibria

### 1. Carbonate System (Aqueous CO2)

**Description**: The carbonate buffer system in aqueous solution, fundamental to ocean and cloud chemistry.

**Reactions**:
```
CO2(aq) + H2O <-> H2CO3 <-> HCO3- + H+ <-> CO3-- + 2H+
```

**Constraints**:
- `K1 = [H+][HCO3-] / [CO2(aq)] ≈ 4.3e-7`
- `K2 = [H+][CO3--] / [HCO3-] ≈ 4.7e-11`
- Conservation: `[CO2] + [HCO3-] + [CO3--] = C_total`
- Charge balance: `[H+] = [HCO3-] + 2[CO3--] + [OH-]`

**Initial Conditions**:
- C_total = 2.0e-3 M (typical ocean)
- pH = 8.1 (initial guess)

**Expected Steady State**:
- [CO2(aq)] ≈ 10 μM
- [HCO3-] ≈ 1.8 mM
- [CO3--] ≈ 0.2 mM

**Challenge**: Multiple coupled equilibria, wide range of K values (4 orders of magnitude), pH-dependent speciation

**Tests**:
- [ ] Correct equilibrium ratios at steady state
- [ ] Conservation of total carbon
- [ ] Charge balance maintained
- [ ] Response to CO2 addition (ocean acidification scenario)

---

### 2. Ammonia-Ammonium Equilibrium

**Description**: Gas-aqueous partitioning with acid-base equilibrium, important for atmospheric aerosol formation.

**Reactions**:
```
NH3(g) <-> NH3(aq)          Henry's law: H = 60 M/atm
NH3(aq) + H+ <-> NH4+       K = 1.7e9
```

**With sulfate** (optional extension):
```
2NH4+ + SO4-- <-> (NH4)2SO4(s)
```

**Constraints**:
- Henry equilibrium: `[NH3(aq)] = H * p_NH3`
- Acid-base: `K = [NH4+] / ([NH3(aq)] * [H+])`
- N conservation: `[NH3(g)] + [NH3(aq)] + [NH4+] = N_total`

**Initial Conditions**:
- N_total = 10 ppb (gas phase equivalent)
- pH = 4.0 (acidic aerosol)

**Expected Behavior**:
- At low pH: mostly NH4+
- At high pH: mostly NH3(aq) and NH3(g)

**Challenge**: Gas-aqueous partitioning coupled with acid-base equilibrium, pH sensitivity

**Tests**:
- [ ] Correct partitioning at different pH values
- [ ] Mass balance across phases
- [ ] Response to pH change (titration curve)

---

### 3. Sulfur Dioxide System

**Description**: SO2 dissolution and speciation, critical for acid rain and cloud chemistry.

**Reactions**:
```
SO2(g) <-> SO2(aq)           H = 1.2 M/atm
SO2(aq) + H2O <-> HSO3- + H+  K1 = 1.3e-2
HSO3- <-> SO3-- + H+          K2 = 6.3e-8
```

**Constraints**:
- Henry: `[SO2(aq)] = H * p_SO2`
- `K1 = [HSO3-][H+] / [SO2(aq)]`
- `K2 = [SO3--][H+] / [HSO3-]`
- S(IV) conservation: `[SO2] + [HSO3-] + [SO3--] = S_total`

**Initial Conditions**:
- p_SO2 = 1 ppb
- pH = 5.0 (cloud water)

**Expected Speciation** (pH 5):
- SO2(aq): ~2%
- HSO3-: ~97%
- SO3--: ~1%

**Challenge**: pH-dependent speciation, moderately stiff (K1/K2 ratio ~10^5)

**Tests**:
- [ ] Correct speciation vs pH curve
- [ ] S(IV) conservation
- [ ] Oxidation to sulfate (add slow reaction)

---

## Multi-Constraint Coupled Systems

### 4. Iron-Oxalate Photochemistry

**Description**: Iron redox cycling catalyzed by organic ligands, important in cloud and fog chemistry.

**Reactions**:
```
Fe(III) + hv -> Fe(II) + products     j = 0.01 s-1
Fe(II) + H2O2 -> Fe(III) + OH + OH-   k = 76 M-1s-1
Fe(III) + Ox <-> Fe(III)-Ox           K = 1e8 M-1
```

**Constraints**:
- Fe conservation: `[Fe(II)] + [Fe(III)] + [Fe-Ox] = Fe_total`
- Oxalate conservation: `[Ox] + [Fe-Ox] = Ox_total`
- Complexation equilibrium: `K = [Fe-Ox] / ([Fe(III)] * [Ox])`

**Initial Conditions**:
- Fe_total = 1 μM
- Ox_total = 10 μM
- [H2O2] = 10 μM

**Challenge**: Conservation constraints + equilibrium constraint together, photolysis-driven redox cycling

**Tests**:
- [ ] Fe and Ox conservation maintained
- [ ] Equilibrium satisfied at each timestep
- [ ] Correct Fe(II)/Fe(III) ratio under illumination

---

### 5. Ozone-NOx Photostationary State

**Description**: The fundamental daytime O3-NOx relationship in the troposphere.

**Reactions**:
```
NO2 + hv -> NO + O       j = 0.01 s-1
O + O2 + M -> O3         k2 = 6e-34 cm6 s-1
NO + O3 -> NO2 + O2      k3 = 1.8e-14 cm3 s-1
```

**Photostationary State Constraint**:
```
[O3] = j(NO2) * [NO2] / (k3 * [NO])
```
Or equivalently: `j*[NO2] - k3*[NO][O3] = 0`

**Initial Conditions**:
- [NOx] = [NO] + [NO2] = 10 ppb
- [O3] = 50 ppb (background)

**Expected Behavior**:
- Leighton ratio: `[O3][NO]/[NO2] = j/k3`
- Fast equilibration (seconds)

**Challenge**: Photolysis-driven equilibrium, diurnal variation (j changes), treating O atom as QSSA

**Tests**:
- [ ] Photostationary state achieved rapidly
- [ ] Correct Leighton ratio
- [ ] Response to NOx perturbation
- [ ] Diurnal cycle with time-varying j

---

## Stiff & Numerically Challenging

### 6. Extremely Stiff Equilibrium

**Description**: Test solver stability with extreme equilibrium constants.

**Reactions**:
```
A <-> B    K_eq = 1e12
B -> C     k = 1e-6 s-1  (slow loss)
```

**Constraint**:
- `[B] / [A] = K_eq = 1e12`
- Conservation: `[A] + [B] + [C] = total`

**Initial Conditions**:
- [A] = 1.0, [B] = 0, [C] = 0

**Expected Behavior**:
- Immediate equilibration: [B]/[A] = 1e12
- Slow drain to C over time

**Challenge**: 18 orders of magnitude range between fast equilibrium and slow kinetics, tests regularization approach

**Tests**:
- [ ] Equilibrium maintained despite extreme K
- [ ] No numerical instability
- [ ] Correct long-term evolution to C
- [ ] Works with different time steps

---

### 7. Temperature-Dependent Equilibrium

**Description**: Equilibrium constant varies with time via temperature dependence.

**Reactions**:
```
A + B <-> C    K(T) = K0 * exp(-Ea/R * (1/T - 1/T0))
```

**Parameters**:
- K0 = 1000 at T0 = 300 K
- Ea = 50 kJ/mol

**Temperature Profile**:
```
T(t) = 300 + 10*sin(2*pi*t/86400)  # Diurnal cycle
```

**Constraint**:
- `[C] / ([A]*[B]) = K(T(t))`

**Challenge**: Time-varying equilibrium constant, constraint changes each timestep

**Tests**:
- [ ] Equilibrium tracks temperature
- [ ] Correct K(T) relationship
- [ ] Smooth behavior through temperature cycle

---

### 8. Competing Equilibria

**Description**: Multiple species competing for binding, tests multi-constraint solver.

**Reactions**:
```
A + B <-> AB     K1 = 1000
A + C <-> AC     K2 = 100
B + C <-> BC     K3 = 10
```

**Constraints** (6 total):
- `[AB] / ([A]*[B]) = K1`
- `[AC] / ([A]*[C]) = K2`
- `[BC] / ([B]*[C]) = K3`
- `[A] + [AB] + [AC] = A_total`
- `[B] + [AB] + [BC] = B_total`
- `[C] + [AC] + [BC] = C_total`

**Initial Conditions**:
- A_total = B_total = C_total = 1.0

**Challenge**: 6 coupled constraints, 6 algebraic variables, potential for near-singular Jacobian

**Tests**:
- [ ] All equilibria satisfied simultaneously
- [ ] All conservation laws satisfied
- [ ] Unique solution found
- [ ] Stable with perturbations

---

## Practical Applications

### 9. ISORROPIA-style Aerosol Thermodynamics

**Description**: Simplified inorganic aerosol equilibrium, similar to the ISORROPIA model.

**Reactions**:
```
NH4NO3(s) <-> NH3(g) + HNO3(g)      Kp1
NH4Cl(s) <-> NH3(g) + HCl(g)        Kp2
(NH4)2SO4(s) <-> 2NH3(g) + H2SO4    Kp3 (solid stable)
```

**Constraints** (conditional on solid presence):
- If NH4NO3(s) present: `p_NH3 * p_HNO3 = Kp1`
- If NH4Cl(s) present: `p_NH3 * p_HCl = Kp2`
- Sulfate conservation: all sulfate as (NH4)2SO4

**Initial Conditions**:
- Total NH4 = 10 μg/m³
- Total NO3 = 5 μg/m³
- Total SO4 = 8 μg/m³

**Challenge**: Phase transitions, conditional constraints (solid present or dissolved), activity coefficients

**Tests**:
- [ ] Correct gas-aerosol partitioning
- [ ] Phase diagram reproduction
- [ ] Response to RH change

---

### 10. Chapman Mechanism + Fast O Equilibrium ⭐ NEXT

**Description**: Classic stratospheric ozone chemistry with atomic oxygen treated as quasi-steady-state algebraic variable.

**Reactions**:
```
O2 + hv -> 2O           j1 = 1e-12 s-1 (stratosphere)
O + O2 + M -> O3        k2 = 6.0e-34*(T/300)^-2.4 cm6 s-1
O3 + hv -> O2 + O       j3 = 1e-3 s-1
O + O3 -> 2O2           k4 = 8.0e-12*exp(-2060/T) cm3 s-1
```

**QSSA Constraint for O atom**:
```
d[O]/dt ≈ 0
2*j1*[O2] + j3*[O3] = k2*[O][O2][M] + k4*[O][O3]
```

Solving for [O]:
```
[O] = (2*j1*[O2] + j3*[O3]) / (k2*[O2][M] + k4*[O3])
```

**Odd Oxygen Conservation**:
```
[Ox] = [O] + [O3] ≈ [O3]  (since [O] << [O3])
```

**Initial Conditions** (20 km altitude):
- [O2] = 4e17 cm-3
- [O3] = 5e12 cm-3
- [M] = 2e18 cm-3
- T = 220 K

**Expected Behavior**:
- [O] ~ 1e7 cm-3 (ppt levels)
- O3 lifetime ~ weeks
- Diurnal cycle in [O] (j1, j3 vary)

**Challenge**:
- QSSA as true algebraic constraint
- Real atmospheric chemistry
- Stiff: [O]/[O3] ratio ~10^-5
- Tests the core DAE capability

**Tests**:
- [ ] QSSA constraint satisfied (d[O]/dt ≈ 0)
- [ ] Correct [O]/[O3] ratio
- [ ] O3 evolution matches analytical solution
- [ ] Diurnal variation with j(t)
- [ ] Comparison with pure ODE solution

**Implementation Notes**:
- Species: O2, O, O3 (O is algebraic)
- Constraint: `2*j1*[O2] + j3*[O3] - k2*[O][O2][M] - k4*[O][O3] = 0`
- This requires a new constraint type or CustomConstraint

---

## Edge Cases

### 11. Constraint Switching

**Description**: A constraint that becomes active/inactive based on system state.

**Reactions**:
```
A + B <-> C    K = 1000 (reversible)
C -> D         k = 0.1 s-1 (when [C] > threshold)
```

**Constraint** (conditional):
- When `[C] <= 0.5`: equilibrium constraint active
- When `[C] > 0.5`: equilibrium breaks, C drains to D

**Challenge**: Constraint becomes active/inactive, discontinuous behavior, event detection

**Tests**:
- [ ] Smooth transition at threshold
- [ ] Correct behavior in both regimes
- [ ] No oscillation at boundary

---

### 12. Near-Singular System

**Description**: Test detection and handling of redundant/singular constraints.

**Reactions**:
```
A <-> B    K1 = 1.0
B <-> C    K2 = 1.0
C <-> A    K3 = 1.0
```

Note: Thermodynamic consistency requires `K1 * K2 * K3 = 1`

**Constraints**:
- `[B]/[A] = K1`
- `[C]/[B] = K2`
- `[A]/[C] = K3` (redundant!)

**Challenge**: Third constraint is linearly dependent, Jacobian is singular, need to detect and handle

**Tests**:
- [ ] Detect redundant constraint
- [ ] Graceful handling (warning or automatic removal)
- [ ] Correct solution despite redundancy

---

## Test Implementation Priority

| Priority | Test | Reason |
|----------|------|--------|
| 1 | Chapman Mechanism (#10) | Core use case, real chemistry, tests QSSA |
| 2 | Extremely Stiff (#6) | Tests numerical robustness |
| 3 | Competing Equilibria (#8) | Tests multi-constraint coupling |
| 4 | O3-NOx PSS (#5) | Common atmospheric application |
| 5 | Carbonate System (#1) | Multi-phase, real chemistry |
| 6 | Temperature-Dependent (#7) | Time-varying constraints |
| 7 | Others | As needed |

---

## Notes

### CSV Output Format

All numerical integration tests should output CSV data to stdout for plotting and validation:

```cpp
// At start of test:
std::cout << "time";
for (const auto& name : state.variable_names_)
  std::cout << "," << name;
std::cout << ",O_analytical,O3_analytical" << std::endl;  // if available

// After each timestep:
std::cout << std::scientific << std::setprecision(6);
std::cout << time;
for (std::size_t i = 0; i < state.variables_[0].size(); ++i)
  std::cout << "," << state.variables_[0][i];
std::cout << "," << O_analytical << "," << O3_analytical << std::endl;
```

**Example output**:
```csv
time,O2,O,O3,O_analytical,O3_analytical
0.000000e+00,4.000000e+17,0.000000e+00,5.000000e+12,1.200000e+07,5.000000e+12
1.000000e-03,4.000000e+17,1.198234e+07,4.999998e+12,1.200000e+07,4.999999e+12
```

**Usage**:
```bash
# Run test and save CSV
./test_chapman > chapman_results.csv

# Plot with Python
python plot_results.py chapman_results.csv

# Quick plot with gnuplot
gnuplot -e "set datafile separator ','; plot 'chapman_results.csv' using 1:4 with lines"
```

### Analytical Solutions

Where possible, tests should compare against:
1. Analytical steady-state solutions
2. Reference numerical solutions (high-precision ODE solver)
3. Published results from literature

### Tolerances

Suggested tolerances for validation:
- Equilibrium ratios: within 1% of K_eq
- Conservation: within 1e-10 of initial total
- Steady state: d[X]/dt < 1e-12 * [X]

### Test File Location

Integration tests should go in:
```
test/integration/test_dae_*.cpp
```

Unit tests for new constraint types:
```
test/unit/constraint/test_*.cpp
```
