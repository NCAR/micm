# DAE Rosenbrock Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone benchmark that compares a full-ODE Robertson solve against a QSSA-reduced DAE solve across a stiffness sweep, emitting deterministic solver-stat counters + accuracy + wall-clock to CSV for paper/docs figures.

**Architecture:** Two reusable headers (Robertson reactions+ICs, and a custom QSSA constraint model) are shared by (a) a GoogleTest correctness guard and (b) a CSV-emitting benchmark executable. Both solver builds use identical Rosenbrock parameters and tolerances; the only difference is that the DAE build marks species B algebraic via the custom constraint that replaces B's ODE row with `dB/dt = 0`. Everything lives under `benchmark/` behind an OFF-by-default `MICM_ENABLE_BENCHMARKS` CMake option.

**Tech Stack:** C++20, MICM header-only API (`micm/CPU.hpp`, `CpuSolverBuilder`, `RosenbrockSolverParameters`), GoogleTest, CMake, Python+matplotlib for plotting.

---

## Background reference (the chemistry — do not re-derive)

Robertson system (rates are user-defined parameters `r1`,`r2`,`r3`):

```
r1:  A      -> B            rate = k1 * [A]            k1 = 0.04
r2:  2B     -> B + C        rate = k2 * [B]^2          k2 = 3e7  (stiffness knob)
r3:  B + C  -> A + C        rate = k3 * [B] * [C]      k3 = 1e4
```

Forcing:
```
dA/dt = -k1*A           + k3*B*C
dB/dt =  k1*A - k2*B^2  - k3*B*C
dC/dt =         k2*B^2
```

QSSA reduction replaces B's ODE row with the algebraic constraint:
```
G(A,B,C) = k1*A - k2*B^2 - k3*B*C = 0      (quadratic in B)
```
Constraint Jacobian (the solver SUBTRACTS dG/dy, matching ProcessSet convention):
```
dG/dA = k1
dG/dB = -2*k2*B - k3*C
dG/dC = -k3*B
```
Consistent initial B (with A=1, C=0): positive root of `k2*B^2 - k1*A = 0` => `B0 = sqrt(k1*A / k2)`.

Initial condition: A=1, B=0 (projected to B0), C=0. Conservation: A+B+C=1.

---

## File Structure

- Create `benchmark/robertson_system.hpp` — builds the gas phase, the three reactions, exposes initial conditions and the consistent-B helper. One source of truth for both solver builds.
- Create `benchmark/robertson_qssa_constraint.hpp` — the constraint-only external model implementing `G` and its Jacobian. No new state variables.
- Create `benchmark/robertson_dae.cpp` — `RunCase`, `ReferenceSolution`, the k2 sweep, and the CSV writer. The benchmark executable.
- Create `benchmark/test_robertson_qssa.cpp` — GoogleTest correctness guard: constraint residual derivative check + DAE-vs-ODE post-transient agreement + residual≈0.
- Create `benchmark/plot_robertson_dae.py` — reads the CSV, writes two figures.
- Create `benchmark/CMakeLists.txt` — the benchmark exe + the gtest, behind `MICM_ENABLE_BENCHMARKS`.
- Create `benchmark/README.md` — how to build/run/plot.
- Modify `CMakeLists.txt` — add the option and `add_subdirectory(benchmark)`.

---

## Task 1: CMake scaffolding behind MICM_ENABLE_BENCHMARKS

**Files:**
- Modify: `CMakeLists.txt` (option block near lines 27-34; subdirectory block near line 75)
- Create: `benchmark/CMakeLists.txt`
- Create: `benchmark/robertson_dae.cpp` (temporary stub for this task)

- [ ] **Step 1: Add the option.** In `CMakeLists.txt`, after the existing `option(MICM_BUILD_SHARED_LIBS ...)` line (line 34), add:

```cmake
option(MICM_ENABLE_BENCHMARKS "Build the DAE/ODE benchmark programs" OFF)
```

- [ ] **Step 2: Wire the subdirectory.** In `CMakeLists.txt`, immediately after the `if(PROJECT_IS_TOP_LEVEL AND MICM_ENABLE_TESTS)` block that ends with `add_subdirectory(test)` (around line 75-77), add a new block:

```cmake
if(PROJECT_IS_TOP_LEVEL AND MICM_ENABLE_BENCHMARKS)
  add_subdirectory(benchmark)
endif()
```

- [ ] **Step 3: Create the temporary stub** `benchmark/robertson_dae.cpp`:

```cpp
// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <iostream>

int main()
{
  std::cout << "robertson_dae benchmark stub\n";
  return 0;
}
```

- [ ] **Step 4: Create** `benchmark/CMakeLists.txt`:

```cmake
################################################################################
# DAE / ODE benchmarks (built only when MICM_ENABLE_BENCHMARKS=ON)

add_executable(robertson_dae robertson_dae.cpp)
target_link_libraries(robertson_dae PUBLIC musica::micm)
target_compile_features(robertson_dae PUBLIC cxx_std_20)
```

- [ ] **Step 5: Configure and build the stub.**

Run:
```bash
cd /Users/fillmore/EarthSystem/MICM && cmake -S . -B build -DMICM_ENABLE_BENCHMARKS=ON -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER >/dev/null && cmake --build build --target robertson_dae 2>&1 | tail -5
```
Expected: builds cleanly, produces `build/benchmark/robertson_dae`.

> NOTE: `-DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER` is required in this repo so FetchContent builds GTest from source — without it tests abort (see memory `build-gtest-conda-rpath`).

- [ ] **Step 6: Run the stub.**

Run: `cd /Users/fillmore/EarthSystem/MICM && ./build/benchmark/robertson_dae`
Expected: prints `robertson_dae benchmark stub`.

- [ ] **Step 7: Commit.**

```bash
cd /Users/fillmore/EarthSystem/MICM
git add CMakeLists.txt benchmark/CMakeLists.txt benchmark/robertson_dae.cpp
git commit -m "Add benchmark subdirectory scaffolding behind MICM_ENABLE_BENCHMARKS"
```

---

## Task 2: Robertson system header

**Files:**
- Create: `benchmark/robertson_system.hpp`

- [ ] **Step 1: Write the header.** This is one source of truth for the chemistry, used by both solver builds. Create `benchmark/robertson_system.hpp`:

```cpp
// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Robertson stiff problem definition, shared by the full-ODE and DAE benchmark
// solver builds so both run provably identical chemistry.
//
//   r1: A      -> B         rate = k1 * [A]
//   r2: 2B     -> B + C     rate = k2 * [B]^2
//   r3: B + C  -> A + C     rate = k3 * [B] * [C]
#pragma once

#include <micm/CPU.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace robertson
{
  /// Default Robertson rate constants (Hairer & Wanner II, p.3).
  inline constexpr double K1_DEFAULT = 0.04;
  inline constexpr double K2_DEFAULT = 3.0e7;  // the stiffness knob
  inline constexpr double K3_DEFAULT = 1.0e4;

  /// The species and reactions of the Robertson system. The gas phase must be
  /// kept alive for the lifetime of the solver, so it is returned alongside the
  /// processes.
  struct System
  {
    micm::Phase gas_phase;
    std::vector<micm::Process> processes;
  };

  /// Build the Robertson reaction system. Rate constants are supplied at solve
  /// time via SetCustomRateParameter("r1"/"r2"/"r3", ...), not baked in here.
  inline System MakeSystem()
  {
    auto a = micm::Species("A");
    auto b = micm::Species("B");
    auto c = micm::Species("C");

    micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b, c } };

    micm::Process r1 = micm::ChemicalReactionBuilder()
                           .SetReactants({ a })
                           .SetProducts({ micm::StoichSpecies(b, 1) })
                           .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r1" })
                           .SetPhase(gas_phase)
                           .Build();

    micm::Process r2 = micm::ChemicalReactionBuilder()
                           .SetReactants({ b, b })
                           .SetProducts({ micm::StoichSpecies(b, 1), micm::StoichSpecies(c, 1) })
                           .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r2" })
                           .SetPhase(gas_phase)
                           .Build();

    micm::Process r3 = micm::ChemicalReactionBuilder()
                           .SetReactants({ b, c })
                           .SetProducts({ micm::StoichSpecies(a, 1), micm::StoichSpecies(c, 1) })
                           .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r3" })
                           .SetPhase(gas_phase)
                           .Build();

    return System{ gas_phase, { r1, r2, r3 } };
  }

  /// Consistent quasi-steady-state value of B given A and C (positive root of
  /// k2*B^2 + k3*C*B - k1*A = 0). Used to project the DAE initial condition onto
  /// the constraint manifold (B(0)=0 is otherwise inconsistent with G=0).
  inline double ConsistentB(double k1, double k2, double k3, double a, double c)
  {
    double disc = k3 * c * k3 * c + 4.0 * k2 * k1 * a;
    return (-k3 * c + std::sqrt(disc)) / (2.0 * k2);
  }
}  // namespace robertson
```

- [ ] **Step 2: Verify it compiles** by including it from the stub. Temporarily replace `benchmark/robertson_dae.cpp` body with:

```cpp
// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "robertson_system.hpp"

#include <iostream>

int main()
{
  auto sys = robertson::MakeSystem();
  std::cout << "reactions: " << sys.processes.size()
            << " B0=" << robertson::ConsistentB(0.04, 3e7, 1e4, 1.0, 0.0) << "\n";
  return 0;
}
```

- [ ] **Step 3: Build and run.**

Run:
```bash
cd /Users/fillmore/EarthSystem/MICM && cmake --build build --target robertson_dae 2>&1 | tail -5 && ./build/benchmark/robertson_dae
```
Expected: prints `reactions: 3 B0=3.65148e-05` (B0 = sqrt(0.04/3e7) ≈ 3.65e-5).

- [ ] **Step 4: Commit.**

```bash
cd /Users/fillmore/EarthSystem/MICM
git add benchmark/robertson_system.hpp benchmark/robertson_dae.cpp
git commit -m "Add Robertson system header for benchmark"
```

---

## Task 3: QSSA constraint model + correctness gtest (TDD)

**Files:**
- Create: `benchmark/robertson_qssa_constraint.hpp`
- Create: `benchmark/test_robertson_qssa.cpp`
- Modify: `benchmark/CMakeLists.txt`

- [ ] **Step 1: Write the failing test first.** Create `benchmark/test_robertson_qssa.cpp`. This builds the DAE solver and asserts the QSSA solution tracks a tight-tolerance full-ODE reference after the transient and keeps the constraint residual near zero. It fails to compile until the constraint header exists.

```cpp
// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "robertson_qssa_constraint.hpp"
#include "robertson_system.hpp"

#include <micm/CPU.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace
{
  // Integrate Robertson from t=0 to t_end, returning {A,B,C} at t_end.
  // If qssa==true, B is algebraic via the QSSA constraint; otherwise full ODE.
  std::vector<double> Integrate(bool qssa, double k1, double k2, double k3, double rtol, double t_end)
  {
    auto sys = robertson::MakeSystem();
    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();

    // The constraint is copied by value into AddExternalModel; keep it in scope
    // through the build expression. Build in one unbroken chain (never store the
    // intermediate builder in `auto` — the chained setters return a base
    // SolverBuilder& and would slice). Both branches produce the SAME solver
    // type because AddExternalModel returns SolverBuilder& (type-erased), so the
    // ternary is well-typed.
    robertson::QssaConstraint constraint(k1, k2, k3);

    auto solver =
        qssa ? micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                   .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = sys.gas_phase }))
                   .SetReactions(sys.processes)
                   .AddExternalModel(constraint)
                   .SetReorderState(false)
                   .Build()
             : micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                   .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = sys.gas_phase }))
                   .SetReactions(sys.processes)
                   .SetReorderState(false)
                   .Build();

    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.variables_[0].size(), rtol * 1e-2));
    state.SetCustomRateParameter("r1", k1);
    state.SetCustomRateParameter("r2", k2);
    state.SetCustomRateParameter("r3", k3);

    auto map = state.variable_map_;
    state.variables_[0][map.at("A")] = 1.0;
    state.variables_[0][map.at("B")] = qssa ? robertson::ConsistentB(k1, k2, k3, 1.0, 0.0) : 0.0;
    state.variables_[0][map.at("C")] = 0.0;
    state.conditions_[0].temperature_ = 272.5;
    state.conditions_[0].pressure_ = 101253.3;
    state.conditions_[0].air_density_ = 1e6;
    solver.UpdateStateParameters(state);

    double done = 0.0;
    while (done < t_end)
    {
      auto result = solver.Solve(t_end - done, state);
      EXPECT_EQ(result.state_, micm::SolverState::Converged);
      done += result.stats_.final_time_;
    }
    return { state.variables_[0][map.at("A")],
             state.variables_[0][map.at("B")],
             state.variables_[0][map.at("C")] };
  }
}  // namespace

// The QSSA-reduced DAE should match the full ODE for A and C after the transient.
TEST(RobertsonQssa, MatchesFullOdePostTransient)
{
  const double k1 = 0.04, k2 = 3e7, k3 = 1e4, t_end = 1.0e3;
  auto ode = Integrate(false, k1, k2, k3, 1e-10, t_end);  // tight-tol reference
  auto dae = Integrate(true, k1, k2, k3, 1e-6, t_end);

  auto rel = [](double got, double ref) { return std::abs(got - ref) / (std::abs(ref) + 1e-30); };
  EXPECT_LT(rel(dae[0], ode[0]), 1e-3) << "A: dae=" << dae[0] << " ode=" << ode[0];
  EXPECT_LT(rel(dae[2], ode[2]), 1e-3) << "C: dae=" << dae[2] << " ode=" << ode[2];
}

// The QSSA constraint residual G = k1*A - k2*B^2 - k3*B*C must be ~0 at the end.
TEST(RobertsonQssa, ConstraintResidualNearZero)
{
  const double k1 = 0.04, k2 = 3e7, k3 = 1e4;
  auto dae = Integrate(true, k1, k2, k3, 1e-6, 1.0e3);
  double g = k1 * dae[0] - k2 * dae[1] * dae[1] - k3 * dae[1] * dae[2];
  EXPECT_NEAR(g, 0.0, 1e-6) << "residual G=" << g;
}
```

- [ ] **Step 2: Add the gtest target.** Append to `benchmark/CMakeLists.txt`:

```cmake
add_executable(test_robertson_qssa test_robertson_qssa.cpp)
target_link_libraries(test_robertson_qssa PUBLIC musica::micm GTest::gtest_main)
target_compile_features(test_robertson_qssa PUBLIC cxx_std_20)
add_test(NAME test_robertson_qssa COMMAND test_robertson_qssa)
```

- [ ] **Step 3: Run to confirm it fails (missing header).**

Run:
```bash
cd /Users/fillmore/EarthSystem/MICM && cmake -S . -B build -DMICM_ENABLE_BENCHMARKS=ON -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER >/dev/null && cmake --build build --target test_robertson_qssa 2>&1 | tail -5
```
Expected: FAIL — `robertson_qssa_constraint.hpp: No such file or directory`.

- [ ] **Step 4: Write the constraint header.** Create `benchmark/robertson_qssa_constraint.hpp`. It mirrors the constraint-only external model interface in `test/integration/test_external_model_constraints.cpp` (EquilibriumConstraintModel). B is the algebraic species; its ODE row is replaced by `G=0`. The Jacobian lambda receives `(state_variables, state_parameters, jac)` and reads B,C from the first argument.

```cpp
// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Constraint-only external model implementing the Robertson QSSA on species B:
//   G(A,B,C) = k1*A - k2*B^2 - k3*B*C = 0   (replaces dB/dt = 0)
//
// B, A, C are pre-existing species, so this model declares no new state
// variables and is wired via SolverBuilder::AddExternalModel only (no entry in
// SystemParameters::external_models_).
#pragma once

#include <micm/system/conditions.hpp>

#include <functional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace robertson
{
  class QssaConstraint
  {
   public:
    QssaConstraint(double k1, double k2, double k3)
        : k1_(k1),
          k2_(k2),
          k3_(k3)
    {
    }

    // B's ODE row becomes algebraic.
    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { "B" };
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      return { "A", "B", "C" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& state_indices) const
    {
      auto i_a = state_indices.at("A");
      auto i_b = state_indices.at("B");
      auto i_c = state_indices.at("C");
      // Constraint row is B; depends on A, B, C.
      return { { i_b, i_a }, { i_b, i_b }, { i_b, i_c } };
    }

    std::set<std::string> ConstraintStateParameterNames() const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>
    ConstraintUpdateStateParametersFunction(const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
    }

    // Residual: G = k1*A - k2*B^2 - k3*B*C  (written into B's forcing row).
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& var) const
    {
      auto i_a = var.at("A");
      auto i_b = var.at("B");
      auto i_c = var.at("C");
      double k1 = k1_, k2 = k2_, k3 = k3_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, DenseMatrixPolicy& forcing)
      {
        for (std::size_t i = 0; i < state.NumRows(); ++i)
        {
          double a = state[i][i_a], b = state[i][i_b], c = state[i][i_c];
          forcing[i][i_b] = k1 * a - k2 * b * b - k3 * b * c;
        }
      };
    }

    // Jacobian (solver subtracts dG/dy):
    //   dG/dA = k1 ; dG/dB = -2*k2*B - k3*C ; dG/dC = -k3*B
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& var,
        const SparseMatrixPolicy&) const
    {
      auto i_a = var.at("A");
      auto i_b = var.at("B");
      auto i_c = var.at("C");
      double k1 = k1_, k2 = k2_, k3 = k3_;
      return [=](const DenseMatrixPolicy& state, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
      {
        for (std::size_t i = 0; i < jac.NumberOfBlocks(); ++i)
        {
          double b = state[i][i_b], c = state[i][i_c];
          jac[i][i_b][i_a] -= k1;
          jac[i][i_b][i_b] -= (-2.0 * k2 * b - k3 * c);
          jac[i][i_b][i_c] -= (-k3 * b);
        }
      };
    }

   private:
    double k1_;
    double k2_;
    double k3_;
  };
}  // namespace robertson
```

- [ ] **Step 5: Build and run the test.**

Run:
```bash
cd /Users/fillmore/EarthSystem/MICM && cmake --build build --target test_robertson_qssa 2>&1 | tail -5 && ./build/benchmark/test_robertson_qssa
```
Expected: PASS — both `RobertsonQssa.MatchesFullOdePostTransient` and `RobertsonQssa.ConstraintResidualNearZero`.

- [ ] **Step 6: Restore the benchmark stub** so the `robertson_dae` target still builds (it will be filled in Task 4). Replace `benchmark/robertson_dae.cpp` body with the Task 1 stub printing `robertson_dae benchmark stub`, then build:

Run: `cd /Users/fillmore/EarthSystem/MICM && cmake --build build --target robertson_dae 2>&1 | tail -3`
Expected: builds cleanly.

- [ ] **Step 7: Commit.**

```bash
cd /Users/fillmore/EarthSystem/MICM
git add benchmark/robertson_qssa_constraint.hpp benchmark/test_robertson_qssa.cpp benchmark/CMakeLists.txt benchmark/robertson_dae.cpp
git commit -m "Add Robertson QSSA constraint with correctness gtest"
```

---

## Task 4: Benchmark core — RunCase, ReferenceSolution, stats aggregation

**Files:**
- Modify: `benchmark/robertson_dae.cpp`

- [ ] **Step 1: Write the full benchmark program.** Replace the entire contents of `benchmark/robertson_dae.cpp`:

```cpp
// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Benchmark: full-ODE vs QSSA-DAE Robertson across a stiffness (k2) sweep.
// Emits deterministic solver-stat counters (primary evidence), median
// wall-clock (secondary, hardware-dependent), and post-transient accuracy of
// the DAE vs a tight-tolerance full-ODE reference. Writes one CSV.
#include "robertson_qssa_constraint.hpp"
#include "robertson_system.hpp"

#include <micm/CPU.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  enum class Method
  {
    FullOde,
    DaeQssa
  };

  // Aggregated solver counters plus the trajectory sampled at the output times.
  struct CaseResult
  {
    std::uint64_t number_of_steps = 0;
    std::uint64_t accepted = 0;
    std::uint64_t rejected = 0;
    std::uint64_t function_calls = 0;
    std::uint64_t jacobian_updates = 0;
    std::uint64_t decompositions = 0;
    std::uint64_t solves = 0;
    bool converged = true;
    double wallclock_median_us = 0.0;
    std::vector<double> a_at_output;  // one entry per output time
    std::vector<double> c_at_output;
  };

  // Logarithmically spaced output times in [t_first, t_last].
  std::vector<double> OutputTimes(double t_first, double t_last, int n)
  {
    std::vector<double> times;
    double log_first = std::log10(t_first);
    double log_last = std::log10(t_last);
    for (int i = 0; i < n; ++i)
      times.push_back(std::pow(10.0, log_first + (log_last - log_first) * i / (n - 1)));
    return times;
  }

  // Run one integration to completion, accumulating solver stats across every
  // Solve() call and recording A and C at each output time. Stats are gathered
  // only on the timed pass (reps==1 for that); extra reps are wall-clock only.
  CaseResult RunCase(Method method, double k1, double k2, double k3, double rtol, const std::vector<double>& output_times, int wallclock_reps)
  {
    auto sys = robertson::MakeSystem();
    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    // Build in one unbroken chain (never store the intermediate builder in
    // `auto` — chained setters return a base SolverBuilder& and would slice).
    // Both branches yield the SAME solver type (AddExternalModel returns
    // SolverBuilder&, type-erased), so the ternary is well-typed.
    robertson::QssaConstraint constraint(k1, k2, k3);

    auto solver =
        (method == Method::DaeQssa)
            ? micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                  .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = sys.gas_phase }))
                  .SetReactions(sys.processes)
                  .AddExternalModel(constraint)
                  .SetReorderState(false)
                  .Build()
            : micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                  .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = sys.gas_phase }))
                  .SetReactions(sys.processes)
                  .SetReorderState(false)
                  .Build();

    auto init_state = [&](auto& state)
    {
      state.SetRelativeTolerance(rtol);
      state.SetAbsoluteTolerances(std::vector<double>(state.variables_[0].size(), rtol * 1e-2));
      state.SetCustomRateParameter("r1", k1);
      state.SetCustomRateParameter("r2", k2);
      state.SetCustomRateParameter("r3", k3);
      auto map = state.variable_map_;
      state.variables_[0][map.at("A")] = 1.0;
      state.variables_[0][map.at("B")] = (method == Method::DaeQssa) ? robertson::ConsistentB(k1, k2, k3, 1.0, 0.0) : 0.0;
      state.variables_[0][map.at("C")] = 0.0;
      state.conditions_[0].temperature_ = 272.5;
      state.conditions_[0].pressure_ = 101253.3;
      state.conditions_[0].air_density_ = 1e6;
      solver.UpdateStateParameters(state);
    };

    CaseResult out;

    // ---- Stats + trajectory pass (untimed) ----
    {
      auto state = solver.GetState(1);
      init_state(state);
      auto map = state.variable_map_;
      double current = 0.0;
      for (double t_out : output_times)
      {
        double done = 0.0;
        double dt = t_out - current;
        while (done < dt)
        {
          auto result = solver.Solve(dt - done, state);
          if (result.state_ != micm::SolverState::Converged)
          {
            out.converged = false;
            break;
          }
          out.number_of_steps += result.stats_.number_of_steps_;
          out.accepted += result.stats_.accepted_;
          out.rejected += result.stats_.rejected_;
          out.function_calls += result.stats_.function_calls_;
          out.jacobian_updates += result.stats_.jacobian_updates_;
          out.decompositions += result.stats_.decompositions_;
          out.solves += result.stats_.solves_;
          done += result.stats_.final_time_;
        }
        out.a_at_output.push_back(state.variables_[0][map.at("A")]);
        out.c_at_output.push_back(state.variables_[0][map.at("C")]);
        current = t_out;
        if (!out.converged)
          break;
      }
    }

    // ---- Wall-clock passes (timed, stats ignored) ----
    std::vector<double> samples;
    for (int rep = 0; rep < wallclock_reps; ++rep)
    {
      auto state = solver.GetState(1);
      init_state(state);
      auto t0 = std::chrono::steady_clock::now();
      double current = 0.0;
      for (double t_out : output_times)
      {
        double done = 0.0;
        double dt = t_out - current;
        while (done < dt)
        {
          auto result = solver.Solve(dt - done, state);
          if (result.state_ != micm::SolverState::Converged)
            break;
          done += result.stats_.final_time_;
        }
        current = t_out;
      }
      auto t1 = std::chrono::steady_clock::now();
      samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    if (!samples.empty())
    {
      std::sort(samples.begin(), samples.end());
      out.wallclock_median_us = samples[samples.size() / 2];
    }
    return out;
  }

  // Max relative error of A and C of `cand` vs `ref`, over output times >= t_skip.
  double MaxRelErrorPostTransient(const CaseResult& cand, const CaseResult& ref, const std::vector<double>& output_times, double t_skip)
  {
    double worst = 0.0;
    for (std::size_t i = 0; i < output_times.size(); ++i)
    {
      if (output_times[i] < t_skip)
        continue;
      if (i >= cand.a_at_output.size() || i >= ref.a_at_output.size())
        break;
      auto rel = [](double got, double r) { return std::abs(got - r) / (std::abs(r) + 1e-30); };
      worst = std::max(worst, rel(cand.a_at_output[i], ref.a_at_output[i]));
      worst = std::max(worst, rel(cand.c_at_output[i], ref.c_at_output[i]));
    }
    return worst;
  }
}  // namespace

int main()
{
  const double k1 = robertson::K1_DEFAULT;
  const double k3 = robertson::K3_DEFAULT;
  const double rtol = 1e-6;
  const double t_skip = 1e-2;       // transient cutoff: B settles well before this
  const int wallclock_reps = 7;     // odd -> well-defined median
  const std::vector<double> k2_sweep = { 3e5, 3e6, 3e7, 3e8, 3e9 };
  const auto output_times = OutputTimes(1e-3, 1e6, 19);  // spans transient + long tail

  std::ofstream csv("robertson_dae_benchmark.csv");
  csv << "method,k2,rtol,number_of_steps,accepted,rejected,function_calls,"
         "jacobian_updates,decompositions,solves,converged,wallclock_median_us,max_rel_err\n";

  auto write_row = [&](const std::string& method, double k2, const CaseResult& r, double max_rel_err)
  {
    csv << method << ',' << k2 << ',' << rtol << ',' << r.number_of_steps << ',' << r.accepted << ',' << r.rejected << ','
        << r.function_calls << ',' << r.jacobian_updates << ',' << r.decompositions << ',' << r.solves << ','
        << (r.converged ? 1 : 0) << ',' << r.wallclock_median_us << ',' << max_rel_err << '\n';
  };

  for (double k2 : k2_sweep)
  {
    auto reference = RunCase(Method::FullOde, k1, k2, k3, 1e-10, output_times, 0);  // tight-tol accuracy reference
    auto ode = RunCase(Method::FullOde, k1, k2, k3, rtol, output_times, wallclock_reps);
    auto dae = RunCase(Method::DaeQssa, k1, k2, k3, rtol, output_times, wallclock_reps);

    double ode_err = MaxRelErrorPostTransient(ode, reference, output_times, t_skip);
    double dae_err = MaxRelErrorPostTransient(dae, reference, output_times, t_skip);

    write_row("full_ode", k2, ode, ode_err);
    write_row("dae_qssa", k2, dae, dae_err);

    std::cout << "k2=" << k2 << "  ODE steps=" << ode.number_of_steps << " (conv=" << ode.converged << ")"
              << "  DAE steps=" << dae.number_of_steps << " (conv=" << dae.converged << ")"
              << "  DAE max_rel_err=" << dae_err << "\n";
  }

  std::cout << "wrote robertson_dae_benchmark.csv\n";
  return 0;
}
```

- [ ] **Step 2: Build.**

Run: `cd /Users/fillmore/EarthSystem/MICM && cmake --build build --target robertson_dae 2>&1 | tail -5`
Expected: builds cleanly.

- [ ] **Step 3: Run and inspect.**

Run: `cd /Users/fillmore/EarthSystem/MICM && cd build/benchmark && ./robertson_dae && head -20 robertson_dae_benchmark.csv`
Expected: prints one line per k2 where **DAE steps are far below ODE steps and the gap widens as k2 grows**, DAE `max_rel_err` stays small (≲1e-3), and a CSV with 10 data rows (5 k2 × 2 methods) is written. Note any non-converged ODE rows at high k2 — that is itself a data point, not a failure.

- [ ] **Step 4: Commit.**

```bash
cd /Users/fillmore/EarthSystem/MICM
git add benchmark/robertson_dae.cpp
git commit -m "Implement Robertson DAE vs ODE stiffness-sweep benchmark"
```

---

## Task 5: Plotting script

**Files:**
- Create: `benchmark/plot_robertson_dae.py`

- [ ] **Step 1: Write the script.** Create `benchmark/plot_robertson_dae.py`:

```python
#!/usr/bin/env python3
# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
"""Plot Robertson DAE-vs-ODE benchmark results.

Usage:
    python3 plot_robertson_dae.py robertson_dae_benchmark.csv
Produces:
    robertson_cost_vs_stiffness.png   (steps & Jacobian updates vs k2)
    robertson_accuracy_vs_stiffness.png
"""
import csv
import sys

import matplotlib.pyplot as plt


def load(path):
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def series(rows, method, x_key, y_key, conv_only=False):
    xs, ys = [], []
    for r in rows:
        if r["method"] != method:
            continue
        if conv_only and r["converged"] != "1":
            continue
        xs.append(float(r[x_key]))
        ys.append(float(r[y_key]))
    pairs = sorted(zip(xs, ys))
    return [p[0] for p in pairs], [p[1] for p in pairs]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "robertson_dae_benchmark.csv"
    rows = load(path)

    # Figure 1: cost vs stiffness
    fig, ax = plt.subplots(figsize=(7, 5))
    for method, label in [("full_ode", "full ODE"), ("dae_qssa", "DAE (QSSA)")]:
        x, y = series(rows, method, "k2", "number_of_steps")
        ax.loglog(x, y, "o-", label=f"{label}: steps")
        x, y = series(rows, method, "k2", "jacobian_updates")
        ax.loglog(x, y, "s--", label=f"{label}: Jacobian updates")
    ax.set_xlabel("k2 (stiffness)")
    ax.set_ylabel("count (summed over integration)")
    ax.set_title("Solver cost vs stiffness: Robertson")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig("robertson_cost_vs_stiffness.png", dpi=150)

    # Figure 2: DAE accuracy vs stiffness
    fig, ax = plt.subplots(figsize=(7, 5))
    x, y = series(rows, "dae_qssa", "k2", "max_rel_err")
    ax.loglog(x, y, "o-", color="C2", label="DAE max post-transient rel. error")
    ax.set_xlabel("k2 (stiffness)")
    ax.set_ylabel("max relative error vs tight-tol ODE")
    ax.set_title("QSSA reduction accuracy vs stiffness")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig("robertson_accuracy_vs_stiffness.png", dpi=150)

    print("wrote robertson_cost_vs_stiffness.png and robertson_accuracy_vs_stiffness.png")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on the CSV.**

Run:
```bash
cd /Users/fillmore/EarthSystem/MICM/build/benchmark && python3 /Users/fillmore/EarthSystem/MICM/benchmark/plot_robertson_dae.py robertson_dae_benchmark.csv
```
Expected: prints the "wrote ..." line and produces two PNGs. (If matplotlib is missing, `pip install matplotlib`; this is a dev-only dependency, not part of the build.)

- [ ] **Step 3: Commit.**

```bash
cd /Users/fillmore/EarthSystem/MICM
git add benchmark/plot_robertson_dae.py
git commit -m "Add plotting script for Robertson DAE benchmark"
```

---

## Task 6: README and final verification

**Files:**
- Create: `benchmark/README.md`

- [ ] **Step 1: Write the README.** Create `benchmark/README.md`:

```markdown
# DAE Rosenbrock Benchmarks

Compares a full-ODE solve of the stiff Robertson problem against a QSSA-reduced
DAE solve (species B made algebraic via a custom constraint that enforces
`dB/dt = 0`), swept across the stiffness knob `k2`.

The headline evidence is the **deterministic solver-stat counters** (steps,
accepted/rejected, function calls, Jacobian updates) — these reproduce exactly
on any machine. Wall-clock (median over repetitions) is reported as secondary,
hardware-dependent context. Accuracy is the DAE's max post-transient relative
error against a tight-tolerance full-ODE reference.

> The QSSA reduction is an approximation: it discards B's initial fast transient.
> The benchmark therefore reports both the cost savings **and** the post-transient
> accuracy, so the tradeoff is explicit.

## Build

```bash
cmake -S . -B build -DMICM_ENABLE_BENCHMARKS=ON -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER
cmake --build build --target robertson_dae
cmake --build build --target test_robertson_qssa   # correctness guard
```

## Run

```bash
cd build/benchmark
./test_robertson_qssa          # verify the QSSA constraint is correct
./robertson_dae                # writes robertson_dae_benchmark.csv
python3 ../../benchmark/plot_robertson_dae.py robertson_dae_benchmark.csv
```

## Files

- `robertson_system.hpp` — Robertson reactions/ICs (shared by both solver builds)
- `robertson_qssa_constraint.hpp` — custom QSSA constraint on B
- `robertson_dae.cpp` — sweep + CSV writer
- `test_robertson_qssa.cpp` — correctness guard (Jacobian/agreement/residual)
- `plot_robertson_dae.py` — figures from the CSV

## Future work

Phase 2 (separate spec): a reversible-equilibrium case where the algebraic
`EquilibriumConstraint` is exact in the stiff limit. The sweep, CSV schema, and
plotting here are intended to be reused for it.
```

- [ ] **Step 2: Full clean verification.** Reconfigure from scratch and run everything to confirm the whole feature builds and passes.

Run:
```bash
cd /Users/fillmore/EarthSystem/MICM && cmake --build build --target test_robertson_qssa robertson_dae 2>&1 | tail -5 && ./build/benchmark/test_robertson_qssa && (cd build/benchmark && ./robertson_dae)
```
Expected: gtest reports all PASS; benchmark prints the per-k2 summary showing DAE steps ≪ ODE steps with the gap widening as k2 increases.

- [ ] **Step 3: Commit.**

```bash
cd /Users/fillmore/EarthSystem/MICM
git add benchmark/README.md
git commit -m "Add benchmark README and finalize Robertson DAE benchmark"
```

---

## Self-Review notes (for the implementer)

- **Fair comparison invariant:** both builds use `FourStageDifferentialAlgebraicRosenbrockParameters`, identical rtol/atol, identical output-time grid, identical ICs. The ONLY difference is `AddExternalModel(constraint)` in the DAE build. Do not change Rosenbrock parameters between the two — that would invalidate the step-count comparison.
- **Rate/constraint consistency:** `k1,k2,k3` are passed to both `SetCustomRateParameter` (kinetics) and the `QssaConstraint` constructor (residual). They must come from the same variables in `RunCase` — they do.
- **No silent truncation:** every k2 produces two CSV rows even on non-convergence (`converged` column flags it); the stats pass breaks out but still records the row.
- **Constraint Jacobian sign:** entries are SUBTRACTED (`jac -= dG/dy`), matching `ProcessSet` and the stub models. Getting a sign wrong will surface as Task 3's `MatchesFullOdePostTransient` failing or non-convergence.
```
