// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Cloud chemistry test over the range of atmospheric conditions.
//
// StubCloudChemistry is an external model in the total formulation. The total A_T is differential,
// and two constraints divide it between A_G and A_AQ with a Henry's law equilibrium. The aqueous
// reaction A_AQ -> P_AQ also writes into the algebraic A_AQ row. This is the structure of gas-aqueous
// cloud chemistry (https://github.com/NCAR/musica/issues/956) and of
// https://github.com/NCAR/micm/issues/1094.
//
// Because the equilibrium holds at all times, the total has the exact solution
//   A_T(t) = A0 * exp(-k * K_p / (1 + K_p) * t)
// with [A_G] = A_T / (1 + K_p), [A_AQ] = A_T * K_p / (1 + K_p), and [P_AQ] = A0 - A_T.
//
// The test sweeps temperature, pressure, and liquid water content over the ranges in the musica
// issue. The partition coefficient K_p then covers about eight orders of magnitude.
//
// Shared by the CPU and Kokkos drivers. The CopyToDevice()/CopyToHost() calls are no-ops on the
// CPU matrix types, so both backends run the identical call sequence.
#pragma once

#include "stub_cloud_chemistry.hpp"

#include <micm/CPU.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <sstream>
#include <type_traits>
#include <vector>

/// @brief Solves the cloud chemistry system at each combination of conditions and compares to the exact solution
/// @param make_builder Returns a solver builder for the given solver parameters
template<class BuilderFactory>
void TestCloudChemistryConditionSweep(BuilderFactory make_builder)
{
  constexpr bool IS_DOUBLE = std::is_same_v<micm::Real, double>;

  // Chemistry similar to H2O2 in cloud water
  const micm::Real k = 10.0;                    // aqueous rate constant (s-1)
  const micm::Real h_ref = 1.0e5;               // Henry's law constant at 298.15 K (mol L-1 atm-1)
  const micm::Real dh_r = 7300.0;               // temperature dependence of the Henry's law constant (K)
  const micm::Real mixing_ratio = 1.0e-9;       // initial total A (mol mol-1)
  const micm::Real gas_constant = 8.314462618;  // J mol-1 K-1

  // Ranges from https://github.com/NCAR/musica/issues/956
  const std::vector<micm::Real> temperatures{ 233.0, 260.0, 280.0, 285.0, 286.0, 300.0 };  // K
  const std::vector<micm::Real> pressures{ 20000.0, 70000.0, 85000.0, 100000.0 };          // Pa
  const std::vector<micm::Real> liquid_water_contents{ 1.0e-7, 3.0e-5, 3.0e-4, 1.0e-2 };   // kg m-3

  const micm::Real time_step = 5.0;  // s
  const micm::Index number_of_steps = 12;

  // Solver tolerances, and the tolerances for the comparison to the exact solution
  const micm::Real solver_rtol = IS_DOUBLE ? 1.0e-6 : 1.0e-4;
  // Float cannot resolve an absolute tolerance near its roundoff level, so relax it in float mode
  const micm::Real solver_atol_factor = IS_DOUBLE ? 1.0e-6 : 1.0e-4;
  const micm::Real compare_rtol = IS_DOUBLE ? 1.0e-4 : 1.0e-2;

  auto A_G = micm::Species("A_G");
  micm::Phase gas_phase{ "gas", { A_G } };

  for (const auto lwc : liquid_water_contents)
  {
    StubCloudChemistry cloud(k, h_ref, dh_r, lwc);

    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    auto solver = make_builder(options)
                      .SetSystem(micm::System(gas_phase))
                      .SetReactions({})
                      .SetReorderState(false)
                      .AddExternalModel(cloud)
                      .Build();

    for (const auto temperature : temperatures)
    {
      for (const auto pressure : pressures)
      {
        std::ostringstream trace;
        trace << "T = " << temperature << " K, P = " << pressure << " Pa, LWC = " << lwc << " kg m-3";
        SCOPED_TRACE(trace.str());

        const micm::Real k_p = cloud.PartitionCoefficient(temperature);
        const micm::Real a0 = mixing_ratio * pressure / (gas_constant * temperature);  // mol m-3
        const micm::Real compare_atol = compare_rtol * a0;

        auto state = solver.GetState(1);
        state.SetRelativeTolerance(solver_rtol);
        state.SetAbsoluteTolerances(std::vector<micm::Real>(state.state_size_, solver_atol_factor * a0));

        const auto i_total = state.variable_map_.at("CLOUD.A_T");
        const auto i_gas = state.variable_map_.at("A_G");
        const auto i_aq = state.variable_map_.at("CLOUD.A_AQ");
        const auto i_p = state.variable_map_.at("CLOUD.P_AQ");

        // Start on the equilibrium manifold, so that the initialization does not move mass
        state.variables_[0][i_total] = a0;
        state.variables_[0][i_gas] = a0 / (1.0 + k_p);
        state.variables_[0][i_aq] = a0 * k_p / (1.0 + k_p);
        state.variables_[0][i_p] = 0.0;
        state.conditions_[0].temperature_ = temperature;
        state.conditions_[0].pressure_ = pressure;

        state.variables_.CopyToDevice();
        state.conditions_.CopyToDevice();
        state.custom_rate_parameters_.CopyToDevice();
        solver.UpdateStateParameters(state);

        const micm::Real effective_rate = k * k_p / (1.0 + k_p);
        micm::Real time = 0.0;
        for (micm::Index step = 0; step < number_of_steps; ++step)
        {
          auto result = solver.Solve(time_step, state);
          state.variables_.CopyToHost();
          ASSERT_EQ(result.state_, micm::SolverState::Converged) << "Step " << step;
          time += time_step;

          const micm::Real total = a0 * std::exp(-effective_rate * time);
          const micm::Real exact_gas = total / (1.0 + k_p);
          const micm::Real exact_aq = total * k_p / (1.0 + k_p);
          const micm::Real exact_p = a0 - total;

          const micm::Real total_value = state.variables_[0][i_total];
          const micm::Real gas = state.variables_[0][i_gas];
          const micm::Real aq = state.variables_[0][i_aq];
          const micm::Real p = state.variables_[0][i_p];

          EXPECT_NEAR(total_value, total, compare_rtol * total + compare_atol) << "A_T at step " << step;
          EXPECT_NEAR(gas, exact_gas, compare_rtol * exact_gas + compare_atol) << "A_G at step " << step;
          EXPECT_NEAR(aq, exact_aq, compare_rtol * exact_aq + compare_atol) << "A_AQ at step " << step;
          EXPECT_NEAR(p, exact_p, compare_rtol * exact_p + compare_atol) << "P_AQ at step " << step;
          EXPECT_NEAR(k_p * gas - aq, 0.0, compare_rtol * aq + compare_atol) << "Equilibrium at step " << step;
          EXPECT_NEAR(gas + aq + p, a0, compare_rtol * a0) << "Mass at step " << step;
          EXPECT_NEAR(total_value + p, a0, compare_rtol * a0) << "Total mass at step " << step;
        }
      }
    }
  }
}
