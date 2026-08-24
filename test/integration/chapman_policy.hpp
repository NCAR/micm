// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "../precision_matchers.hpp"

#include <micm/CPU.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

/// @brief Integrate the Chapman ozone mechanism and check the invariants it must obey.
///
/// The seven reactions built below all conserve oxygen atoms:
///
///   O1D + N2   -> O + N2      (1 -> 1)
///   O1D + O2   -> O + O2      (3 -> 3)
///   O   + O3   -> 2 O2        (4 -> 4)
///   O   + O2 + M -> O3 + M    (3 -> 3)
///   O2         -> 2 O         (2 -> 2)
///   O3         -> O1D + O2    (3 -> 3)
///   O3         -> O + O2      (3 -> 3)
///
/// so N_O = [O] + [O1D] + 2*[O2] + 3*[O3] is invariant. Every other species appears with
/// identical stoichiometry on both sides (N2 and M) or in no reaction at all (Ar, H2O, CO2),
/// so those concentrations must not move either.
///
/// Species indices come from state.variable_map_ rather than from the declaration order,
/// because the solver is free to reorder the state.
template<class BuilderPolicy>
void TestChapmanSystem(BuilderPolicy builder)
{
  auto o = micm::Species("O");
  auto o1d = micm::Species("O1D");
  auto o2 = micm::Species("O2");
  auto o3 = micm::Species("O3");
  auto m = micm::Species("M");
  auto ar = micm::Species("Ar");
  auto n2 = micm::Species("N2");
  auto h2o = micm::Species("H2O");
  auto co2 = micm::Species("CO2");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ o, o1d, o2, o3, m, ar, n2, h2o, co2 } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ o1d, n2 })
                         .SetProducts({ micm::StoichSpecies(o, 1), micm::StoichSpecies(n2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 })
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ o1d, o2 })
                         .SetProducts({ micm::StoichSpecies(o, 1), micm::StoichSpecies(o2, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 })
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ o, o3 })
                         .SetProducts({ micm::StoichSpecies(o2, 2) })
                         .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 8e-12, .B_ = 0, .C_ = -2060 })
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r4 = micm::ChemicalReactionBuilder()
                         .SetReactants({ o, o2, m })
                         .SetProducts({ micm::StoichSpecies(o3, 1), micm::StoichSpecies(m, 1) })
                         .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .B_ = 0, .C_ = 2.4 })
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process photo_1 = micm::ChemicalReactionBuilder()
                              .SetReactants({ o2 })
                              .SetProducts({ micm::StoichSpecies(o, 2) })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jO2" })
                              .SetPhase(gas_phase)
                              .Build();

  micm::Process photo_2 = micm::ChemicalReactionBuilder()
                              .SetReactants({ o3 })
                              .SetProducts({ micm::StoichSpecies(o1d, 1), micm::StoichSpecies(o2, 1) })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jO3a" })
                              .SetPhase(gas_phase)
                              .Build();

  micm::Process photo_3 = micm::ChemicalReactionBuilder()
                              .SetReactants({ o3 })
                              .SetProducts({ micm::StoichSpecies(o, 1), micm::StoichSpecies(o2, 1) })
                              .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jO3b" })
                              .SetPhase(gas_phase)
                              .Build();

  auto solver = builder.SetSystem(micm::System(gas_phase))
                    .SetReactions({ r1, r2, r3, r4, photo_1, photo_2, photo_3 })
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  auto state = solver.GetState();

  std::vector<micm::Real> concentrations{ 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
  state.variables_[0] = concentrations;
  std::vector<micm::Real> photo_rates{ 0.1, 0.2, 0.3 };
  state.custom_rate_parameters_[0] = photo_rates;
  state.conditions_[0].temperature_ = 2;
  state.conditions_[0].pressure_ = 3;

  const auto& variable_map = state.variable_map_;
  // Oxygen-atom count, weighting each species by the oxygen atoms it carries.
  auto oxygen_atoms = [&variable_map](const auto& state_ref)
  {
    micm::Real total = 0.0;
    for (const auto& [name, weight] : { std::pair<std::string, micm::Real>{ "O", 1.0 },
                                        std::pair<std::string, micm::Real>{ "O1D", 1.0 },
                                        std::pair<std::string, micm::Real>{ "O2", 2.0 },
                                        std::pair<std::string, micm::Real>{ "O3", 3.0 } })
    {
      total += weight * state_ref.variables_[0][variable_map.at(name)];
    }
    return total;
  };

  const micm::Real oxygen_initial = oxygen_atoms(state);
  const micm::Real o_initial = state.variables_[0][variable_map.at("O")];
  const micm::Real o2_initial = state.variables_[0][variable_map.at("O2")];

  // Species that no reaction can change: identical stoichiometry on both sides, or absent
  // from every reaction. Only assert on the ones the solver actually kept in the state.
  std::vector<std::pair<std::string, micm::Real>> spectators;
  for (const auto& name : { "M", "Ar", "N2", "H2O", "CO2" })
  {
    auto it = variable_map.find(name);
    if (it != variable_map.end())
    {
      spectators.emplace_back(name, state.variables_[0][it->second]);
    }
  }

  state.variables_.CopyToDevice();
  state.custom_rate_parameters_.CopyToDevice();
  state.conditions_.CopyToDevice();

  for (micm::Index t{}; t < 100; ++t)
  {
    state.custom_rate_parameters_[0] = photo_rates;
    state.custom_rate_parameters_.CopyToDevice();
    solver.UpdateStateParameters(state);
    auto result = solver.Solve(30.0, state);
    ASSERT_EQ(result.state_, micm::SolverState::Converged)
        << "step " << t << " did not converge: " << micm::SolverStateToString(result.state_);

    state.variables_.CopyToHost();

    for (const auto& name : state.variable_names_)
    {
      const micm::Real value = state.variables_[0][variable_map.at(name)];
      ASSERT_TRUE(std::isfinite(static_cast<double>(value))) << name << " not finite at step " << t;
      EXPECT_GE(value, micm::Real{ 0.0 }) << name << " went negative at step " << t;
    }

    EXPECT_REAL_REL(oxygen_atoms(state), oxygen_initial, 1.0e-8) << "oxygen atoms not conserved at step " << t;

    for (const auto& [name, expected] : spectators)
    {
      EXPECT_REAL_REL(state.variables_[0][variable_map.at(name)], expected, 1.0e-8)
          << name << " participates in no net reaction but changed at step " << t;
    }
  }

  // The photolysis rates are user-defined inputs, so a state that never reached the device
  // would sit at its initial values. O2 is destroyed by jO2 = 0.1 /s over 100 steps of 30 s;
  // essentially all of the oxygen ends up as atomic O.
  const micm::Real o_final = state.variables_[0][variable_map.at("O")];
  const micm::Real o2_final = state.variables_[0][variable_map.at("O2")];
  EXPECT_LT(o2_final, micm::Real{ 0.5 } * o2_initial) << "O2 was not photolysed -- did the rates reach the device?";
  EXPECT_GT(o_final, o_initial) << "O was not produced";
  EXPECT_GT(o_final, micm::Real{ 0.5 } * oxygen_initial) << "oxygen did not end up as atomic O";
}
