// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// The 7-reaction Chapman mechanism, as a benchmark mechanism.
//
// This is the same system as test/integration/test_chapman_integration.cpp,
// scaled to many grid cells so per-Solve() work dominates over per-call
// overhead.
#pragma once

#include <micm/CPU.hpp>

#include <string_view>
#include <vector>

namespace bench
{
  struct Chapman
  {
    static constexpr std::string_view kName = "chapman";

    template<class Builder>
    static auto Build(Builder builder)
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

      return builder.SetSystem(micm::System(gas_phase))
          .SetReactions({ r1, r2, r3, r4, photo_1, photo_2, photo_3 })
          .SetIgnoreUnusedSpecies(true)
          .Build();
    }

    template<class State>
    static void InitState(State& state, micm::Index num_cells)
    {
      // Initial concentrations: matches the integration test's cell 0 for consistency.
      std::vector<micm::Real> concentrations{ 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
      std::vector<micm::Real> photo_rates{ 0.1, 0.2, 0.3 };
      for (micm::Index c = 0; c < num_cells; ++c)
      {
        state.variables_[c] = concentrations;
        state.custom_rate_parameters_[c] = photo_rates;
        state.conditions_[c].temperature_ = 273.15 + 25.0;
        state.conditions_[c].pressure_ = 101325.0;
      }
    }
  };
}  // namespace bench
