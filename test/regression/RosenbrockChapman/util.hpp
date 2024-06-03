#pragma once

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <gtest/gtest.h>

#include <utility>
#include <vector>

using yields = std::pair<micm::Species, double>;

micm::Phase createGasPhase()
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

  return micm::Phase{ std::vector<micm::Species>{ m, ar, co2, h2o, n2, o1d, o, o2, o3 } };
}

std::vector<micm::Process> createProcesses(const micm::Phase& gas_phase)
{
  micm::Process r1 =
      micm::Process::Create()
          .SetReactants({ micm::Species("O1D"), micm::Species("N2") })
          .SetProducts({ Yields(micm::Species("O"), 1), Yields(micm::Species("N2"), 1) })
          .SetRateConstant(micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2.15e-11, .C_ = 110 }))
          .SetPhase(gas_phase);

  micm::Process r2 =
      micm::Process::Create()
          .SetReactants({ micm::Species("O1D"), micm::Species("O2") })
          .SetProducts({ Yields(micm::Species("O"), 1), Yields(micm::Species("O2"), 1) })
          .SetRateConstant(micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .C_ = 55 }))
          .SetPhase(gas_phase);

  micm::Process r3 =
      micm::Process::Create()
          .SetReactants({ micm::Species("O"), micm::Species("O3") })
          .SetProducts({ Yields(micm::Species("O2"), 2) })
          .SetRateConstant(micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 8e-12, .C_ = -2060 }))
          .SetPhase(gas_phase);

  micm::Process r4 =
      micm::Process::Create()
          .SetReactants({ micm::Species("O"), micm::Species("O2"), micm::Species("M") })
          .SetProducts({ Yields(micm::Species("O3"), 1), Yields(micm::Species("M"), 1) })
          .SetRateConstant(micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .B_ = 2.4 }))
          .SetPhase(gas_phase);

  micm::Process photo_1 = micm::Process::Create()
                              .SetReactants({ micm::Species("O2") })
                              .SetProducts({ Yields(micm::Species("O"), 2) })
                              .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "jO2" }))
                              .SetPhase(gas_phase);

  micm::Process photo_2 = micm::Process::Create()
                              .SetReactants({ micm::Species("O3") })
                              .SetProducts({ Yields(micm::Species("O1D"), 1), Yields(micm::Species("O2"), 1) })
                              .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "jO3a" }))
                              .SetPhase(gas_phase);

  micm::Process photo_3 = micm::Process::Create()
                              .SetReactants({ micm::Species("O3") })
                              .SetProducts({ Yields(micm::Species("O"), 1), Yields(micm::Species("O2"), 1) })
                              .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "jO3b" }))
                              .SetPhase(gas_phase);

  return { photo_1, photo_2, photo_3, r1, r2, r3, r4 };
}

template<class SolverBuilderPolicy>
auto getChapmanSolver(
    SolverBuilderPolicy& builder,
    const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  return builder
        .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
        .SetReactions(std::move(processes))
        .SetNumberOfGridCells(number_of_grid_cells)
        .Build();
}