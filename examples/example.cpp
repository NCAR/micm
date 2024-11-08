// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "util_example.hpp"

#include <micm/configure/solver_config.hpp>
#include <micm/solver/solver_builder.hpp>

#include <iostream>
#include <string>

namespace fs = std::filesystem;
using namespace micm;

int main(const int argc, const char* argv[])
{
  if (argc < 3)
  {
    std::cout << "  path                  Path to the folder holding the configuration data" << std::endl;
    std::cout << "  initial condition     csv file contanining initial condtions" << std::endl;
    return 1;
  }

  fs::path config_path{ argv[1] };
  std::string initial_conditions_file{ argv[2] };

  // Read csv file containing initial conditions
  auto dataMap = readCSVToMap(initial_conditions_file);
  if (dataMap.empty())
  {
    std::cout << "Error: Failed in reading CSV file. Please verify file path and/or file name of \""
              << initial_conditions_file << "\"" << std::endl;
    return 1;
  }
  else if (dataMap.is_incomplete())
  {
    std::cout << "Error: Parsing csv file is incomplete . Please verify content formatting of \"" << initial_conditions_file
              << "\"" << std::endl;
    return 1;
  }

  // Read configure files
  SolverConfig solverConfig;
  solverConfig.ReadAndParse(config_path);

  SolverParameters solver_params = solverConfig.GetSolverParams();

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  auto solver = CpuSolverBuilder_DoolittleLU<micm::RosenbrockSolverParameters>(solver_params.parameters_)
                   .SetSystem(chemical_system)
                   .SetReactions(reactions)
                   .Build();

  State state = solver.GetState();

  state.conditions_[0].temperature_ = dataMap.environments["temperature"];  // K
  state.conditions_[0].pressure_ = dataMap.environments["pressure"];        // Pa
  state.conditions_[0].air_density_ = dataMap.environments["air_density"];  // mol m-3

  state.SetConcentrations(dataMap.concentrations);

  if (!dataMap.custom_rate_params.empty())
  {
    for (auto& [rate_constant_key, rate_constant_val] : dataMap.custom_rate_params)
    {
      state.SetCustomRateParameter(rate_constant_key, rate_constant_val[0]);
    }
  }

  // choose a timestep and print the initial state
  double time_step = dataMap.environments["time_step"];  // s

  state.PrintHeader();
  state.PrintState(0);

  // Depending on how stiff the system is
  // the solver integration step may not be able to solve for the full time step
  // so we need to track how much time the solver was able to integrate for and continue
  // solving until we finish
  double elapsed_solve_time = 0;

  while (elapsed_solve_time < time_step)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(time_step - elapsed_solve_time, state);
    elapsed_solve_time = result.final_time_;
    if (result.state_ != SolverState::Converged)
    {
      std::cout << "solver failed to converge" << std::endl;
      return 1;
    }
  }
  state.PrintState(time_step);

  return 0;
}