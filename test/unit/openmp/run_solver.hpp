#pragma once

#include <micm/util/types.hpp>

#include <omp.h>

#include <iostream>
#include <vector>

std::vector<micm::Real> run_solver_on_thread_with_own_state(auto& solver, auto& state)
{
  std::cout << "Running solver on thread " << omp_get_thread_num() << std::endl;

  // mol m-3
  state.variables_[0] = { 1, 0, 0 };

  micm::Real k1 = 0.04;
  micm::Real k2 = 3e7;
  micm::Real k3 = 1e4;
  state.SetCustomRateParameter("r1", k1);
  state.SetCustomRateParameter("r2", k2);
  state.SetCustomRateParameter("r3", k3);

  micm::Real temperature = 272.5;  // [K]
  micm::Real pressure = 101253.3;  // [Pa]
  micm::Real air_density = 1e6;    // [mol m-3]

  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  micm::Real time_step = 200;  // s

  for (micm::Index i = 0; i < 10; ++i)
  {
    micm::Real elapsed_solve_time = 0;
    solver.UpdateStateParameters(state);

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.stats_.final_time_;
      ;
    }
  }

  return state.variables_.AsVector();
}