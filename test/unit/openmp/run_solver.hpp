#pragma once

std::vector<double> run_solver_on_thread_with_own_state(auto& solver, auto& state)
{
  std::cout << "Running solver on thread " << omp_get_thread_num() << std::endl;

  // mol m-3
  state.variables_[0] = { 1, 0, 0 };

  double k1 = 0.04;
  double k2 = 3e7;
  double k3 = 1e4;
  state.SetCustomRateParameter("PHOTO.r1", k1);
  state.SetCustomRateParameter("PHOTO.r2", k2);
  state.SetCustomRateParameter("PHOTO.r3", k3);

  double temperature = 272.5;  // [K]
  double pressure = 101253.3;  // [Pa]
  double air_density = 1e6;    // [mol m-3]

  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  double time_step = 200;  // s

  for (int i = 0; i < 10; ++i)
  {
    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_ = result.result_;
    }
  }

  return state.variables_.AsVector();
}