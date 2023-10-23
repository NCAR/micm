#include <micm/configure/solver_config.hpp>
#include <micm/solver/rosenbrock.hpp>

using namespace micm;

void print_header()
{
  std::cout << std::setw(5) << "time"
            << "," << std::setw(10) << "A"
            << "," << std::setw(10) << "B"
            << "," << std::setw(10) << "C" << std::endl;
}

template<template<class> class T>
void print_state(double time, State<T>& state)
{
  std::ios oldState(nullptr);
  oldState.copyfmt(std::cout);

  std::cout << std::setw(5) << time << ",";
  std::cout << std::scientific << std::setprecision(2) << std::setw(10) << state.variables_[0][state.variable_map_["A"]]
            << "," << std::setw(10) << state.variables_[0][state.variable_map_["B"]] << "," << std::setw(10)
            << state.variables_[0][state.variable_map_["C"]] << std::endl;

  std::cout.copyfmt(oldState);
}

int main()
{
  constexpr size_t n_threads = 3;

  SolverConfig solverConfig;

  std::string config_path = "./configs/robertson";
  ConfigParseStatus status = solverConfig.ReadAndParse(config_path);
  if (status != micm::ConfigParseStatus::Success)
  {
    throw "Parsing failed";
  }

  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  RosenbrockSolver<> solver{ chemical_system,
                             reactions,
                             RosenbrockSolverParameters::three_stage_rosenbrock_parameters(1, false) };
  State<Matrix> state = solver.GetState();

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

  print_header();
  print_state(0, state);
  for (int i = 0; i < 10; ++i)
  {
    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_ = result.result_;
    }
    print_state(time_step * (i + 1), state);
  }
}