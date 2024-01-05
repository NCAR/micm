#include <omp.h>

#include <micm/configure/solver_config.hpp>
#include <micm/solver/rosenbrock.hpp>

using namespace micm;

void print_header()
{
  std::cout << std::setw(10) << "A"
            << "," << std::setw(10) << "B"
            << "," << std::setw(10) << "C" << std::endl;
}

void print_results(std::vector<double> results)
{
  std::ios oldState(nullptr);
  oldState.copyfmt(std::cout);

  std::cout << std::scientific << std::setprecision(2) << std::setw(10) << results[0] << "," << std::setw(10) << results[1]
            << "," << std::setw(10) << results[2] << std::endl;

  std::cout.copyfmt(oldState);
}

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

  std::vector<std::vector<double>> results(n_threads);

  RosenbrockSolver<> solver{ chemical_system, reactions, RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };

#pragma omp parallel num_threads(n_threads)
  {
    // each thread should use its own state
    auto state = solver.GetState();
    std::vector<double> result = run_solver_on_thread_with_own_state(solver, state);
    results[omp_get_thread_num()] = result;
#pragma omp barrier
  }

  std::cout << "Thread 1" << std::endl;
  print_header();
  print_results(results[0]);

  std::cout << "Thread 2" << std::endl;
  print_header();
  print_results(results[1]);

  std::cout << "Thread 3" << std::endl;
  print_header();
  print_results(results[2]);
}