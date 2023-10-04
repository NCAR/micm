#include <gtest/gtest.h>
#include <omp.h>

#include <micm/configure/solver_config.hpp>
#include <micm/solver/rosenbrock.hpp>

using namespace micm;

template<class T>
using SparseMatrixPolicy = SparseMatrix<T>;

std::vector<double> test_solver_on_thread(System chemical_system, std::vector<Process> reactions)
{
  std::cout << "Running solver on thread " << omp_get_thread_num() << std::endl;
  RosenbrockSolver<> solver{ chemical_system,
                                  reactions,
                                  RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };
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

TEST(OpenMP, OneFileReadThreeThreads)
{
  constexpr size_t n_threads = 3;

  SolverConfig solverConfig;

  std::string config_path = "./unit_configs/robertson";
  ConfigParseStatus status = solverConfig.ReadAndParse(config_path);
  if (status != micm::ConfigParseStatus::Success)
  {
    throw "Parsing failed";
  }

  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  std::vector<std::vector<double>> results(3);

#pragma omp parallel num_threads(n_threads)
  {
    std::vector<double> result = test_solver_on_thread(chemical_system, reactions);
    results[omp_get_thread_num()] = result;
#pragma omp barrier 
  }

  for(int i = 0; i < results[0].size(); ++i) {
    EXPECT_EQ(results[0][i], results[1][i]);
    EXPECT_EQ(results[0][i], results[2][i]);
    EXPECT_EQ(results[1][i], results[2][i]);
  }
}