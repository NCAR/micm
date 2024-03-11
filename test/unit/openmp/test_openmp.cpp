#include <gtest/gtest.h>
#include <omp.h>

#include <micm/configure/solver_config.hpp>
#include <micm/solver/rosenbrock.hpp>

#include "run_solver.hpp"

using namespace micm;

TEST(OpenMP, OneSolverManyStates)
{
  constexpr size_t n_threads = 8;

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

  std::vector<std::vector<double>> results(n_threads);

  RosenbrockSolver<> solver{ chemical_system, reactions, RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };

#pragma omp parallel num_threads(n_threads)
  {
    auto state = solver.GetState();
    std::vector<double> result = run_solver_on_thread_with_own_state(solver, state);
    results[omp_get_thread_num()] = result;
#pragma omp barrier
  }

  // compare each thread to thread 1
  for (int i = 1; i < n_threads; ++i)
  {
    for (int j = 0; j < results[0].size(); ++j)
    {
      EXPECT_EQ(results[0][j], results[i][j]);
    }
  }
}
