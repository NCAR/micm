#include <gtest/gtest.h>
#include <omp.h>

#include <micm/configure/solver_config.hpp>
#include <micm/jit/jit_compiler.hpp>
#include <micm/solver/jit_rosenbrock.hpp>
#include <micm/util/sparse_matrix.hpp>

using namespace micm;

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;

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

TEST(OpenMP, JITOneSolverManyStates)
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

  auto jit{ micm::JitCompiler::create() };
  JitRosenbrockSolver<
    Group1VectorMatrix, 
    Group1SparseVectorMatrix, 
    JitLinearSolver<1, Group1SparseVectorMatrix>
    > solver(
      jit.get(), 
      chemical_system, 
      reactions, 
      RosenbrockSolverParameters::three_stage_rosenbrock_parameters() 
    );

#pragma omp parallel num_threads(n_threads)
  {
    auto state = solver.GetState();
    std::vector<double> result = run_solver_on_thread_with_own_state(solver, state);
    results[omp_get_thread_num()] = result;
#pragma omp barrier
  }

  // compare each thread to thread 1
  for (int i = 1; i < n_threads; ++i) {
    for (int j = 0; j < results[0].size(); ++j) {
      EXPECT_EQ(results[0][j], results[i][j]);
    }
  }
}
