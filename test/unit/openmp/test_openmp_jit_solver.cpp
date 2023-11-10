#include <gtest/gtest.h>
#include <omp.h>

#include <micm/configure/solver_config.hpp>
#include <micm/jit/jit_compiler.hpp>
#include <micm/solver/jit_rosenbrock.hpp>
#include <micm/util/sparse_matrix.hpp>

#include "run_solver.hpp"

using namespace micm;

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;

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
      JitLinearSolver<1, Group1SparseVectorMatrix>,
      micm::JitProcessSet<1>>
      solver(jit.get(), chemical_system, reactions, RosenbrockSolverParameters::three_stage_rosenbrock_parameters());

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
