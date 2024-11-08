#include "../../solver/test_rosenbrock_solver_policy.hpp"

#include <micm/configure/solver_config.hpp>
#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/solver/jit_rosenbrock.hpp>
#include <micm/jit/solver/jit_solver_builder.hpp>
#include <micm/jit/solver/jit_solver_parameters.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

void run_solver(auto& solver)
{
  auto state = solver.GetState();

  state.variables_[0] = { 1, 0, 0 };

  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa
  state.conditions_[0].air_density_ = 1e6;     // mol m-3

  double time_step = 500;  // s

  for (int i = 0; i < 10; ++i)
  {
    double elapsed_solve_time = 0;
    solver.CalculateRateConstants(state);

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
    }
  }
}

template<std::size_t L>
using JitBuilder = micm::JitSolverBuilder_Doolittle<micm::JitRosenbrockSolverParameters, L>;

TEST(JitRosenbrockSolver, AlphaMinusJacobian)
{
  testAlphaMinusJacobian(JitBuilder<1>(micm::JitRosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 1);
  testAlphaMinusJacobian(JitBuilder<2>(micm::JitRosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 2);
  testAlphaMinusJacobian(JitBuilder<3>(micm::JitRosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 3);
  testAlphaMinusJacobian(JitBuilder<4>(micm::JitRosenbrockSolverParameters::ThreeStageRosenbrockParameters()), 4);
}

TEST(JitRosenbrockSolver, MultipleInstances)
{
  micm::SolverConfig solverConfig;
  std::string config_path = "./unit_configs/robertson";
  solverConfig.ReadAndParse(config_path);

  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  auto builder = JitBuilder<1>(micm::JitRosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                     .SetSystem(chemical_system)
                     .SetReactions(reactions);

  auto solver1 = builder.Build();
  auto solver2 = builder.Build();
  auto solver3 = builder.Build();

  run_solver(solver1);
  run_solver(solver2);
  run_solver(solver3);
}