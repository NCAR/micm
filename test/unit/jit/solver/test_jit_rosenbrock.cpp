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
using JitBuilder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>;

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

TEST(JitRosenbrockSolver, SingularSystemZeroInBottomRightOfU)
{
  auto params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  params.check_singularity_ = true;
  auto vector = JitBuilder<4>(params);

  auto vector_solver = getSingularSystemZeroInBottomRightOfU(vector).SetNumberOfGridCells(4).Build();

  auto vector_state = vector_solver.GetState();

  double k1 = -2;
  double k2 = 1.0;

  vector_state.SetCustomRateParameter("r1", { k1, k1, k1, k1 });
  vector_state.SetCustomRateParameter("r2", { k2, k2, k2, k2 });

  vector_state.variables_[0] = { 1.0, 1.0 };
  vector_state.variables_[1] = { 1.0, 1.0 };
  vector_state.variables_[2] = { 1.0, 1.0 };
  vector_state.variables_[3] = { 1.0, 1.0 };

  // to get a jacobian with an LU factorization that contains a zero on the diagonal
  // of U, we need det(alpha * I - jacobian) = 0
  // for the system above, that means we have to have alpha + k1 + k2 = 0
  // in this case, one of the reaction rates will be negative but it's good enough to
  // test the singularity check
  // alpha is 1 / (H * gamma), where H is the time step and gamma is the gamma value from
  // the rosenbrock paramters
  // so H needs to be 1 / ( (-k1 - k2) * gamma)
  // since H is positive we need -k1 -k2 to be positive, hence the smaller, negative value for k1
  double H = 1 / ((-k1 - k2) * params.gamma_[0]);
  params.h_start_ = H;

  vector_solver.CalculateRateConstants(vector_state);

  auto vector_result = vector_solver.Solve(2 * H, vector_state, params);
  EXPECT_NE(vector_result.stats_.singular_, 0);
}

TEST(JitRosenbrockSolver, SingularSystemZeroAlongDiagonalNotBottomRight)
{
  auto params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  double k1 = -1.0;
  double k2 = -1.0;
  double k3 = 1.0;

  // to get a jacobian with an LU factorization that contains a zero on the diagonal
  // of U, we need det(alpha * I - jacobian) = 0
  // for the system above, that means we have to set alpha = -k1, or alpha=-k2, or alpha=k3
  double H = 1 / (-k1 * params.gamma_[0]);

  params.check_singularity_ = true;
  params.h_start_ = H;

  auto vector = JitBuilder<4>(params);

  auto vector_solver = getSolverForSingularSystemOnDiagonal(vector).SetNumberOfGridCells(4).Build();

  auto vector_state = vector_solver.GetState();

  vector_state.SetCustomRateParameter("r1", { k1, k1, k1, k1 });
  vector_state.SetCustomRateParameter("r2", { k2, k2, k2, k2 });
  vector_state.SetCustomRateParameter("r3", { k3, k3, k3, k3 });

  vector_state.variables_[0] = { 1.0, 1.0, 1.0 };
  vector_state.variables_[1] = { 1.0, 1.0, 1.0 };
  vector_state.variables_[2] = { 1.0, 1.0, 1.0 };
  vector_state.variables_[3] = { 1.0, 1.0, 1.0 };

  vector_solver.CalculateRateConstants(vector_state);

  auto vector_result = vector_solver.Solve(2 * H, vector_state);
  EXPECT_NE(vector_result.stats_.singular_, 0);
}