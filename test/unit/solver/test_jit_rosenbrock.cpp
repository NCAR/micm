#include <micm/configure/solver_config.hpp>
#include <micm/jit/jit_compiler.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/jit_rosenbrock.hpp>
#include <micm/solver/jit_solver_parameters.hpp>
#include <micm/solver/jit_solver_builder.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

template<class SolverBuilderPolicy>
SolverBuilderPolicy getSolver(SolverBuilderPolicy builder)
{
  // ---- foo  bar  baz  quz  quuz
  // foo   0    1    2    -    -
  // bar   3    4    5    -    -
  // baz   6    -    7    -    -
  // quz   -    8    -    9    -
  // quuz 10    -   11    -    12

  auto foo = micm::Species("foo");
  auto bar = micm::Species("bar");
  auto baz = micm::Species("baz");
  auto quz = micm::Species("quz");
  auto quuz = micm::Species("quuz");

  micm::Phase gas_phase{ std::vector<micm::Species>{ foo, bar, baz, quz, quuz } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ foo, baz })
                         .SetProducts({ Yields(bar, 1), Yields(quuz, 2.4) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.0e-11, .B_ = 0, .C_ = 110 }));

  micm::Process r2 = micm::Process::Create()
                         .SetReactants({ bar })
                         .SetProducts({ Yields(foo, 1), Yields(quz, 1.4) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 1.0e-6 }));

  micm::Process r3 = micm::Process::Create().SetReactants({ quz }).SetProducts({}).SetPhase(gas_phase).SetRateConstant(
      micm::ArrheniusRateConstant({ .A_ = 3.5e-6 }));

  return builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
      .SetReactions(std::vector<micm::Process>{ r1, r2, r3 })
      .SetReorderState(false);
}

template<class SolverBuilderPolicy>
void testAlphaMinusJacobian(SolverBuilderPolicy builder, std::size_t number_of_grid_cells)
{
  builder = getSolver(builder);
  auto solver = builder.SetNumberOfGridCells(number_of_grid_cells).Build();
  auto jacobian = solver.GetState().jacobian_;

  EXPECT_EQ(jacobian.NumberOfBlocks(), number_of_grid_cells);
  EXPECT_EQ(jacobian.NumRows(), 5);
  EXPECT_EQ(jacobian.NumColumns(), jacobian.NumRows());
  EXPECT_EQ(jacobian[0].Size(), 5);
  EXPECT_EQ(jacobian[0][0].Size(), 5);
  EXPECT_GE(jacobian.AsVector().size(), 13 * number_of_grid_cells);

  // Generate a negative Jacobian matrix
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    jacobian[i_cell][0][0] = -12.2;
    jacobian[i_cell][0][1] = -24.3 * (i_cell + 2);
    jacobian[i_cell][0][2] = -42.3;
    jacobian[i_cell][1][0] = -0.43;
    jacobian[i_cell][1][1] = -23.4;
    jacobian[i_cell][1][2] = -83.4 / (i_cell + 3);
    jacobian[i_cell][2][0] = -4.74;
    jacobian[i_cell][2][2] = -6.91;
    jacobian[i_cell][3][1] = -59.1;
    jacobian[i_cell][3][3] = -83.4;
    jacobian[i_cell][4][0] = -78.5;
    jacobian[i_cell][4][2] = -53.6;
    jacobian[i_cell][4][4] = -1.0;
  }
  solver.solver_.AlphaMinusJacobian(jacobian, 42.042);
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    EXPECT_NEAR(jacobian[i_cell][0][0], 42.042 - 12.2, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][0][1], -24.3 * (i_cell + 2), 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][0][2], -42.3, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][1][0], -0.43, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][1][1], 42.042 - 23.4, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][1][2], -83.4 / (i_cell + 3), 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][2][0], -4.74, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][2][2], 42.042 - 6.91, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][3][1], -59.1, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][3][3], 42.042 - 83.4, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][4][0], -78.5, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][4][2], -53.6, 1.0e-5);
    EXPECT_NEAR(jacobian[i_cell][4][4], 42.042 - 1.0, 1.0e-5);
  }
}

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

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
    }
  }
}

template <std::size_t L>
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
