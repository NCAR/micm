#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <cstddef>

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

template<class SolverBuilderPolicy>
SolverBuilderPolicy getSingularSystemZeroInBottomRightOfU(SolverBuilderPolicy builder)
{
  // A -> B
  // B -> A

  auto a = micm::Species("a");
  auto b = micm::Species("b");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ a })
                         .SetProducts({ Yields(b, 1) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r1" }));

  micm::Process r2 = micm::Process::Create()
                         .SetReactants({ b })
                         .SetProducts({ Yields(a, 1) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r2" }));

  return builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
      .SetReactions(std::vector<micm::Process>{ r1, r2 })
      .SetReorderState(false);
}

template<class SolverBuilderPolicy>
SolverBuilderPolicy getSolverForSingularSystemOnDiagonal(SolverBuilderPolicy builder)
{
  // A -> B, k1
  // B -> C, k2
  // C -> A, k3

  auto a = micm::Species("a");
  auto b = micm::Species("b");
  auto c = micm::Species("c");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c } };

  micm::Process r1 = micm::Process::Create()
                         .SetReactants({ a })
                         .SetProducts({ Yields(b, 1) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r1" }));

  micm::Process r2 = micm::Process::Create()
                         .SetReactants({ b })
                         .SetProducts({ Yields(c, 1) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r2" }));

  micm::Process r3 = micm::Process::Create()
                         .SetReactants({ c })
                         .SetProducts({ Yields(a, 1) })
                         .SetPhase(gas_phase)
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r3" }));

  return builder.SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
      .SetReactions(std::vector<micm::Process>{ r1, r2, r3 })
      .SetReorderState(false);
}