#include <gtest/gtest.h>

#include <micm/solver/chapman_ode_solver.hpp>
#include <micm/solver/rosenbrock.hpp>

#include "util.hpp"

TEST(RegressionRosenbrock, rate_constants)
{
  micm::ChapmanODESolver fixed_solver{};
  auto solver = getMultiCellChapmanSolver(3);

  auto state = solver.GetState();
  auto fixed_state = fixed_solver.GetState();
  const std::vector<std::vector<double>> photo_rates{
    { 1.0e-4, 1.0e-5, 1.0e-6 },
    { 3.2e-4, 7.3e-5, 3.2e-6 },
    { 5.2e-4, 8.2e-5, 4.6e-6 }
  };
  state.custom_rate_parameters_ = photo_rates;
  state.conditions_[0].temperature_ = 284.19;    // [K]
  state.conditions_[0].pressure_    = 101245.0;  // [Pa]
  state.conditions_[1].temperature_ = 215.02;    // [K]
  state.conditions_[1].pressure_    = 100789.2;  // [Pa]
  state.conditions_[2].temperature_ = 299.31;    // [K]
  state.conditions_[2].pressure_    = 101398.0;  // [Pa]

  solver.UpdateState(state);
  
  for (size_t i{}; i < 3; ++i) {
    fixed_state.conditions_[0].temperature_ = state.conditions_[i].temperature_;
    fixed_state.conditions_[0].pressure_ = state.conditions_[i].pressure_;
    fixed_state.custom_rate_parameters_[0] = photo_rates[i];
    fixed_solver.UpdateState(fixed_state);
  
    EXPECT_EQ(state.rate_constants_[i].size(), fixed_state.rate_constants_[0].size());
    for (size_t j{}; j < state.rate_constants_[i].size(); ++j)
    {
      EXPECT_EQ(state.rate_constants_[i][j], fixed_state.rate_constants_[0][j]);
    }
  }
}