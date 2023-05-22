#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>
#include <micm/solver/rosenbrock.hpp>
#include "util.hpp"

TEST(RegressionRosenbrock, rate_constants) {
  micm::ChapmanODESolver fixed_solver{};
  auto solver = getChapmanSolver();

  auto state = solver.GetState();
  auto fixed_state = fixed_solver.GetState();
  const std::vector<double> concentrations{ 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
  state.concentrations_ = concentrations;
  const std::vector<double> photo_rates{1.0e-4, 1.0e-5, 1.0e-6};
  state.custom_rate_parameters_ = photo_rates;
  state.temperature_ = 284.19;   // [K]
  state.pressure_    = 101245.0; // [Pa]
  fixed_state.temperature_ = state.temperature_;
  fixed_state.pressure_ = state.pressure_;

  solver.calculate_rate_constants(state);
  fixed_solver.calculate_rate_constants(fixed_state);
  
  EXPECT_EQ(state.rate_constants_.size(), fixed_state.rate_constants_.size());
  for(size_t i{}; i < state.rate_constants_.size(); ++i) {
    EXPECT_EQ(state.rate_constants_[i], fixed_state.rate_constants_[i]);
  }
}