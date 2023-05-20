#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>
#include <micm/solver/rosenbrock.hpp>
#include "util.hpp"

TEST(RegressionRosenbrock, rate_constants) {
  micm::ChapmanODESolver fixed_solver{};
  auto solver = getChapmanSolver();

  auto state = solver.GetState();
  const std::vector<double> concentrations{ 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
  state.concentrations_ = concentrations;
  const std::vector<double> photo_rates{1.0e-4, 1.0e-5, 1.0e-6};
  state.custom_rate_parameters_ = photo_rates;
  state.temperature_ = 284.19;   // [K]
  state.pressure_    = 101245.0; // [Pa]

  solver.calculate_rate_constants(state);
  fixed_solver.calculate_rate_constants(state);
  
  EXPECT_EQ(solver.rate_constants_.size(), fixed_solver.rate_constants_.size());
  for(size_t i{}; i < solver.rate_constants_.size(); ++i) {
    EXPECT_EQ(solver.rate_constants_[i], fixed_solver.rate_constants_[i]);
  }
}