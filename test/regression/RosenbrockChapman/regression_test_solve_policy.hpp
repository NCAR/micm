#include "chapman_ode_solver.hpp"
#include "util.hpp"

#include <random>

template<class SolverPolicy>
void testSolve(SolverPolicy& solver, double relative_tolerance = 1.0e-8)
{
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());
  micm::ChapmanODESolver fixed_solver{};

  auto state = solver.GetState();
  auto fixed_state = fixed_solver.GetState();

  // set conditions
  const std::vector<std::vector<double>> photo_rates{ { 1.0e-4, 1.0e-5, 1.0e-6 },
                                                      { 3.2e-4, 7.3e-5, 3.2e-6 },
                                                      { 5.2e-4, 8.2e-5, 4.6e-6 } };
  state.custom_rate_parameters_ = photo_rates;
  state.conditions_[0].temperature_ = 284.19;  // [K]
  state.conditions_[0].pressure_ = 101245.0;   // [Pa]
  state.conditions_[1].temperature_ = 215.02;  // [K]
  state.conditions_[1].pressure_ = 100789.2;   // [Pa]
  state.conditions_[2].temperature_ = 299.31;  // [K]
  state.conditions_[2].pressure_ = 101398.0;   // [Pa]
  std::vector<std::vector<double>> variables(3, std::vector<double>(fixed_state.variables_.NumColumns(), 0.0));
  double abs_tol = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < fixed_state.variables_.NumColumns(); ++j)
    {
      variables[i][j] = get_double();
      abs_tol = std::max(abs_tol, variables[i][j]);
    }
  state.variables_ = variables;
  abs_tol *= 1.0e-12;

  // run solvers
  solver.CalculateRateConstants(state);
  auto results = solver.Solve(500.0, state);
  micm::ChapmanODESolver::SolverResult fixed_results[3];
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
      fixed_state.custom_rate_parameters_[0][j] = photo_rates[i][j];
    fixed_state.conditions_[0].temperature_ = state.conditions_[i].temperature_;
    fixed_state.conditions_[0].pressure_ = state.conditions_[i].pressure_;
    for (int j = 0; j < fixed_state.variables_.NumColumns(); ++j)
      fixed_state.variables_[0][j] = variables[i][state.variable_map_[fixed_solver.species_names()[j]]];
    fixed_solver.UpdateState(fixed_state);
    fixed_results[i] = fixed_solver.Solve(0.0, 500.0, fixed_state);
  }

  // compare results
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < fixed_results[i].result_.size(); ++j)
    {
      double a = state.variables_[i][state.variable_map_[fixed_solver.species_names()[j]]];
      double b = fixed_results[i].result_[j];
      EXPECT_NEAR(a, b, (std::abs(a) + std::abs(b)) * relative_tolerance + abs_tol);
    }
}