#pragma once

#include "chapman_ode_solver.hpp"
#include "util.hpp"

#include <random>

template<class SolverPolicy>
void testJacobian(SolverPolicy& solver)
{
  std::random_device rnd_device;
  std::mt19937 engine{ rnd_device() };
  std::lognormal_distribution dist{ -2.0, 4.0 };

  micm::ChapmanODESolver fixed_solver{};

  auto state = solver.GetState(3);

  auto& state_vec = state.variables_.AsVector();
  std::generate(state_vec.begin(), state_vec.end(), [&]() { return dist(engine); });
  std::generate(state_vec.begin(), state_vec.end(), [&]() { return dist(engine); });

  auto& rate_constants_vec = state.rate_constants_.AsVector();
  std::generate(rate_constants_vec.begin(), rate_constants_vec.end(), [&]() { return dist(engine); });

  auto& jacobian = state.jacobian_;
  jacobian.Fill(0.0);
  solver.solver_.rates_.SubtractJacobianTerms(state, state.variables_, jacobian);

  // The hard-coded solver stores dF/dy in a 23-element CSC array. The modern solver stores -dF/dy
  // in a 27-element sparse matrix (CSR order). The 4 extra entries in the modern solver correspond
  // to species with net-zero stoichiometry in certain reactions (M in r4, N2 in r1, O2 in r2);
  // those entries are always 0. The (row, col) pairs below map each of the 23 hard-coded entries
  // using species indices M=0, Ar=1, CO2=2, H2O=3, N2=4, O1D=5, O=6, O2=7, O3=8.
  static const std::vector<std::pair<std::size_t, std::size_t>> hardcoded_positions = {
    { 0, 0 }, { 6, 0 }, { 7, 0 }, { 8, 0 },  // col M
    { 1, 1 },                                  // col Ar
    { 2, 2 },                                  // col CO2
    { 3, 3 },                                  // col H2O
    { 4, 4 }, { 5, 4 }, { 6, 4 },              // col N2
    { 5, 5 }, { 6, 5 },                        // col O1D
    { 6, 6 }, { 7, 6 }, { 8, 6 },              // col O
    { 5, 7 }, { 6, 7 }, { 7, 7 }, { 8, 7 },   // col O2
    { 5, 8 }, { 6, 8 }, { 7, 8 }, { 8, 8 }    // col O3
  };

  for (std::size_t i{}; i < 3; ++i)
  {
    double number_density_air = 1.0;
    std::vector<double> rate_constants = state.rate_constants_[i];
    std::vector<double> variables(state.variables_.NumColumns());
    for (std::size_t j{}; j < state.variables_.NumColumns(); ++j)
      variables[j] = state.variables_[i][state.variable_map_[fixed_solver.species_names()[j]]];
    std::vector<double> fixed_jacobian = fixed_solver.dforce_dy(rate_constants, variables, number_density_air);

    EXPECT_EQ(fixed_jacobian.size(), hardcoded_positions.size());

    for (std::size_t j{}; j < fixed_jacobian.size(); ++j)
    {
      auto [hc_row, hc_col] = hardcoded_positions[j];
      std::size_t modern_row = state.variable_map_[fixed_solver.species_names()[hc_row]];
      std::size_t modern_col = state.variable_map_[fixed_solver.species_names()[hc_col]];
      // SubtractJacobianTerms stores -dF/dy; dforce_dy returns +dF/dy
      double a = jacobian[i][modern_row][modern_col];
      double b = fixed_jacobian[j];
      EXPECT_NEAR(a, -b, (std::abs(a) + std::abs(b)) * 1.0e-8 + 1.0e-12);
    }
  }
}
