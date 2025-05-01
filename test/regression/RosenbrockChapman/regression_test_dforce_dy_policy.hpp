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
  auto& rate_const_vec = state.rate_constants_.AsVector();
  std::generate(state_vec.begin(), state_vec.end(), [&]() { return dist(engine); });

  auto& jacobian = state.jacobian_;
  jacobian.Fill(0.0);
  solver.solver_.rates_.SubtractJacobianTerms(state.rate_constants_, state.variables_, jacobian);

  for (std::size_t i{}; i < 3; ++i)
  {
    double number_density_air = 1.0;
    std::vector<double> rate_constants = state.rate_constants_[i];
    std::vector<double> variables(state.variables_.NumColumns());
    for (std::size_t j{}; j < state.variables_.NumColumns(); ++j)
      variables[j] = state.variables_[i][state.variable_map_[fixed_solver.species_names()[j]]];
    std::vector<double> fixed_jacobian = fixed_solver.dforce_dy(rate_constants, variables, number_density_air);

    // TODO: The sparse matrix data ordering in the hard-coded solver is different (maybe because of pivoting?)
    //       As the remaining linear solver functions are generalized, use the logic in the preprocessor to
    //       decipher the data elements in the hard-coded solver sparse matrix to finish this test.
    // EXPECT_EQ(jacobian.FlatBlockSize(), fixed_jacobian.size());
    for (std::size_t j{}; j < fixed_jacobian.size(); ++j)
    {
      // EXPECT_NEAR(jacobian.AsVector()[i * jacobian.FlatBlockSize() + j], fixed_jacobian[j], 1.0e-10);
    }
  }
}