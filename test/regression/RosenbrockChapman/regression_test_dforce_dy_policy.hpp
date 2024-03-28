#pragma once

#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include <random>

#include "chapman_ode_solver.hpp"
#include "util.hpp"

template<class OdeSolverPolicy>
void testJacobian(OdeSolverPolicy& solver)
{
  std::random_device rnd_device;
  std::mt19937 engine{ rnd_device() };
  std::lognormal_distribution dist{ -2.0, 4.0 };

  micm::ChapmanODESolver fixed_solver{};

  auto state = solver.GetState();

  auto& state_vec = state.variables_.AsVector();
  std::generate(state_vec.begin(), state_vec.end(), [&]() { return dist(engine); });
  auto& rate_const_vec = state.rate_constants_.AsVector();
  std::generate(state_vec.begin(), state_vec.end(), [&]() { return dist(engine); });

  auto& jacobian = state.jacobian_;
  solver.CalculateNegativeJacobian(state.rate_constants_, state.variables_, jacobian);

  for (std::size_t i{}; i < 3; ++i)
  {
    double number_density_air = 1.0;
    std::vector<double> rate_constants = state.rate_constants_[i];
    std::vector<double> variables(state.variables_[i].size());
    for (std::size_t j{}; j < state.variables_[i].size(); ++j)
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

template<class T>
using DenseMatrix = micm::Matrix<T>;
template<class T>
using SparseMatrix = micm::SparseMatrix<T>;

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;