#pragma once

#include <micm/solver/rosenbrock.hpp>

template<template<class> class MatrixPolicy = micm::Matrix, template<class> class SparseMatrixPolicy = micm::SparseMatrix, class LinearSolverPolicy = micm::LinearSolver<double, SparseMatrixPolicy>>
class E5 : public micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>
{
 public:
  /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
  /// @param system The chemical system to create the solver for
  /// @param processes The collection of chemical processes that will be applied during solving
  E5(const micm::System& system,
     const std::vector<micm::Process>& processes,
     const micm::RosenbrockSolverParameters& parameters)
      : micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>()
  {
    this->system_ = system;
    this->processes_ = processes;
    this->parameters_ = parameters;
    this->N_ = this->system_.StateSize() * this->parameters_.number_of_grid_cells_;
    auto builder = SparseMatrixPolicy<double>::create(4).number_of_blocks(1).initial_value(0.0);
    for (int i = 0; i < 4; ++i)
    {
      for (int j = 0; j < 4; ++j)
      {
        builder = builder.with_element(i, j);
      }
    }
    this->jacobian_ = builder;
    for (std::size_t i = 0; i < this->jacobian_[0].size(); ++i)
      this->jacobian_diagonal_elements_.push_back(this->jacobian_.VectorIndex(0, i, i));
    this->linear_solver_ = LinearSolverPolicy(this->jacobian_, 1.0e-30);
  }

  ~E5()
  {
  }

  /// @brief Calculate a chemical forcing
  /// @param rate_constants List of rate constants for each needed species
  /// @param number_densities The number density of each species
  /// @param forcing Vector of forcings for the current conditions
  void CalculateForcing(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      MatrixPolicy<double>& forcing) override
  {
    std::fill(forcing.AsVector().begin(), forcing.AsVector().end(), 0.0);
    this->stats_.function_calls += 1;

    auto data = number_densities.AsVector();

    double prod1 = 7.89e-10 * data[0];
    double prod2 = 1.1e7 * data[0] * data[2];
    double prod3 = 1.13e9 * data[1] * data[2];
    double prod4 = 1.13e3 * data[3];
    forcing[0][0] = -prod1 - prod2;
    forcing[0][1] = prod1 - prod3;
    forcing[0][3] = prod2 - prod4;
    forcing[0][2] = forcing[0][1] - forcing[0][3];
  }

  /// @brief Compute the derivative of the forcing w.r.t. each chemical, the jacobian
  /// @param rate_constants List of rate constants for each needed species
  /// @param number_densities The number density of each species
  /// @param jacobian The matrix of partial derivatives
  void CalculateJacobian(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      SparseMatrixPolicy<double>& jacobian) override
  {
    std::fill(jacobian.AsVector().begin(), jacobian.AsVector().end(), 0.0);
    auto data = number_densities.AsVector();
    this->stats_.jacobian_updates += 1;

    double A = 7.89e-10;
    double B = 1.1e7;
    double CM = 1.13e9;
    double C = 1.13e3;
    jacobian[0][0][0] = -A - B * data[2];
    jacobian[0][0][1] = 0.0;
    jacobian[0][0][2] = -B * data[0];
    jacobian[0][0][3] = 0.0;
    jacobian[0][1][0] = A;
    jacobian[0][1][1] = -CM * data[2];
    jacobian[0][1][2] = -CM * data[1];
    jacobian[0][1][3] = 0.0;
    jacobian[0][2][0] = A - B * data[2];
    jacobian[0][2][1] = -CM * data[2];
    jacobian[0][2][2] = -B * data[0] - CM * data[1];
    jacobian[0][2][3] = C;
    jacobian[0][3][0] = B * data[2];
    jacobian[0][3][1] = 0.0;
    jacobian[0][3][2] = B * data[0];
    jacobian[0][3][3] = -C;
  }
};
