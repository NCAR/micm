#pragma once

#include <micm/solver/rosenbrock.hpp>

template<template<class> class MatrixPolicy = micm::Matrix, template<class> class SparseMatrixPolicy = micm::SparseMatrix, class LinearSolverPolicy = micm::LinearSolver<double, SparseMatrixPolicy>>
class Oregonator : public micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>
{
 public:
  /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
  /// @param system The chemical system to create the solver for
  /// @param processes The collection of chemical processes that will be applied during solving
  Oregonator(
      const micm::System& system,
      const std::vector<micm::Process>& processes,
      const micm::RosenbrockSolverParameters& parameters)
      : micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>()
  {
    this->system_ = system;
    this->processes_ = processes;
    this->parameters_ = parameters;
    this->N_ = this->system_.StateSize() * this->parameters_.number_of_grid_cells_;
    auto builder = SparseMatrixPolicy<double>::create(3).number_of_blocks(1).initial_value(0.0);
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        builder = builder.with_element(i, j);
      }
    }
    this->jacobian_ = builder;
    for (std::size_t i = 0; i < this->jacobian_[0].size(); ++i)
      this->jacobian_diagonal_elements_.push_back(this->jacobian_.VectorIndex(0, i, i));
    this->linear_solver_ = LinearSolverPolicy(this->jacobian_, 1.0e-30);
  }

  ~Oregonator()
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

    auto data = number_densities.AsVector();

    forcing[0][0] = 77.27 * (data[1] + data[0] * (1.0 - 8.375e-6 * data[0] - data[1]));
    forcing[0][1] = (data[2] - (1.0 + data[0]) * data[1]) / 77.27;
    forcing[0][2] = 0.161 * (data[0] - data[2]);
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
    auto data = number_densities.AsVector();

    jacobian[0][0][0] = 77.27 * (1. - 2. * 8.375e-6 * data[0] - data[1]);
    jacobian[0][0][1] = 77.27 * (1. - data[0]);
    jacobian[0][0][2] = 0;

    jacobian[0][1][0] = -data[1] / 77.27;
    jacobian[0][1][1] = -(1. + data[0]) / 77.27;
    jacobian[0][1][2] = 1. / 77.27;

    jacobian[0][2][0] = .161;
    jacobian[0][2][1] = 0;
    jacobian[0][2][2] = -.161;
  }
};
