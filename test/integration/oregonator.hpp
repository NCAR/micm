#pragma once

#include <micm/solver/rosenbrock.hpp>

template<template<class> class MatrixPolicy = micm::Matrix, template<class> class SparseMatrixPolicy = micm::SparseMatrix>
class Oregonator : public micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>
{
 public:
  /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
  /// @param system The chemical system to create the solver for
  /// @param processes The collection of chemical processes that will be applied during solving
  Oregonator(
      const micm::System& system,
      const std::vector<micm::Process>& processes,
      const micm::RosenbrockSolverParameters& parameters)
      : micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>()
  {
    this->system_ = system;
    this->processes_ = processes;
    this->parameters_ = parameters;
    this->N_ = this->system_.StateSize() * this->parameters_.number_of_grid_cells_;
    this->jacobian_ = SparseMatrixPolicy<double>::create(9).number_of_blocks(1);
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
    auto force = forcing.AsVector();

    force[0] = 77.27 * (data[1] + data[0] * (1.0 - 8.375e-6 * data[0] - data[1]));
    force[1] = (data[2] - (1.0 + data[0]) * data[1]) / 77.27;
    force[2] = 0.161 * (data[0] - data[2]);
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
    auto jac = jacobian.AsVector();

    jac[0] = 77.27 * (1. - 2. * 8.375e-6 * data[0] - data[1]);
    jac[1] = 77.27 * (1. - data[0]);
    jac[2] = 0;

    jac[3] = -data[1] / 77.27;
    jac[4] = -(1. + data[0]) / 77.27;
    jac[5] = 1. / 77.27;

    jac[6] = .161;
    jac[6] = 0;
    jac[6] = -.161;
  }
};