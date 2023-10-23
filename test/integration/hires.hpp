#pragma once

#include <micm/solver/rosenbrock.hpp>

template<
    template<class> class MatrixPolicy = micm::Matrix,
    template<class> class SparseMatrixPolicy = micm::SparseMatrix,
    class LinearSolverPolicy = micm::LinearSolver<double, SparseMatrixPolicy>>
class HIRES : public micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>
{
 public:
  /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
  /// @param system The chemical system to create the solver for
  /// @param processes The collection of chemical processes that will be applied during solving
  HIRES(
      const micm::System& system,
      const std::vector<micm::Process>& processes,
      const micm::RosenbrockSolverParameters& parameters)
      : micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>()
  {
    this->system_ = system;
    this->processes_ = processes;
    this->parameters_ = parameters;
    this->N_ = this->system_.StateSize() * this->parameters_.number_of_grid_cells_;
    auto builder = SparseMatrixPolicy<double>::create(8).number_of_blocks(1).initial_value(0.0);
    for (int i = 0; i < 8; ++i)
    {
      for (int j = 0; j < 8; ++j)
      {
        builder = builder.with_element(i, j);
      }
    }
    this->jacobian_ = builder;
    for (std::size_t i = 0; i < this->jacobian_[0].size(); ++i)
      this->jacobian_diagonal_elements_.push_back(this->jacobian_.VectorIndex(0, i, i));
    this->linear_solver_ = LinearSolverPolicy(this->jacobian_, 1.0e-30);
  }

  ~HIRES()
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

    forcing[0][0] = -1.71 * data[0] + 0.43 * data[1] + 8.32 * data[2] + 0.0007;
    forcing[0][1] = 1.71 * data[0] - 8.75 * data[1];
    forcing[0][2] = -10.03 * data[2] + 0.43 * data[3] + 0.035 * data[4];
    forcing[0][3] = 8.32 * data[1] + 1.71 * data[2] - 1.12 * data[3];
    forcing[0][4] = -1.745 * data[4] + 0.43 * data[5] + 0.43 * data[6];
    forcing[0][5] = -280.0 * data[5] * data[7] + 0.69 * data[3] + 1.71 * data[4] - 0.43 * data[5] + 0.69 * data[6];
    forcing[0][6] = 280.0 * data[5] * data[7] - 1.81 * data[6];
    forcing[0][7] = -forcing[0][6];
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

    jacobian[0][0][0] = -1.71;
    jacobian[0][0][1] = 0.43;
    jacobian[0][0][2] = 8.32;

    jacobian[0][1][0] = 1.71;
    jacobian[0][1][1] = -8.75;

    jacobian[0][2][2] = -10.03;
    jacobian[0][2][3] = 0.43;
    jacobian[0][2][4] = 0.035;

    jacobian[0][3][1] = 8.32;
    jacobian[0][3][2] = 1.71;
    jacobian[0][3][3] = -1.12;

    jacobian[0][4][4] = -1.745;
    jacobian[0][4][5] = 0.43;
    jacobian[0][4][6] = 0.43;

    jacobian[0][5][3] = 0.69;
    jacobian[0][5][4] = 1.71;
    jacobian[0][5][5] = -0.43 - 280.0 * data[7];
    jacobian[0][5][6] = 0.69;
    jacobian[0][5][7] = -280.0 * data[6];

    jacobian[0][6][5] = 280.0 * data[7];
    jacobian[0][6][6] = -1.81;
    jacobian[0][6][7] = 280.0 * data[6];

    jacobian[0][7][5] = -280.0 * data[7];
    jacobian[0][7][6] = 1.81;
    jacobian[0][7][7] = -280.0 * data[6];
  }
};
