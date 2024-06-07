#pragma once

#include <micm/solver/rosenbrock.hpp>

template<class MatrixPolicy, class SparseMatrixPolicy>
class Oregonator
{
  std::size_t number_of_grid_cells_;
  std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements_;
  const std::vector<std::string> variable_names_ = { "A", "B", "C" };

 public:
  
  Oregonator() = delete;

  Oregonator(std::size_t number_of_grid_cells, std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements)
    : number_of_grid_cells_(number_of_grid_cells),
      nonzero_jacobian_elements_(nonzero_jacobian_elements)
  {
  }

  /// @brief Creates a new solver for the Oregonator system
  template<class SolverPolicy, class LinearSolverPolicy>
  static auto CreateSolver(auto parameters, std::size_t number_of_grid_cells)
  {
    
    parameters.relative_tolerance_ = 1e-4;
    parameters.absolute_tolerance_ = std::vector<double>(3, 1.0e-10);

    std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements;
    auto jacobian_builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(number_of_grid_cells).InitialValue(0.0);
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        jacobian_builder = jacobian_builder.WithElement(i, j);
        nonzero_jacobian_elements.insert(std::make_pair(i, j));
      }
    }
    SparseMatrixPolicy jacobian = SparseMatrixPolicy(jacobian_builder);

    return SolverPolicy(
      parameters,
      LinearSolverPolicy(jacobian, 1.0e-30),
      Oregonator<MatrixPolicy, SparseMatrixPolicy>(number_of_grid_cells, nonzero_jacobian_elements),
      jacobian);
  }

  ~Oregonator()
  {
  }

  micm::State<MatrixPolicy, SparseMatrixPolicy> GetState() const
  {
    auto state =
        micm::State<MatrixPolicy, SparseMatrixPolicy>{ { .number_of_grid_cells_ = number_of_grid_cells_,
                                                         .number_of_species_ = 3,
                                                         .number_of_rate_constants_ = 0,
                                                         .variable_names_ = variable_names_,
                                                         .nonzero_jacobian_elements_ = nonzero_jacobian_elements_ } };

    state.jacobian_ = micm::BuildJacobian<SparseMatrixPolicy>(
        nonzero_jacobian_elements_,
        number_of_grid_cells_,
        variable_names_.size());

    auto lu = micm::LuDecomposition::GetLUMatrices(state.jacobian_, 1.0e-30);
    auto lower_matrix = std::move(lu.first);
    auto upper_matrix = std::move(lu.second);
    state.lower_matrix_ = lower_matrix;
    state.upper_matrix_ = upper_matrix;

    return state;
  }

  /// @brief Calculate a chemical forcing
  /// @param rate_constants List of rate constants for each needed species
  /// @param number_densities The number density of each species
  /// @param forcing Vector of forcings for the current conditions
  void AddForcingTerms(
      const MatrixPolicy& rate_constants,
      const MatrixPolicy& number_densities,
      MatrixPolicy& forcing)
  {
    auto data = number_densities.AsVector();

    forcing[0][0] += 77.27 * (data[1] + data[0] * (1.0 - 8.375e-6 * data[0] - data[1]));
    forcing[0][1] += (data[2] - (1.0 + data[0]) * data[1]) / 77.27;
    forcing[0][2] += 0.161 * (data[0] - data[2]);
  }

  /// @brief Compute the derivative of the forcing w.r.t. each chemical, and return the negative jacobian
  /// @param rate_constants List of rate constants for each needed species
  /// @param number_densities The number density of each species
  /// @param jacobian The matrix of negative partial derivatives
  void SubtractJacobianTerms(
      const MatrixPolicy& rate_constants,
      const MatrixPolicy& number_densities,
      SparseMatrixPolicy& jacobian)
  {
    auto data = number_densities.AsVector();

    jacobian[0][0][0] -= 77.27 * (1. - 2. * 8.375e-6 * data[0] - data[1]);
    jacobian[0][0][1] -= 77.27 * (1. - data[0]);
    jacobian[0][0][2] -= 0;

    jacobian[0][1][0] -= -data[1] / 77.27;
    jacobian[0][1][1] -= -(1. + data[0]) / 77.27;
    jacobian[0][1][2] -= 1. / 77.27;

    jacobian[0][2][0] -= .161;
    jacobian[0][2][1] -= 0;
    jacobian[0][2][2] -= -.161;
  }
};
