#pragma once

#include <micm/solver/rosenbrock.hpp>
#include <micm/util/sparse_matrix.hpp>

template<class MatrixPolicy, class SparseMatrixPolicy>
class E5
{
  std::size_t number_of_grid_cells_;
  std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements_;
  const std::vector<std::string> variable_names_ = { "y1", "y2", "y3", "y4" };

 public:

  E5() = delete;

  E5(std::size_t number_of_grid_cells, std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements)
    : number_of_grid_cells_(number_of_grid_cells),
      nonzero_jacobian_elements_(nonzero_jacobian_elements)
  {
  }

  /// @brief Creates a new solver for the Oregonator system
  template<class SolverPolicy, class LinearSolverPolicy>
  static auto CreateSolver(auto parameters, std::size_t number_of_grid_cells)
  {
    parameters.relative_tolerance_ = 1e-2;
    parameters.absolute_tolerance_ = std::vector<double>(4, 1.7e-24);

    std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements;
    auto jacobian_builder = SparseMatrixPolicy::Create(4).SetNumberOfBlocks(number_of_grid_cells).InitialValue(0.0);
    for (int i = 0; i < 4; ++i)
    {
      for (int j = 0; j < 4; ++j)
      {
        jacobian_builder = jacobian_builder.WithElement(i, j);
        nonzero_jacobian_elements.insert(std::make_pair(i, j));
      }
    }
    SparseMatrixPolicy jacobian = SparseMatrixPolicy(jacobian_builder);

    return SolverPolicy(
      parameters,
      LinearSolverPolicy(jacobian, 1.0e-30),
      E5<MatrixPolicy, SparseMatrixPolicy>(number_of_grid_cells, nonzero_jacobian_elements),
      jacobian);
  }

  ~E5()
  {
  }

  micm::State<MatrixPolicy, SparseMatrixPolicy> GetState() const
  {
    auto state = micm::State<MatrixPolicy, SparseMatrixPolicy>{ {
      .number_of_grid_cells_ = number_of_grid_cells_,
      .number_of_species_ = 4,
      .number_of_rate_constants_ = 0,
      .variable_names_ = variable_names_,
      .nonzero_jacobian_elements_ = nonzero_jacobian_elements_
    } };

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

    double prod1 = 7.89e-10 * data[0];
    double prod2 = 1.1e7 * data[0] * data[2];
    double prod3 = 1.13e9 * data[1] * data[2];
    double prod4 = 1.13e3 * data[3];
    forcing[0][0] += -prod1 - prod2;
    forcing[0][1] += prod1 - prod3;
    forcing[0][3] += prod2 - prod4;
    forcing[0][2] += forcing[0][1] - forcing[0][3];
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

    double A = 7.89e-10;
    double B = 1.1e7;
    double CM = 1.13e9;
    double C = 1.13e3;
    jacobian[0][0][0] -= -A - B * data[2];
    jacobian[0][0][1] -= 0.0;
    jacobian[0][0][2] -= -B * data[0];
    jacobian[0][0][3] -= 0.0;
    jacobian[0][1][0] -= A;
    jacobian[0][1][1] -= -CM * data[2];
    jacobian[0][1][2] -= -CM * data[1];
    jacobian[0][1][3] -= 0.0;
    jacobian[0][2][0] -= A - B * data[2];
    jacobian[0][2][1] -= -CM * data[2];
    jacobian[0][2][2] -= -B * data[0] - CM * data[1];
    jacobian[0][2][3] -= C;
    jacobian[0][3][0] -= B * data[2];
    jacobian[0][3][1] -= 0.0;
    jacobian[0][3][2] -= B * data[0];
    jacobian[0][3][3] -= -C;
  }
};
