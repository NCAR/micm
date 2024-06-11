#pragma once

#include <micm/solver/rosenbrock.hpp>
#include <micm/util/sparse_matrix.hpp>

template<class MatrixPolicy, class SparseMatrixPolicy>
class HIRES
{
  std::size_t number_of_grid_cells_;
  std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements_;
  const std::vector<std::string> variable_names_ = { "y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8" };

 public:
  HIRES() = delete;

  HIRES(std::size_t number_of_grid_cells, std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements)
      : number_of_grid_cells_(number_of_grid_cells),
        nonzero_jacobian_elements_(nonzero_jacobian_elements)
  {
  }

  /// @brief Creates a new solver for the Oregonator system
  template<class SolverPolicy, class LinearSolverPolicy>
  static auto CreateSolver(auto parameters, std::size_t number_of_grid_cells)
  {
    parameters.relative_tolerance_ = 1e-3;
    parameters.absolute_tolerance_ = std::vector<double>(8, 1.0e-4 * parameters.relative_tolerance_);

    std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements;
    auto jacobian_builder = SparseMatrixPolicy::Create(8).SetNumberOfBlocks(number_of_grid_cells).InitialValue(0.0);
    for (int i = 0; i < 8; ++i)
    {
      for (int j = 0; j < 8; ++j)
      {
        jacobian_builder = jacobian_builder.WithElement(i, j);
        nonzero_jacobian_elements.insert(std::make_pair(i, j));
      }
    }
    SparseMatrixPolicy jacobian = SparseMatrixPolicy(jacobian_builder);
    LinearSolverPolicy linear_solver(jacobian, 1.0e-30);
    HIRES<MatrixPolicy, SparseMatrixPolicy> hires(number_of_grid_cells, nonzero_jacobian_elements);

    return SolverPolicy(parameters, std::move(linear_solver), std::move(hires), jacobian);
  }

  ~HIRES()
  {
  }

  micm::State<MatrixPolicy, SparseMatrixPolicy> GetState() const
  {
    auto state =
        micm::State<MatrixPolicy, SparseMatrixPolicy>{ { .number_of_grid_cells_ = number_of_grid_cells_,
                                                         .number_of_species_ = 8,
                                                         .number_of_rate_constants_ = 0,
                                                         .variable_names_ = variable_names_,
                                                         .nonzero_jacobian_elements_ = nonzero_jacobian_elements_ } };

    state.jacobian_ =
        micm::BuildJacobian<SparseMatrixPolicy>(nonzero_jacobian_elements_, number_of_grid_cells_, variable_names_.size());

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
  void AddForcingTerms(const MatrixPolicy& rate_constants, const MatrixPolicy& number_densities, MatrixPolicy& forcing)
  {
    auto data = number_densities.AsVector();

    forcing[0][0] += -1.71 * data[0] + 0.43 * data[1] + 8.32 * data[2] + 0.0007;
    forcing[0][1] += 1.71 * data[0] - 8.75 * data[1];
    forcing[0][2] += -10.03 * data[2] + 0.43 * data[3] + 0.035 * data[4];
    forcing[0][3] += 8.32 * data[1] + 1.71 * data[2] - 1.12 * data[3];
    forcing[0][4] += -1.745 * data[4] + 0.43 * data[5] + 0.43 * data[6];
    forcing[0][5] += -280.0 * data[5] * data[7] + 0.69 * data[3] + 1.71 * data[4] - 0.43 * data[5] + 0.69 * data[6];
    forcing[0][6] += 280.0 * data[5] * data[7] - 1.81 * data[6];
    forcing[0][7] += -forcing[0][6];
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

    jacobian[0][0][0] -= -1.71;
    jacobian[0][0][1] -= 0.43;
    jacobian[0][0][2] -= 8.32;

    jacobian[0][1][0] -= 1.71;
    jacobian[0][1][1] -= -8.75;

    jacobian[0][2][2] -= -10.03;
    jacobian[0][2][3] -= 0.43;
    jacobian[0][2][4] -= 0.035;

    jacobian[0][3][1] -= 8.32;
    jacobian[0][3][2] -= 1.71;
    jacobian[0][3][3] -= -1.12;

    jacobian[0][4][4] -= -1.745;
    jacobian[0][4][5] -= 0.43;
    jacobian[0][4][6] -= 0.43;

    jacobian[0][5][3] -= 0.69;
    jacobian[0][5][4] -= 1.71;
    jacobian[0][5][5] -= -0.43 - 280.0 * data[7];
    jacobian[0][5][6] -= 0.69;
    jacobian[0][5][7] -= -280.0 * data[6];

    jacobian[0][6][5] -= 280.0 * data[7];
    jacobian[0][6][6] -= -1.81;
    jacobian[0][6][7] -= 280.0 * data[6];

    jacobian[0][7][5] -= -280.0 * data[7];
    jacobian[0][7][6] -= 1.81;
    jacobian[0][7][7] -= -280.0 * data[6];
  }
};
