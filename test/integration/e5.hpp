#pragma once

#include <micm/solver/rosenbrock.hpp>

template<
    template<class> class MatrixPolicy = micm::Matrix,
    template<class> class SparseMatrixPolicy = micm::SparseMatrix,
    class LinearSolverPolicy = micm::LinearSolver<double, SparseMatrixPolicy>>
class E5 : public micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>
{
  std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements_;

 public:
  /// @brief Builds a Rosenbrock solver for the given system, processes, and solver parameters
  /// @param system The chemical system to create the solver for
  /// @param processes The collection of chemical processes that will be applied during solving
  E5(const micm::System& system,
     const std::vector<micm::Process>& processes,
     const micm::RosenbrockSolverParameters& parameters)
      : micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>()
  {
    this->processes_ = processes;
    this->parameters_ = parameters;

    auto builder = SparseMatrixPolicy<double>::create(4).number_of_blocks(1).initial_value(0.0);
    for (int i = 0; i < 4; ++i)
    {
      for (int j = 0; j < 4; ++j)
      {
        builder = builder.with_element(i, j);
        nonzero_jacobian_elements_.insert(std::make_pair(i, j));
      }
    }
    SparseMatrixPolicy<double> jacobian = SparseMatrixPolicy<double>(builder);

    std::vector<std::size_t> jacobian_diagonal_elements;
    for (std::size_t i = 0; i < jacobian[0].size(); ++i)
      jacobian_diagonal_elements.push_back(jacobian.VectorIndex(0, i, i));

    std::vector<std::string> param_labels{};
    for (const auto& process : processes)
      if (process.rate_constant_)
        for (auto& label : process.rate_constant_->CustomParameters())
          param_labels.push_back(label);

    std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> state_reordering;
    this->state_parameters_ = {
      .number_of_grid_cells_ = 1,
      .number_of_rate_constants_ = processes.size(),
      .variable_names_ = system.UniqueNames(state_reordering),
      .custom_rate_parameter_labels_ = param_labels,
      .jacobian_diagonal_elements_ = jacobian_diagonal_elements,
    };

    this->linear_solver_ = LinearSolverPolicy(jacobian, 1.0e-30);
  }

  ~E5()
  {
  }

  micm::State<MatrixPolicy, SparseMatrixPolicy> GetState() const override
  {
    auto state = micm::State<MatrixPolicy, SparseMatrixPolicy>{ this->state_parameters_ };

    state.jacobian_ = micm::build_jacobian<SparseMatrixPolicy>(
        nonzero_jacobian_elements_,
        this->state_parameters_.number_of_grid_cells_,
        this->state_parameters_.variable_names_.size());

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
  void CalculateForcing(
      const MatrixPolicy<double>& rate_constants,
      const MatrixPolicy<double>& number_densities,
      MatrixPolicy<double>& forcing) override
  {
    std::fill(forcing.AsVector().begin(), forcing.AsVector().end(), 0.0);

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
