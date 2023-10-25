#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <micm/system/conditions.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace micm
{

  /// @brief Invariants that can be used to construct a state
  struct StateParameters
  {
    std::size_t number_of_grid_cells_{ 1 };
    std::size_t number_of_rate_constants_{ 0 };
    std::map<std::string, std::size_t> custom_rate_parameter_map_;
    std::vector<std::string> variable_names_{};
    std::vector<std::string> custom_rate_parameter_labels_{};
    std::vector<std::size_t> jacobian_diagonal_elements_;
  };

  template<
    template<class> class MatrixPolicy = Matrix, 
    template<class> class SparseMatrixPolicy = StandardSparseMatrix>
  struct State
  {
    /// @brief The concentration of chemicals, varies through time
    MatrixPolicy<double> variables_;
    /// @brief Rate paramters particular to user-defined rate constants, may vary in time 
    MatrixPolicy<double> custom_rate_parameters_;
    /// @brief The reaction rates, may vary in time
    MatrixPolicy<double> rate_constants_;
    /// @brief Atmospheric conditions, varies in time
    std::vector<Conditions> conditions_;
    /// @brief The jacobian structure, varies for each solve
    SparseMatrixPolicy<double> jacobian_;
    /// @brief Immutable data required for the state
    std::map<std::string, std::size_t> variable_map_;
    std::map<std::string, std::size_t> custom_rate_parameter_map_;
    std::vector<std::string> variable_names_{};
    SparseMatrixPolicy<double> lower_matrix_;
    SparseMatrixPolicy<double> upper_matrix_;

    /// @brief
    State();

    /// @brief
    /// @param state_size The number of System state variables
    /// @param custom_parameters_size The number of custom rate parameters
    /// @param process_size The number of processes to store rate constants for
    State(const std::size_t state_size, const std::size_t custom_parameters_size, const std::size_t process_size);

    /// @brief
    /// @param parameters State dimension information
    State(const StateParameters& parameters);

    /// @brief Set species' concentrations
    /// @param species_to_concentration
    void SetConcentrations(const std::unordered_map<std::string, std::vector<double>>& species_to_concentration);

    /// @brief Set a single species concentration
    /// @param species the species to set the concentration for
    /// @param concentration concentration(s) [mol m-3]
    void SetConcentration(const Species& species, double concentration);
    void SetConcentration(const Species& species, const std::vector<double>& concentration);

    /// @brief Set custom parameters for rate constant calculations by label
    /// @param parameters map of custom rate parameters
    void SetCustomRateParameters(const std::unordered_map<std::string, std::vector<double>>& parameters);

    /// @brief Set a single custom rate constant parameter
    /// @param label parameter label
    /// @param value new parameter value
    void SetCustomRateParameter(const std::string& label, double value);
    void SetCustomRateParameter(const std::string& label, const std::vector<double>& values);
  };

}  // namespace micm

#include "state.inl"