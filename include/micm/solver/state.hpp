// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/solver/lu_decomposition.hpp>
#include <micm/solver/temporary_variables.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <map>
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
    std::size_t number_of_species_{ 0 };
    std::size_t number_of_rate_constants_{ 0 };
    std::vector<std::string> variable_names_{};
    std::vector<std::string> custom_rate_parameter_labels_{};
    std::set<std::pair<std::size_t, std::size_t>> nonzero_jacobian_elements_{};
  };

  template<class DenseMatrixPolicy = StandardDenseMatrix, class SparseMatrixPolicy = StandardSparseMatrix>
  struct State
  {
    /// Type of the DenseMatrixPolicy
    using DenseMatrixPolicyType = DenseMatrixPolicy;

    /// @brief The concentration of chemicals, varies through time
    DenseMatrixPolicy variables_;
    /// @brief Rate paramters particular to user-defined rate constants, may vary in time
    DenseMatrixPolicy custom_rate_parameters_;
    /// @brief The reaction rates, may vary in time
    DenseMatrixPolicy rate_constants_;
    /// @brief Atmospheric conditions, varies in time
    std::vector<Conditions> conditions_;
    /// @brief The jacobian structure, varies for each solve
    SparseMatrixPolicy jacobian_;
    /// @brief Immutable data required for the state
    std::map<std::string, std::size_t> variable_map_;
    std::map<std::string, std::size_t> custom_rate_parameter_map_;
    std::vector<std::string> variable_names_{};
    SparseMatrixPolicy lower_matrix_;
    SparseMatrixPolicy upper_matrix_;
    std::size_t state_size_;
    std::size_t number_of_grid_cells_;
    std::unique_ptr<TemporaryVariables> temporary_variables_;

    /// @brief Constructor with parameters
    /// @param parameters State dimension information
    State(const StateParameters& parameters);

    /// @brief Copy constructor
    /// @param other The state object to be copied
    State(const State& other)
    {
      variables_ = other.variables_;
      custom_rate_parameters_ = other.custom_rate_parameters_;
      rate_constants_ = other.rate_constants_;
      conditions_ = other.conditions_;
      jacobian_ = other.jacobian_;
      variable_map_ = other.variable_map_;
      custom_rate_parameter_map_ = other.custom_rate_parameter_map_;
      variable_names_ = other.variable_names_;
      lower_matrix_ = other.lower_matrix_;
      upper_matrix_ = other.upper_matrix_;
      state_size_ = other.state_size_;
      number_of_grid_cells_ = other.number_of_grid_cells_;
      temporary_variables_ = std::make_unique<TemporaryVariables>(*other.temporary_variables_);
    }

    /// @brief Assignment operator
    /// @param other The state object to be assigned
    /// @return Reference to the assigned state object
    State& operator=(const State& other)
    {
      if (this != &other)
      {
        variables_ = other.variables_;
        custom_rate_parameters_ = other.custom_rate_parameters_;
        rate_constants_ = other.rate_constants_;
        conditions_ = other.conditions_;
        jacobian_ = other.jacobian_;
        variable_map_ = other.variable_map_;
        custom_rate_parameter_map_ = other.custom_rate_parameter_map_;
        variable_names_ = other.variable_names_;
        lower_matrix_ = other.lower_matrix_;
        upper_matrix_ = other.upper_matrix_;
        state_size_ = other.state_size_;
        number_of_grid_cells_ = other.number_of_grid_cells_;
        temporary_variables_ = std::make_unique<TemporaryVariables>(*other.temporary_variables_);
      }
      return *this;
    }

    /// @brief Move constructor
    /// @param other The state object to be moved
    State(State&& other) noexcept
        : variables_(std::move(other.variables_)),
          custom_rate_parameters_(std::move(other.custom_rate_parameters_)),
          rate_constants_(std::move(other.rate_constants_)),
          conditions_(std::move(other.conditions_)),
          jacobian_(std::move(other.jacobian_)),
          variable_map_(std::move(other.variable_map_)),
          custom_rate_parameter_map_(std::move(other.custom_rate_parameter_map_)),
          variable_names_(std::move(other.variable_names_)),
          lower_matrix_(std::move(other.lower_matrix_)),
          upper_matrix_(std::move(other.upper_matrix_)),
          state_size_(other.state_size_),
          number_of_grid_cells_(other.number_of_grid_cells_),
          temporary_variables_(std::move(other.temporary_variables_))
    {
    }

    /// @brief Move assignment operator
    /// @param other The state object to be moved
    /// @return Reference to the moved state object
    State& operator=(State&& other) noexcept
    {
      if (this != &other)
      {
        variables_ = std::move(other.variables_);
        custom_rate_parameters_ = std::move(other.custom_rate_parameters_);
        rate_constants_ = std::move(other.rate_constants_);
        conditions_ = std::move(other.conditions_);
        jacobian_ = std::move(other.jacobian_);
        variable_map_ = std::move(other.variable_map_);
        custom_rate_parameter_map_ = std::move(other.custom_rate_parameter_map_);
        variable_names_ = std::move(other.variable_names_);
        lower_matrix_ = std::move(other.lower_matrix_);
        upper_matrix_ = std::move(other.upper_matrix_);
        state_size_ = other.state_size_;
        number_of_grid_cells_ = other.number_of_grid_cells_;
        temporary_variables_ = std::move(other.temporary_variables_);

        other.state_size_ = 0;
        other.number_of_grid_cells_ = 0;
      }
      return *this;
    }

    /// @brief Set species' concentrations
    /// @param species_to_concentration
    void SetConcentrations(const std::unordered_map<std::string, std::vector<double>>& species_to_concentration);

    /// @brief Set a single species concentration
    /// @param species the species to set the concentration for
    /// @param concentration concentration(s) [mol m-3]
    void SetConcentration(const Species& species, double concentration);
    void SetConcentration(const Species& species, const std::vector<double>& concentration);

    /// @brief Set custom parameters assuming the values are properly ordered
    /// @param parameters map of custom rate parameters
    void UnsafelySetCustomRateParameters(const std::vector<std::vector<double>>& parameters);

    /// @brief Set custom parameters for rate constant calculations by label
    /// @param parameters map of custom rate parameters
    void SetCustomRateParameters(const std::unordered_map<std::string, std::vector<double>>& parameters);

    /// @brief Set a single custom rate constant parameter
    /// @param label parameter label
    /// @param value new parameter value
    void SetCustomRateParameter(const std::string& label, double value);
    void SetCustomRateParameter(const std::string& label, const std::vector<double>& values);

    /// @brief Print a header of species to display concentrations with respect to time
    void PrintHeader();

    /// @brief Print state (concentrations) at the given time
    /// @param time solving time
    void PrintState(double time);

    /// @brief Default constructor
    /// Only defined to be used to create default values in types, but a default constructed state is not useable
    State() : {};
  };

}  // namespace micm

#include "state.inl"
