#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <micm/system/conditions.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace micm
{

  struct StateParameters
  {
    std::vector<std::string> state_variable_names_{};
    std::vector<std::string> custom_rate_parameter_labels_{};
    std::size_t number_of_grid_cells_{ 1 };
    std::size_t number_of_rate_constants_{ 0 };
  };

  template<template<class> class MatrixPolicy = Matrix>
  struct State
  {
    std::vector<Conditions> conditions_;
    std::map<std::string, std::size_t> variable_map_;
    std::map<std::string, std::size_t> custom_rate_parameter_map_;
    std::vector<std::string> variable_names_{};
    MatrixPolicy<double> variables_;
    MatrixPolicy<double> custom_rate_parameters_;
    MatrixPolicy<double> rate_constants_;

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
    void SetConcentrations(
        const std::unordered_map<std::string, std::vector<double>>& species_to_concentration);

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