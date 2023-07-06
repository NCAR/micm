#pragma once

#include <cstddef>
#include <map>
#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/solver/conditions.hpp>
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
    std::size_t number_of_grid_cells_{ 1 };
    std::size_t number_of_custom_parameters_{ 0 };
    std::size_t number_of_rate_constants_{ 0 };
  };

  template<template<class> class MatrixPolicy = Matrix>
  struct State
  {
    std::vector<Conditions> conditions_;
    std::map<std::string, std::size_t> variable_map_;
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
    State(const StateParameters parameters);

    /// @brief Set species' concentrations
    /// @param species_to_concentration
    void SetConcentrations(
        const System& system,
        const std::unordered_map<std::string, double>& species_to_concentration);

    /// @brief Set photolysis rate constants
    /// @param photolysis rate
    void SetPhotolysisRate(
        const std::vector<PhotolysisRateConstant>& photolysis_rate_arr,
        const std::unordered_map<std::string, double>& photolysis_rate);
  };

  template<template<class> class MatrixPolicy>
  inline State<MatrixPolicy>::State()
      : conditions_(),
        variable_map_(),
        variables_(),
        custom_rate_parameters_(),
        rate_constants_()
  {
  }
  template<template<class> class MatrixPolicy>
  inline State<MatrixPolicy>::State(
      const std::size_t state_size,
      const std::size_t custom_parameters_size,
      const std::size_t process_size)
      : conditions_(1),
        variable_map_(),
        variables_(1, state_size, 0.0),
        custom_rate_parameters_(1, custom_parameters_size, 0.0),
        rate_constants_(1, process_size, 0.0)
  {
  }

  template<template<class> class MatrixPolicy>
  inline State<MatrixPolicy>::State(const StateParameters parameters)
      : conditions_(parameters.number_of_grid_cells_),
        variable_map_(),
        variables_(parameters.number_of_grid_cells_, parameters.state_variable_names_.size(), 0.0),
        custom_rate_parameters_(parameters.number_of_grid_cells_, parameters.number_of_custom_parameters_, 0.0),
        rate_constants_(parameters.number_of_grid_cells_, parameters.number_of_rate_constants_, 0.0)
  {
    std::size_t index = 0;
    for (auto& name : parameters.state_variable_names_)
      variable_map_[name] = index++;
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::SetConcentrations(
      const System& system,
      const std::unordered_map<std::string, double>& species_to_concentration)
  {
    std::vector<double> concentrations;
    concentrations.reserve(system.gas_phase_.species_.size());

    for (auto& species : system.gas_phase_.species_)
    {
      auto species_ptr = species_to_concentration.find(species.name_);
      if (species_ptr == species_to_concentration.end())
      {
        throw std::invalid_argument("Concentration value for '" + species.name_ + "' must be given.");
      }
      concentrations.push_back(species_ptr->second);
    }

    variables_[0] = concentrations;
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::SetPhotolysisRate(
      const std::vector<PhotolysisRateConstant>& photolysis_rate_arr,
      const std::unordered_map<std::string, double>& photolysis_rate)
  {
    std::vector<double> photo_rates;

    photo_rates.reserve(photolysis_rate_arr.size());

    for (auto& elem : photolysis_rate_arr)
    {
      auto rate_ptr = photolysis_rate.find(elem.name_);
      if (rate_ptr == photolysis_rate.end())
      {
        throw std::invalid_argument("Photolysis rate constant for '" + elem.name_ + "' must be given.");
      }

      photo_rates.push_back(rate_ptr->second);
    }

    custom_rate_parameters_[0] = photo_rates;
  }
}  // namespace micm