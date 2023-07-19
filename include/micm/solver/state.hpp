#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <micm/process/photolysis_rate_constant.hpp>
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
        const std::unordered_map<std::string, std::vector<double>>& species_to_concentration);

    /// @brief Set photolysis rate constants
    /// @param photolysis rate
    void SetPhotolysisRate(
        const std::vector<PhotolysisRateConstant>& photolysis_rate_arr,
        const std::unordered_map<std::string, std::vector<double>>& photolysis_rate);
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
      const std::unordered_map<std::string, std::vector<double>>& species_to_concentration)
  {
    int num_set_grid_cells = 0;
    unsigned num_species = system.gas_phase_.species_.size();

    std::vector<int> num_concentrations_per_species;
    num_concentrations_per_species.reserve(num_species);

    // Iterate map to store the number of concentration values corresponding to the number of set of grid cells
    for (auto& species : system.gas_phase_.species_)
    {
      auto species_ptr = species_to_concentration.find(species.name_);
      if (species_ptr == species_to_concentration.end())
      {
        throw std::invalid_argument("Concentration value(s) for '" + species.name_ + "' must be given.");
      }
      num_concentrations_per_species.push_back(species_ptr->second.size());
    }

    // Check if number of concentraiton inputs are the same for all species
    if (!std::all_of(
            num_concentrations_per_species.begin(),
            num_concentrations_per_species.end(),
            [&](int& i) { return i == num_concentrations_per_species.front(); }))
    {
      throw std::invalid_argument(
          "Concentration value must be given to all sets of grid cells.");  // TODO: jiwon 7/10 - error message
    }

    num_set_grid_cells = num_concentrations_per_species[0];

    // Find species and iterate through the keys to store concentrations for each set of grid cells
    // 'concentrations' represents an N-D array in contiguous memory (N = num_set_grid_cells)
    std::vector<double> concentrations;
    concentrations.resize(num_species * num_set_grid_cells);

    for (int i = 0; i < num_species; i++)
    {
      auto species_ptr = species_to_concentration.find(system.gas_phase_.species_[i].name_);

      for (int j = 0; j < num_set_grid_cells; j++)
      {
        concentrations[i + num_species * j] = species_ptr->second[j];
      }
    }

    // Extract sub vector to assign to the corresponding set of grid cells.
    // TODO: jiwon 7/12 - I think we want to reduce copy operations here
    std::vector<double> sub_concentrations;
    sub_concentrations.reserve(num_species);

    for (int i = 0; i < num_set_grid_cells; i++)
    {
      sub_concentrations = { concentrations.begin() + (num_species * i),
                             concentrations.begin() + (num_species * i) + num_species };
      variables_[i] = sub_concentrations;
    }
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::SetPhotolysisRate(
      const std::vector<PhotolysisRateConstant>& photolysis_rate_arr,
      const std::unordered_map<std::string, std::vector<double>>& photolysis_rate)
  {
    int num_set_grid_cells = 0;
    unsigned num_photo_values = photolysis_rate_arr.size();

    std::vector<int> num_values_per_key;
    num_values_per_key.reserve(num_photo_values);

    // Iterate map to store the number of rate constant values corresponding to the number of set of grid cells
    for (auto& elem : photolysis_rate_arr)
    {
      auto rate_ptr = photolysis_rate.find(elem.name_);
      if (rate_ptr == photolysis_rate.end())
      {
        throw std::invalid_argument("Photolysis rate constant(s) for '" + elem.name_ + "' must be given.");
      }
      num_values_per_key.push_back(rate_ptr->second.size());
    }

    // Check if number of rate constants inputs are the same for all photolysis rate constant objects.
    if (!std::all_of(
            num_values_per_key.begin(), num_values_per_key.end(), [&](int& i) { return i == num_values_per_key.front(); }))
    {
      throw std::invalid_argument(
          "Photolysis rate constant value must be given to all sets of grid cells.");  // TODO: jiwon 7/10 - error message
    }

    num_set_grid_cells = num_values_per_key[0];

    // Find rate constants and iterate through the keys to store the rate constants for each set of grid cells
    // 'photo_rates' represents an N-D array in contiguous memory (N = num_set_grid_cells)
    std::vector<double> photo_rates;
    photo_rates.resize(num_photo_values * num_set_grid_cells);

    for (int i = 0; i < num_photo_values; i++)
    {
      auto rate_ptr = photolysis_rate.find(photolysis_rate_arr[i].name_);
      for (int j = 0; j < num_set_grid_cells; j++)
      {
        photo_rates[i + num_photo_values * j] = rate_ptr->second[j];
      }
    }

    // Extract sub vector to assign to the corresponding set of grid cells.
    // TODO: jiwon 7/12 - I think we want to reduce copy operations here
    std::vector<double> sub_photo_rates;
    sub_photo_rates.reserve(num_photo_values);

    for (int i = 0; i < num_set_grid_cells; i++)
    {
      sub_photo_rates = { photo_rates.begin() + (num_photo_values * i),
                          photo_rates.begin() + (num_photo_values * i) + num_photo_values };
      custom_rate_parameters_[i] = sub_photo_rates;
    }
  }
}  // namespace micm