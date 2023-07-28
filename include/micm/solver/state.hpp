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

  template<template<class> class MatrixPolicy>
  inline State<MatrixPolicy>::State()
      : conditions_(),
        variable_map_(),
        custom_rate_parameter_map_(),
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
        custom_rate_parameter_map_(),
        variables_(1, state_size, 0.0),
        custom_rate_parameters_(1, custom_parameters_size, 0.0),
        rate_constants_(1, process_size, 0.0)
  {
  }

  template<template<class> class MatrixPolicy>
  inline State<MatrixPolicy>::State(const StateParameters parameters)
      : conditions_(parameters.number_of_grid_cells_),
        variable_map_(),
        custom_rate_parameter_map_(),
        variables_(parameters.number_of_grid_cells_, parameters.state_variable_names_.size(), 0.0),
        custom_rate_parameters_(parameters.number_of_grid_cells_, parameters.custom_rate_parameter_labels_.size(), 0.0),
        rate_constants_(parameters.number_of_grid_cells_, parameters.number_of_rate_constants_, 0.0)
  {
    std::size_t index = 0;
    for (auto& name : parameters.state_variable_names_)
      variable_map_[name] = index++;
    index = 0;
    for (auto& label : parameters.custom_rate_parameter_labels_)
      custom_rate_parameter_map_[label] = index++;
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
        throw std::runtime_error("Concentration value(s) for '" + species.name_ + "' must be given.");
      }
      num_concentrations_per_species.push_back(species_ptr->second.size());
    }

    // Check if number of concentraiton inputs are the same for all species
    if (!std::all_of(
            num_concentrations_per_species.begin(),
            num_concentrations_per_species.end(),
            [&](int& i) { return i == num_concentrations_per_species.front(); }))
    {
      throw std::runtime_error("Concentration value must be given to all sets of grid cells.");
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
  void State<MatrixPolicy>::SetConcentration(const Species& species, double concentration)
  {
    if (variables_.size() != 1)
      throw std::invalid_argument("Incorrect number of concentration values passed to multi-gridcell State");
    variables_[0][variable_map_[species.name_]] = concentration;
  }

  template<template<class> class MatrixPolicy>
  void State<MatrixPolicy>::SetConcentration(const Species& species, const std::vector<double>& concentration)
  {
    if (variables_.size() != concentration.size())
      throw std::invalid_argument("Incorrect number of concentration values passed to multi-gridcell State");
    std::size_t i_species = variable_map_[species.name_];
    for (std::size_t i = 0; i < variables_.size(); ++i)
      variables_[i][i_species] = concentration[i];
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::SetCustomRateParameters(
      const std::unordered_map<std::string, std::vector<double>>& parameters)
  {
    for (auto& pair : parameters)
      SetCustomRateParameter(pair.first, pair.second);
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::SetCustomRateParameter(const std::string& label, double value)
  {
    auto param = custom_rate_parameter_map_.find(label);
    if (param == custom_rate_parameter_map_.end())
      throw std::invalid_argument("Unkonwn rate constant parameter '" + label + "'");
    if (custom_rate_parameters_.size() != 1)
      throw std::invalid_argument("Incorrect number of custom rate parameter values passed to multi-gridcell State");
    custom_rate_parameters_[0][param->second] = value;
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::SetCustomRateParameter(const std::string& label, const std::vector<double>& values)
  {
    auto param = custom_rate_parameter_map_.find(label);
    if (param == custom_rate_parameter_map_.end())
      throw std::invalid_argument("Unkonwn rate constant parameter '" + label + "'");
    if (custom_rate_parameters_.size() != values.size())
      throw std::invalid_argument("Incorrect number of custom rate parameter values passed to multi-gridcell State");
    for (std::size_t i = 0; i < custom_rate_parameters_.size(); ++i)
      custom_rate_parameters_[i][param->second] = values[i];
  }
}  // namespace micm