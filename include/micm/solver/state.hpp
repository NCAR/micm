#pragma once

#include <cstddef>
#include <map>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// custom hash function for 'Species'
namespace std
{
  template<>
  struct hash<micm::Species>
  {
    size_t operator()(const micm::Species& key)
    {
      return hash<std::string>()(key.name_);
    }
  };
}  // namespace std

namespace micm
{

  struct StateParameters
  {
    std::vector<std::string> state_variable_names_{};
    std::size_t number_of_grid_cells_{ 1 };
    std::size_t number_of_custom_parameters_{ 0 };
    std::size_t number_of_rate_constants_{ 0 };
  };

  struct Conditions
  {
    double temperature_{ 0.0 };
    double pressure_{ 0.0 };
    double air_density_{ 1.0 };
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

    /// @brief Set concentrations
    /// @param species_to_concentration
    void SetConcentrations(
        const micm::System& system,
        const std::unordered_map<std::string, double>& species_to_concentration);

    // TODO: jiwon - 6/22 - can 'MUSICA name' be used as a key? Are they hashable (unique)?
    // or do we just want to use index? 
    // 
    void SetCustomRateParams();
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
  inline State<MatrixPolicy>::State(const std::size_t state_size, const std::size_t custom_parameters_size, const std::size_t process_size)
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
      const micm::System& system,
      const std::unordered_map<std::string, double>& species_to_concentration)
  {
    std::vector<double> concentrations;
    concentrations.reserve(system.gas_phase_.species_.size());

    for (auto& species : system.gas_phase_.species_)
    {
      auto species_ptr = species_to_concentration.find(species.name_);
      if (species_ptr == species_to_concentration.end())
      {
        throw std::invalid_argument("Invalid species: " + species.name_);
      }
      concentrations.push_back(species_ptr->second);
    }
 
    variables_[0] = concentrations;
  }
}  // namespace micm