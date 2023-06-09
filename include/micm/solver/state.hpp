#pragma once

#include <cstddef>
#include <map>
#include <micm/util/matrix.hpp>
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

  struct Conditions
  {
    double temperature_{ 0.0 };
    double pressure_{ 0.0 };
    double air_density_{ 1.0 };
  };

  template<template<class> class M = Matrix>
  struct State
  {
    std::vector<Conditions> conditions_;
    std::map<std::string, std::size_t> variable_map_;
    M<double> variables_;
    M<double> custom_rate_parameters_;
    M<double> rate_constants_;

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
  };

  template<template<class> class M>
  inline State<M>::State()
      : conditions_(),
        variable_map_(),
        variables_(),
        custom_rate_parameters_(),
        rate_constants_()
  {
  }
  template<template<class> class M>
  inline State<M>::State(const std::size_t state_size, const std::size_t custom_parameters_size, const std::size_t process_size)
      : conditions_(1),
        variable_map_(),
        variables_(1, state_size, 0.0),
        custom_rate_parameters_(1, custom_parameters_size, 0.0),
        rate_constants_(1, process_size, 0.0)
  {
  }

  template<template<class> class M>
  inline State<M>::State(const StateParameters parameters)
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
}  // namespace micm