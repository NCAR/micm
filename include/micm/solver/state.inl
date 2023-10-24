// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  template<template<class> class MatrixPolicy>
  inline State<MatrixPolicy>::State()
      : conditions_(),
        variable_map_(),
        custom_rate_parameter_map_(),
        variable_names_(),
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
        variable_names_(),
        variables_(1, state_size, 0.0),
        custom_rate_parameters_(1, custom_parameters_size, 0.0),
        rate_constants_(1, process_size, 0.0)
  {
  }

  template<template<class> class MatrixPolicy>
  inline State<MatrixPolicy>::State(const StateParameters& parameters)
      : conditions_(parameters.number_of_grid_cells_),
        variable_map_(),
        custom_rate_parameter_map_(),
        variable_names_(parameters.state_variable_names_),
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
      const std::unordered_map<std::string, std::vector<double>>& species_to_concentration)
  {
    const int num_grid_cells = conditions_.size();
    for (const auto& pair : species_to_concentration)
      SetConcentration({ pair.first }, pair.second);
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::SetConcentration(const Species& species, double concentration)
  {
    auto var = variable_map_.find(species.name_);
    if (var == variable_map_.end())
      throw std::invalid_argument("Unknown variable '" + species.name_ + "'");
    if (variables_.size() != 1)
      throw std::invalid_argument("Incorrect number of concentration values passed to multi-gridcell State");
    variables_[0][variable_map_[species.name_]] = concentration;
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::SetConcentration(const Species& species, const std::vector<double>& concentration)
  {
    auto var = variable_map_.find(species.name_);
    if (var == variable_map_.end())
      throw std::invalid_argument("Unknown variable '" + species.name_ + "'");
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
      throw std::invalid_argument("Unknown rate constant parameter '" + label + "'");
    if (custom_rate_parameters_.size() != 1)
      throw std::invalid_argument("Incorrect number of custom rate parameter values passed to multi-gridcell State");
    custom_rate_parameters_[0][param->second] = value;
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::SetCustomRateParameter(const std::string& label, const std::vector<double>& values)
  {
    auto param = custom_rate_parameter_map_.find(label);
    if (param == custom_rate_parameter_map_.end())
      throw std::invalid_argument("Unknown rate constant parameter '" + label + "'");
    if (custom_rate_parameters_.size() != values.size())
      throw std::invalid_argument("Incorrect number of custom rate parameter values passed to multi-gridcell State");
    for (std::size_t i = 0; i < custom_rate_parameters_.size(); ++i)
      custom_rate_parameters_[i][param->second] = values[i];
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::PrintHeader()
  {
    auto largest_str_iter = std::max_element(variable_names_.begin(), variable_names_.end(),
                                  [](const auto& a, const auto& b) {
                                      return a.size() < b.size();});
    int largest_str_size = largest_str_iter->size();
    int width = (largest_str_size < 10) ? 11 : largest_str_size + 2;

    std::cout << std::setw(5) << "time";
    for(auto& species : variable_names_)
    {
      std::cout << "," << std::setw(width) << species;
    }
    std::cout << std::endl;
  }

  template<template<class> class MatrixPolicy>
  inline void State<MatrixPolicy>::PrintState(double time)
  {
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);

    auto largest_str_iter = std::max_element(variable_names_.begin(), variable_names_.end(),
                                  [](const auto& a, const auto& b) {
                                      return a.size() < b.size();});
    int largest_str_size = largest_str_iter->size();
    int width = (largest_str_size < 10) ? 11 : largest_str_size + 2;

    std::cout << std::setw(5) << time << std::flush;
    for(auto& species : variable_names_)
    {
      std::cout << std::scientific << "," << std::setw(width) << std::setprecision(2) << variables_[0][variable_map_[species]];
    }
    std::cout << std::endl;
    std::cout.copyfmt(oldState);
  }

}  // namespace micm
