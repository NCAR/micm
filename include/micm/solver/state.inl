// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline State<MatrixPolicy, SparseMatrixPolicy>::State()
      : conditions_(),
        variables_(),
        custom_rate_parameters_(),
        rate_constants_()
  {
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline State<MatrixPolicy, SparseMatrixPolicy>::State(
      const std::size_t state_size,
      const std::size_t custom_parameters_size,
      const std::size_t process_size)
      : conditions_(1),
        variables_(1, state_size, 0.0),
        custom_rate_parameters_(1, custom_parameters_size, 0.0),
        rate_constants_(1, process_size, 0.0),
        jacobian_()
  {
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline State<MatrixPolicy, SparseMatrixPolicy>::State(const StateParameters& parameters)
      : conditions_(parameters.number_of_grid_cells_),
        variables_(parameters.number_of_grid_cells_, parameters.variable_names_.size(), 0.0),
        custom_rate_parameters_(parameters.number_of_grid_cells_, parameters.custom_rate_parameter_labels_.size(), 0.0),
        rate_constants_(parameters.number_of_grid_cells_, parameters.number_of_rate_constants_, 0.0),
        variable_map_(),
        custom_rate_parameter_map_(),
        variable_names_(parameters.variable_names_)
  {
    std::size_t index = 0;
    for (auto& name : variable_names_)
      variable_map_[name] = index++;
    index = 0;
    for (auto& label : parameters.custom_rate_parameter_labels_)
      custom_rate_parameter_map_[label] = index++;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::SetConcentrations(
      const std::unordered_map<std::string, std::vector<double>>& species_to_concentration)
  {
    const int num_grid_cells = conditions_.size();
    for (const auto& pair : species_to_concentration)
      SetConcentration({ pair.first }, pair.second);
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::SetConcentration(const Species& species, double concentration)
  {
    auto var = variable_map_.find(species.name_);
    if (var == variable_map_.end())
      throw std::invalid_argument("Unknown variable '" + species.name_ + "'");
    if (variables_.size() != 1)
      throw std::invalid_argument("Incorrect number of concentration values passed to multi-gridcell State");
    variables_[0][variable_map_[species.name_]] = concentration;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::SetConcentration(const Species& species, const std::vector<double>& concentration)
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

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::SetCustomRateParameters(
      const std::unordered_map<std::string, std::vector<double>>& parameters)
  {
    for (auto& pair : parameters)
      SetCustomRateParameter(pair.first, pair.second);
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::SetCustomRateParameter(const std::string& label, double value)
  {
    auto param = custom_rate_parameter_map_.find(label);
    if (param == custom_rate_parameter_map_.end())
      throw std::invalid_argument("Unkonwn rate constant parameter '" + label + "'");
    if (custom_rate_parameters_.size() != 1)
      throw std::invalid_argument("Incorrect number of custom rate parameter values passed to multi-gridcell State");
    custom_rate_parameters_[0][param->second] = value;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::SetCustomRateParameter(const std::string& label, const std::vector<double>& values)
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
