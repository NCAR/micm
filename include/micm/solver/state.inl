// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline State<MatrixPolicy, SparseMatrixPolicy>::State()
      : conditions_(),
        variables_(),
        custom_rate_parameters_(),
        rate_constants_(),
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
        variable_names_(parameters.variable_names_),
        jacobian_(),
        lower_matrix_(),
        upper_matrix_(),
        state_size_(parameters.variable_names_.size()),
        number_of_grid_cells_(parameters.number_of_grid_cells_)
  {
    std::size_t index = 0;
    for (auto& name : variable_names_)
      variable_map_[name] = index++;
    index = 0;
    for (auto& label : parameters.custom_rate_parameter_labels_)
      custom_rate_parameter_map_[label] = index++;

    jacobian_ = build_jacobian<SparseMatrixPolicy>(
      parameters.nonzero_jacobian_elements_,
      parameters.number_of_grid_cells_,
      state_size_
    );
    
    auto lu =  LuDecomposition::GetLUMatrices(jacobian_, 1.0e-30);
    auto lower_matrix = std::move(lu.first);
    auto upper_matrix = std::move(lu.second);
    lower_matrix_ = lower_matrix;
    upper_matrix_ = upper_matrix;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::SetConcentrations(
      const std::unordered_map<std::string, std::vector<double>>& species_to_concentration)
  {
    const std::size_t num_grid_cells = conditions_.size();
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
  inline void State<MatrixPolicy, SparseMatrixPolicy>::UnsafelySetCustomRateParameters(
      const std::vector<std::vector<double>>& parameters)
  {
    if (parameters.size() != variables_.size())
      throw std::invalid_argument("The number of grid cells configured for micm does not match the number of custom rate parameter values passed to multi-gridcell State");

    if (parameters[0].size() != custom_rate_parameters_[0].size())
      throw std::invalid_argument("The number of custom rate parameters configured for micm does not match the provided number of custom rate parameter values");

    for(size_t i = 0; i < number_of_grid_cells_; ++i) {
      custom_rate_parameters_[i] = parameters[i];
    }
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
      throw std::invalid_argument("Unknown rate constant parameter '" + label + "'");
    if (custom_rate_parameters_.size() != 1)
      throw std::invalid_argument("Incorrect number of custom rate parameter values passed to multi-gridcell State");
    custom_rate_parameters_[0][param->second] = value;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::SetCustomRateParameter(const std::string& label, const std::vector<double>& values)
  {
    auto param = custom_rate_parameter_map_.find(label);
    if (param == custom_rate_parameter_map_.end())
      throw std::invalid_argument("Unknown rate constant parameter '" + label + "'");
    if (custom_rate_parameters_.size() != values.size())
      throw std::invalid_argument("Incorrect number of custom rate parameter values passed to multi-gridcell State");
    for (std::size_t i = 0; i < custom_rate_parameters_.size(); ++i)
      custom_rate_parameters_[i][param->second] = values[i];
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::PrintHeader()
  {
    auto largest_str_iter = std::max_element(variable_names_.begin(), variable_names_.end(),
                                  [](const auto& a, const auto& b) {
                                      return a.size() < b.size();});
    int largest_str_size = largest_str_iter->size();
    int width = (largest_str_size < 10) ? 11 : largest_str_size + 2;

    std::cout << std::setw(6) << "time";
    if (variables_.size() > 1) {
      std::cout << "," << std::setw(6) << "grid";
    }

    for(const auto& [species, index] : variable_map_)
    {
      std::cout << "," << std::setw(width) << species;
    }
    std::cout << std::endl;
  }

  template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
  inline void State<MatrixPolicy, SparseMatrixPolicy>::PrintState(double time)
  {
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);

    auto largest_str_iter = std::max_element(variable_names_.begin(), variable_names_.end(),
                                  [](const auto& a, const auto& b) {
                                      return a.size() < b.size();});
    int largest_str_size = largest_str_iter->size();
    int width = (largest_str_size < 10) ? 11 : largest_str_size + 2;


    for(size_t i = 0; i < variables_.size(); ++i) {
      std::cout << std::setw(6) << time << ",";

      if (variables_.size() > 1) {
        std::cout << std::setw(6) << i << ",";
      }

      bool first = true;
      for(const auto& [species, index] : variable_map_)
      {
        if (!first) {
          std::cout << ",";
        }
        std::cout << std::scientific << std::setw(width) << std::setprecision(2) << variables_[i][index];
        first = false;
      }
      std::cout << std::endl;
      std::cout.copyfmt(oldState);
    }
  }

}  // namespace micm
