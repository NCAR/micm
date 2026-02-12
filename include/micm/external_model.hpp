// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <functional>
#include <memory>
#include <string>
#include <set>

namespace micm
{
  // @brief Represents an external model (e.g., aerosol model) that can be integrated into the MICM system
  struct ExternalModel
  {
    /// @brief Function to get the state size (variables, parameters) of the external model
    std::function<std::tuple<std::size_t, std::size_t>()> state_size_func_;
    /// @brief Function to get the state variable names of the external model
    std::function<std::vector<std::string>()> variable_names_func_;
    /// @brief Function to get the state parameter names of the external model
    std::function<std::vector<std::string>()> parameter_names_func_;
    /// @brief Function to get the non-zero Jacobian elements of the external model
    std::function<std::set<std::pair<std::size_t, std::size_t>>(const std::unordered_map<std::string, std::size_t>&)> non_zero_jacobian_elements_func_;
    
    // Default constructor
    ExternalModel() = default;
    
    // Copy constructor
    ExternalModel(const ExternalModel&) = default;
    
    // Move constructor  
    ExternalModel(ExternalModel&&) = default;
    
    // Copy assignment
    ExternalModel& operator=(const ExternalModel&) = default;
    
    // Move assignment
    ExternalModel& operator=(ExternalModel&&) = default;
    
    /// @brief Constructor from an external model instance
    /// @tparam ModelType Type of the external model
    /// @param model Instance of the external model
    template<typename ModelType,
             typename = std::enable_if_t<!std::is_same_v<std::decay_t<ModelType>, ExternalModel>>>
    ExternalModel(ModelType&& model)
    {
      auto shared_model = std::make_shared<std::decay_t<ModelType>>(std::forward<ModelType>(model));
      state_size_func_ = [shared_model]() { return shared_model->StateSize(); };
      variable_names_func_ = [shared_model]() { return shared_model->StateVariableNames(); };
      parameter_names_func_ = [shared_model]() { return shared_model->StateParameterNames(); };
      non_zero_jacobian_elements_func_ = [shared_model](const std::unordered_map<std::string, std::size_t>& species_map) {
        return shared_model->NonZeroJacobianElements(species_map);
      };
    }
  };
  
}