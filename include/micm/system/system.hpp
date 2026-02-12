// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/utils.hpp>

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>

namespace micm
{
  // Helper class for collecting lambda functions from external models
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
  
  struct SystemParameters
  {
    /// @brief  @brief The gas phase
    Phase gas_phase_{};
    /// @brief External models (e.g., aerosol models) that provide additional components to the system
    std::vector<ExternalModel> external_models_{};
  };

  /// @brief Represents the complete chemical state of a grid cell
  ///        Includes the gas phase and other associated phases, each with their own set of species.
  class System
  {
   public:
    /// @brief The gas phase, defining a set of species present in the system
    Phase gas_phase_;
    /// @brief External models (e.g., aerosol models) that provide additional components to the system
    std::vector<ExternalModel> external_models_;

    /// @brief Default constructor
    System() = default;

    /// @brief Parameterized constructor
    System(
        const Phase& gas_phase)
        : gas_phase_(gas_phase),
          external_models_()
    {
    }

    /// @brief Constructor with external models
    template<typename... ExternalModels>
    System(
        const Phase& gas_phase,
        ExternalModels&&... external_models)
        : gas_phase_(gas_phase),
          external_models_{ ExternalModel{ std::forward<ExternalModels>(external_models) }... }
    {
      if (StateSize() != UniqueNames().size())
      {
        throw std::invalid_argument("Mismatch between system state size and number of unique names. Likely duplicate species names.");
      }
    }

    /// @brief Copy constructor
    System(const System&) = default;

    /// @brief Move constructor
    System(System&&) = default;

    /// @brief Constructor from SystemParameters
    System(const SystemParameters& parameters)
        : gas_phase_(parameters.gas_phase_),
          external_models_(parameters.external_models_)
    {
      if (StateSize() != UniqueNames().size())
      {
        throw std::invalid_argument("Mismatch between system state size and number of unique names. Likely duplicate species names.");
      }
    }

    /// @brief Copy assignment operator
    System& operator=(const System&) = default;

    /// @brief Move assignment operator
    System& operator=(System&&) = default;

    /// @brief Returns the number of doubles required to store the system state
    size_t StateSize() const;

    /// @brief Returns a set of unique species names
    /// @return vector of unique state variable names
    std::vector<std::string> UniqueNames() const;

    /// @brief Returns a set of unique species names
    /// @param f Function used to apply specific order to unique names
    /// @return vector of unique state variable names
    std::vector<std::string> UniqueNames(
        const std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> f) const;
  };

  inline size_t System::StateSize() const
  {
    std::size_t state_size = gas_phase_.StateSize();
    for (const auto& model : external_models_)
    {
      state_size += std::get<0>(model.state_size_func_());
    }

    return state_size;
  }

  inline std::vector<std::string> System::UniqueNames() const
  {
    return UniqueNames(nullptr);
  }

  inline std::vector<std::string> System::UniqueNames(
      const std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> f) const
  {
    std::vector<std::string> names;
    names.reserve(StateSize());

    // Exclude phase name for gas phase species to maintain consistency with prior behavior
    // e.g., "O3" instead of "GAS.O3"
    auto gas_names = gas_phase_.SpeciesNames();
    names.insert(names.end(), std::make_move_iterator(gas_names.begin()), std::make_move_iterator(gas_names.end()));

    // Include names from external models
    for (const auto& model : external_models_)
    {
      auto model_names = model.variable_names_func_();
      names.insert(names.end(), model_names.begin(), model_names.end());
    }

    if (f)
    {
      std::vector<std::string> reordered;
      reordered.reserve(names.size());
      for (std::size_t i = 0; i < names.size(); ++i)
        reordered.push_back(f(names, i));
      return reordered;
    }

    return names;
  }

}  // namespace micm
