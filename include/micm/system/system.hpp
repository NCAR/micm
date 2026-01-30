// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/utils.hpp>

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace micm
{
  // Helper class for collecting lambda functions from external models
  struct ExternalModel
  {
    std::function<size_t()> state_size_func_;
    std::function<std::vector<std::string>()> unique_names_func_;
    
    template<typename ModelType>
    ExternalModel(const ModelType& model)
    : state_size_func_([model]() { return model.StateSize(); }),
      unique_names_func_([model]() { return model.UniqueNames(); })
    {
    }
  };
  
  struct SystemParameters
  {
    Phase gas_phase_{};
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
        ExternalModels... external_models)
        : gas_phase_(gas_phase),
          external_models_{ ExternalModel{ external_models }... }
    {
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
      state_size += model.state_size_func_();
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
    auto gas_names = gas_phase_.UniqueSpeciesNames();
    names.insert(names.end(), std::make_move_iterator(gas_names.begin()), std::make_move_iterator(gas_names.end()));

    // Include names from external models
    for (const auto& model : external_models_)
    {
      auto model_names = model.unique_names_func_();
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
