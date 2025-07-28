// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace micm
{
  struct SystemParameters
  {
    Phase gas_phase_{};
    std::unordered_map<std::string, Phase> phases_{};
    std::vector<std::string> custom_variables_{};
  };

  /// @brief Represents the complete chemical state of a grid cell
  ///        Includes the gas phase and other associated phases, each with their own set of species.
  class System
  {
   public:
    /// @brief The gas phase, defining a set of species present in the system
    Phase gas_phase_;
    /// @brief Additional phases (e.g., aqueous, aerosol), mapped by name and representing non-gas phase
    std::unordered_map<std::string, Phase> phases_;
    /// @brief Non-phase elements (e.g., number concentrations) that are time-varying and 
    ///        involved in the model's numerical solution.
    std::vector<std::string> custom_variables_;

    /// @brief Default constructor
    System() = default;

    /// @brief Parameterized constructor
    System(
        const Phase& gas_phase,
        const std::unordered_map<std::string, Phase>& phases,
        const std::vector<std::string>& custom_variables)
        : gas_phase_(gas_phase),
          phases_(phases),
          custom_variables_(custom_variables)
    {
    }

    /// @brief Parameterized constructor
    System(const Phase& gas_phase, const std::unordered_map<std::string, Phase>& phases)
        : gas_phase_(gas_phase),
          phases_(phases),
          custom_variables_()
    {
    }

    /// @brief Parameterized constructor with move semantics
    System(
        Phase&& gas_phase,
        std::unordered_map<std::string, Phase>&& phases,
        std::vector<std::string>&& custom_variables)
        : gas_phase_(std::move(gas_phase)),
          phases_(std::move(phases)),
          custom_variables_(std::move(custom_variables))
    {
    }

    /// @brief Parameterized constructor with move semantics
    System(Phase&& gas_phase, std::unordered_map<std::string, Phase>&& phases)
        : gas_phase_(std::move(gas_phase)),
          phases_(std::move(phases)),
          custom_variables_()
    {
    }

    /// @brief Copy constructor
    System(const System&) = default;

    /// @brief Move constructor
    System(System&&) = default;

    /// @brief Constructor from SystemParameters
    System(const SystemParameters& parameters)
        : gas_phase_(parameters.gas_phase_),
          phases_(parameters.phases_),
          custom_variables_(parameters.custom_variables_)
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
    size_t state_size = gas_phase_.StateSize();
    for (const auto& phase : phases_)
    {
      state_size += phase.second.StateSize();
    }
    state_size += custom_variables_.size();

    return state_size;
  }

  inline std::vector<std::string> System::UniqueNames() const
  {
    return UniqueNames(nullptr);
  }

  inline std::vector<std::string> System::UniqueNames(
      const std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> f) const
  {
    std::size_t num_elements = gas_phase_.species_.size();
    for (const auto& phase : phases_)
    {
      num_elements += phase.second.species_.size();
    }
    num_elements += custom_variables_.size();

    // Collect all names: gas phase, phase species, and custom variables
    std::vector<std::string> names;
    names.reserve(num_elements);
  
    // Add gas phase species
    std::vector<std::string> gas_names = gas_phase_.UniqueNames();
    names.insert(names.end(), gas_names.begin(), gas_names.end());

    // Add phase species with phase prefix
    for (const auto& phase : phases_)
    {
      for (const auto& species_name : phase.second.UniqueNames())
        names.push_back(phase.first + "." + species_name);
    }

    // Add custom variables
    names.insert(names.end(), custom_variables_.begin(), custom_variables_.end());

    // Apply custom ordering or transformation if provided
    if (f)
    {
      const auto orig_names = names;
      for (std::size_t i = 0; i < orig_names.size(); ++i)
        names[i] = f(orig_names, i);
    }
    return names;
  }

}  // namespace micm
