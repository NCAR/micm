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
    std::unordered_map<std::string, std::string> others_{};
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
    /// @brief Tracks non-phase elements (e.g., number concentrations) associated with a model.
    ///        Elements are mapped using a prefix specific to the model's name and representation.
    std::unordered_map<std::string, std::string> others_;

    /// @brief Default constructor
    System() = default;

    /// @brief Parameterized constructor
    System(const Phase& gas_phase, const std::unordered_map<std::string, Phase>& phases, 
           const std::unordered_map<std::string, std::string>& others)
        : gas_phase_(gas_phase),
          phases_(phases),
          others_(others)
    {
    }

    /// @brief Parameterized constructor
    System(const Phase& gas_phase, const std::unordered_map<std::string, Phase>& phases)
        : gas_phase_(gas_phase),
          phases_(phases),
          others_()
    {
    }

    /// @brief Parameterized constructor with move semantics
    System(Phase&& gas_phase, std::unordered_map<std::string, Phase>&& phases,
          const std::unordered_map<std::string, std::string>& others)
        : gas_phase_(std::move(gas_phase)),
          phases_(std::move(phases)),
          others_(std::move(others))
    {
    }

    /// @brief Parameterized constructor with move semantics
    System(Phase&& gas_phase, std::unordered_map<std::string, Phase>&& phases)
        : gas_phase_(std::move(gas_phase)),
          phases_(std::move(phases)),
          others_()
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
          others_(parameters.others_)
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
    return state_size;
  }

  inline std::vector<std::string> System::UniqueNames() const
  {
    return UniqueNames(nullptr);
  }

  inline std::vector<std::string> System::UniqueNames(
      const std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> f) const
  {
    std::vector<std::string> names = gas_phase_.UniqueNames();
    for (auto& phase : phases_)
    {
      for (auto& species_name : phase.second.UniqueNames())
        names.push_back(phase.first + "." + species_name);
    }
    if (f)
    {
      const auto orig_names = names;
      for (std::size_t i = 0; i < orig_names.size(); ++i)
        names[i] = f(orig_names, i);
    }
    return names;
  }

}  // namespace micm
