/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <string>
#include <vector>
#include <unordered_map>

namespace micm
{
  struct SystemParameters
  {
    Phase gas_phase_{};
    std::unordered_map<std::string, Phase> phases_{};
  };

  /**
   * @brief A `System` holds all physical information that represents a grid cell.
   *
   */
  class System
  {
   public:
    /// @brief The gas phase is a micm::Phase and determines what species are present.
    const Phase gas_phase_;
    /// @brief This is a catchall for anything that is not the gas phase.
    const std::unordered_map<std::string, Phase> phases_;

   public:
    /// @brief Default constructor
    System();

    /// @brief
    /// @param gas_phase
    System(const SystemParameters& parameters);

    /// @brief Returns the number of doubles required to store the system state
    size_t StateSize() const;

    /// @brief Returns a set of unique species names
    std::vector<std::string> UniqueNames() const;
  };

  inline micm::System::System()
      : gas_phase_(),
        phases_()
  {
  }

  inline System::System(const SystemParameters& parameters)
      : gas_phase_(parameters.gas_phase_),
        phases_(parameters.phases_)
  {
  }

  inline size_t System::StateSize() const
  {
    size_t state_size = gas_phase_.species_.size();
    for (const auto& phase : phases_)
    {
      state_size += phase.second.species_.size();
    }
    return state_size;
  }

  inline std::vector<std::string> System::UniqueNames() const {
    std::vector<std::string> names{};
    for (auto& species: gas_phase_.species_) {
      names.push_back(species.name_);
    }
    for (auto& phase: phases_) {
      for (auto &species: phase.second.species_) {
        names.push_back(phase.first + "." + species.name_);
      }
    }
    return names;
  }

}  // namespace micm
