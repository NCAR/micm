/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <string>
#include <vector>

namespace micm
{
  struct SystemParameters
  {
    Phase gas_phase_{};
    std::vector<Phase> phases_{};
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
    const std::vector<Phase> phases_;

   public:
    /// @brief Default constructor
    System();

    /// @brief
    /// @param gas_phase
    System(SystemParameters parameters);

    /// @brief Returns the number of doubles required to store the system state
    size_t StateSize() const;
  };

  inline micm::System::System()
      : gas_phase_(),
        phases_()
  {
  }

  inline System::System(SystemParameters parameters)
      : gas_phase_(parameters.gas_phase_),
        phases_(parameters.phases_)
  {
  }

  inline size_t System::StateSize() const {
    size_t state_size = gas_phase_.species_.size();
    for (const auto& phase : phases_) {
      state_size += phase.species_.size();
    }
    return state_size; 
  }

}  // namespace micm
