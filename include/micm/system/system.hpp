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
    System(const SystemParameters& parameters);
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

}  // namespace micm
