/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/condition.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <string>
#include <vector>

namespace micm
{

  /**
   * @brief A `System` holds all physical information that represents a grid cell.
   *
   */
  class System
  {
   private:
    /// @brief The gas phase is a micm::Phase and determines what species are present.
    const Phase gas_phase_;
    /// @brief This is a catchall for anything that is not the gas phase.
    const std::vector<Phase> phases_;
    /// @brief This represents any physical measurement of a grid cell.
    const std::vector<Condition> conditions_;

   public:
    /// @brief Default constructor
    System();

    /// @brief The size of the system
    /// @return The number of reacting species
    std::size_t Size();
    /// @brief Locate a particular phase held by the system
    /// @param name The identifier for a phase
    /// @return A non-owning pointer to a micm::Phase or `nullptr`.
    const Phase* FindPhase(const std::string& name);
    /// @brief Locate a species inside of a phase
    /// @param phase The identifier for a phase
    /// @param name The name of a particular species in a particular phase.
    /// @return A non-owning pointer to a micm::Phase or `nullptr`.
    const Species* FindSpecies(const Phase& phase, const std::string& name);
  };

  inline micm::System::System()
      : gas_phase_(),
        phases_(),
        conditions_()
  {
  }

}  // namespace micm
