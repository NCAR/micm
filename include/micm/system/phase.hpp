/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/species.hpp>
#include <vector>

namespace micm
{

  /**
   * @brief A class which represnts a certain chemical phase (gaseous, aquous)
   *
   * Each phase represents a set of species that participate in chemical reactions in that phase.
   */
  class Phase
  {
   public:
    /// @brief The list of species
    std::vector<Species> species_;

   public:
    /// @brief Default constructor
    Phase() = default;

    /// @brief Create a phase with a set of species
    /// @param species A unique list of species
    Phase(std::vector<Species> species);
  };

  Phase::Phase(std::vector<Species> species)
      : species_(std::move(species))
  {
  }

}  // namespace micm