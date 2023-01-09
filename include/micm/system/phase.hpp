/* Copyright (C) 2022 National Center for Atmospheric Research,
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
  template<typename T>
  class Phase
  {
   private:
    /// @brief The list of species
    const std::vector<Species<T>> species_;

   public:
    /// @brief Default constructor
    Phase() = default;
    /// @brief Phase copy constructor
    /// @param other Another micm::Phase
    Phase(const Phase& other);
    /// @brief Create a phase with a set of species
    /// @param species A unique list of species
    Phase(std::vector<Species<T>> species);

    /// @brief Assignment operator. This assumes ownership of the other's list of species
    /// @param other Another micm::Phase
    /// @return a reference to the created phase
    Phase& operator=(const Phase& other);

    /// @brief The size of this phase
    /// @return The number of participating species
    size_t Size();
  };

}  // namespace micm