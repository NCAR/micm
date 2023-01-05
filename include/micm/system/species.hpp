/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/property.hpp>
#include <string>
#include <vector>

namespace micm
{

  /**
   * @brief A representation of a chemcial species
   * 
   */
  template<typename T>
  class Species
  {
   private:
    /// @brief The name of this species
    const std::string name_;
    /// @brief A list of properties of this species
    const std::vector<Property<T>> properties_;

   public:
    /// @brief Default constructor
    Species() = default;
    /// @brief Construct a species by name only
    /// @param name The name of the species
    Species(std::string name);
    /// @brief 
    /// @param name The name of the species
    /// @param properties The properties of teh species
    Species(std::string name, std::vector<Property<T>> properties);
  };

}  // namespace micm
