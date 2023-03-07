/* Copyright (C) 2023 National Center for Atmospheric Research,
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
  class Species
  {
   public:
    /// @brief The name of this species
    const std::string name_;
    /// @brief A list of properties of this species
    const std::vector<Property> properties_;

    /// @brief Default constructor
    Species() = default;
    /// @brief Construct a species by name only
    /// @param name The name of the species
    Species(const std::string name) : name_(name) {};
    /// @brief
    /// @param name The name of the species
    /// @param properties The properties of the species
    Species(const std::string name, const std::vector<Property> properties) :
      name_(name), properties_(properties) {};
  };

}  // namespace micm
