/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <string>

namespace micm
{

  /**
   * @brief A value representing some aspect of a species
   *
   */
  class Property
  {
   public:
    /// @brief The name of this property
    std::string name_;
    /// @brief The units
    std::string units_;
    /// @brief The value of this property
    double value_;

   public:
    /// @brief Constructs a property
    /// @param name The name of this property
    /// @param units The units of the value
    /// @param value The value of the property
    Property(std::string name, std::string units, const double value)
        : name_(name),
          units_(units),
          value_(value)
    {
    }
  };

}  // namespace micm
