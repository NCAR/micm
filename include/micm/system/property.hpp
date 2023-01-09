/* Copyright (C) 2022 National Center for Atmospheric Research,
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
   * @tparam T the underlying datatype of the value
   */
  template<typename T>
  class Property
  {
   private:
    /// @brief The name of the species
    const std::string name_;
    /// @brief The units
    const std::string units_;
    /// @brief The value of this property
    const T value_;

   public:
    /// @brief Constructs a property
    /// @param name The name of the species
    /// @param units The units of the value
    /// @param value The value of the property
    Property(const std::string& name, const std::string& units, const T value);
  };

}  // namespace micm
