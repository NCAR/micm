/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <string>

namespace micm
{

  /**
   * @brief An environemental condition
   *
   */
  class Condition
  {
   public:
    /// @brief The name of this condition
    const std::string name_;
    /// @brief The units of this condition
    const std::string units_;

   public:
    /// @brief Define an environmental condition
    /// @param name The name
    /// @param units The units
    Condition(std::string name, std::string units)
      : name_(std::move(name)),
        units_(std::move(units)){};
  };
}  // namespace micm
