/* Copyright (C) 2022 National Center for Atmospheric Research,
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
   private:
   public:
    /// @brief The name of this condition
    const std::string name_;
    /// @brief The units of this condition
    const std::string units_;

   public:
    /// @brief Default constructor
    Condition();
    /// @brief Define an environmental condition
    /// @param name The name
    /// @param units The units
    Condition(std::string name, std::string units);
  };

}  // namespace micm
