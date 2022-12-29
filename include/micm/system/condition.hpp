/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <string>

namespace micm
{

  class Condition
  {
   private:
   public:
    const std::string name_;
    const std::string units_;

   public:
    Condition();
    Condition(std::string name, std::string units);
  };

}  // namespace micm
