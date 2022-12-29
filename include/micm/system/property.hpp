/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <string>

namespace micm
{

  template<typename T>
  class Property
  {
   private:
    const std::string name_;
    const std::string units_;
    const T value_;

   public:
    Property(const std::string& name, const std::string& units, const T value);
  };

}  // namespace micm
