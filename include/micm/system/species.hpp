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

  class Species
  {
   private:
    const std::string name_;
    const std::vector<Property<double>> properties_;

   public:
    Species();
    Species(std::string name);
    Species(std::string name, std::vector<Property<double>> properties);
  };

}  // namespace micm
