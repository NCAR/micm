/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_SPECIES_H
#define MICM_SPECIES_H

#include <micm/system/property.h>

#include <string>
#include <vector>

class Species {
  private:
    const std::string name_;
    const std::vector<Property<double>> properties_;

  public:
};

#endif