/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_PHASE_H
#define MICM_PHASE_H

#include <micm/system/species.h>

#include <vector>

class Phase {
  private:
    const std::vector<Species> species_;

  public:
    std::size_t Size();
};

#endif