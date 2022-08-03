/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_PHASE_H
#define MICM_PHASE_H

#include <micm/species.h>

#include <vector>

class Phase {
  private:
    const std::vector<Species> species_;

  public:
    std::size_t Size();
};

#endif