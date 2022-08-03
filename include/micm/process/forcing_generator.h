/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_FORCING_GENERATOR_H
#define MICM_FORCING_GENERATOR_H

#include <micm/jit.h>

class ForcingGenerator {
  public:
    virtual void GenerateForcing(JIT) = 0;
};

#endif