/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_RATE_CONSTANT_H
#define MICM_RATE_CONSTANT_H

#include <micm/jit.h>

class RateConstant {
  public:
    virtual void Generate(JIT jit) = 0;
};

#endif