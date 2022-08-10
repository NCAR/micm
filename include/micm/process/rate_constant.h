/* Copyright (C) 2022 National Center for Atmospheric Research,
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