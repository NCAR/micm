/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_ARRHENIUS_RATE_CONSTANT_H
#define MICM_ARRHENIUS_RATE_CONSTANT_H

#include <micm/jit.h>
#include <micm/process/rate_constant.h>

template <typename T> class ArrheniusRateConstant : public RateConstant {
  private:
    const T A_;
    const T B_;
    const T C_;
    const T D_;
    const T E_;

  public:
};

#endif