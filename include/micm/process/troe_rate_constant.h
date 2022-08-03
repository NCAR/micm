/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_TROE_RATE_CONSTANT_H
#define MICM_TROE_RATE_CONSTANT_H

#include <micm/jit.h>
#include <micm/process/rate_constant.h>

template <typename T> class TroeRateConstant : public RateConstant {
  private:
    const T k_0;
    const T k_inf;
    const T Fc_;
    const T N_;

  public:
};

#endif