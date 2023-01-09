/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

namespace micm
{

  class IntraPhaseProcessBuilder
  {
   private:
   public:
    IntraPhaseProcessBuilder& For(const Phase& phase);
    IntraPhaseProcessBuilder& Reacting(const Species& reactant);
    IntraPhaseProcessBuilder& Producing(const Species& product);
    IntraPhaseProcessBuilder& WithYield(const double yield);
    IntraPhaseProcessBuilder& WithRateConstant(const RateConstant& rate_constant);
    IntraPhaseProcessBuilder& Build();
  };

}  // namespace micm