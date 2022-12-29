/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

namespace micm
{

  class IntraPhaseProcessBuilder
  {
   private:
   public:
    IntraPhaseProcessBuilder& For(const Phase& phase);
    IntraPhaseProcessBuilder& With(const Species& phase);
    IntraPhaseProcessBuilder& Build();
  };

}  // namespace micm