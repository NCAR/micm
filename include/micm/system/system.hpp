/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/condition.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <string>
#include <vector>

namespace micm
{

  class System
  {
   private:
    const Phase gas_phase_;
    const std::vector<Phase> phases_;
    const std::vector<Condition> conditions_;

   public:
    System();

    size_t Size();
    const Phase& FindPhase(const std::string& name);
    const Species& FindSpecies(const Phase& phase, const std::string& name);
  };

}  // namespace micm
