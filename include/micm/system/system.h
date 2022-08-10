/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_SYSTEM_H
#define MICM_SYSTEM_H

#include <micm/system/aerosol_representation.h>
#include <micm/system/condition.h>
#include <micm/system/phase.h>
#include <micm/system/species.h>

#include <vector>
#include <string>

class System {
  private:
    const Phase gas_phase_;
    const std::vector<AerosolRepresentation> aerosols_;
    const std::vector<Phase> phases_;
    const std::vector<Condition> conditions_;

  public:
    std::size_t Size();
    const Phase& FindPhase(const std::string& name);
    const Species& FindSpecies(const Phase& phase, const std::string& name);
};

#endif