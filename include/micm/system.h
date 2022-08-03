/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_SYSTEM_H
#define MICM_SYSTEM_H

#include <micm/aerosol_representation.h>
#include <micm/condition.h>
#include <micm/phase.h>
#include <micm/species.h>

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