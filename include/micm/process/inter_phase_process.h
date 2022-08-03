/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_INTER_PHASE_PROCESS_H
#define MICM_INTER_PHASE_PROCESS_H

#include <micm/system/species.h>

#include <vector>

class InterPhaseProcess {
  private:
    const std::vector<const Species&> phase_1;
    const std::vector<const Species&> phase_2;

  public:
};

#endif