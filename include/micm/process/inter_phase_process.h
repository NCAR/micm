/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_INTER_PHASE_PROCESS_H
#define MICM_INTER_PHASE_PROCESS_H

#include <micm/process/process.h>
#include <micm/system/species.h>

#include <vector>

class InterPhaseProcess : public Process {
  private:
    const std::vector<const Species&> phase_1;
    const std::vector<const Species&> phase_2;

  public:
};

#endif