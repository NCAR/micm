/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_INTER_PHASE_PROCESS_BUILDER_H
#define MICM_INTER_PHASE_PROCESS_BUILDER_H


#include <micm/process/inter_phase_process.h>
#include <micm/system/phase.h>
#include <micm/system/species.h>

class InterPhaseProcessBuilder {
  private:
    const InterPhaseProcess process_;
    const int current_phase_;

  public:
    InterPhaseProcessBuilder();

    InterPhaseProcessBuilder For(const Phase& phase);
    InterPhaseProcessBuilder With(const Species& phase);
    InterPhaseProcess Build();
};

#endif