/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_INTRA_PHASE_PROCESS_BUILDER_H
#define MICM_INTRA_PHASE_PROCESS_BUILDER_H

#include <micm/process/intra_phase_process.h>
#include <micm/system/phase.h>
#include <micm/system/species.h>

class IntraPhaseProcessBuilder {
  private:
    const IntraPhaseProcess process_;

  public:
    IntraPhaseProcessBuilder();

    IntraPhaseProcessBuilder For(const Phase& phase);
    IntraPhaseProcessBuilder With(const Species& phase);
    IntraPhaseProcess Build();
};

#endif