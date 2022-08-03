/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_HENRYS_LAW_INTER_PHASE_PROCESS_BUILDER_H
#define MICM_HENRYS_LAW_INTER_PHASE_PROCESS_BUILDER_H

#include <micm/process/inter_phase_process_builder.h>
#include <micm/process/rate_constant.h>

#include <vector>

class HenrysLawInterPhaseProcessBuilder : public InterPhaseProcessBuilder {
  public:
    InterPhaseProcessBuilder ForwardRateConstant(const RateConstant& rate_constant);
    InterPhaseProcessBuilder EquilibriumRateConstant(const RateConstant& rate_constant);
    InterPhaseProcessBuilder Build();
};

#endif