/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_HENRYS_LAW_INTER_PHASE_PROCESS_H
#define MICM_HENRYS_LAW_INTER_PHASE_PROCESS_H

#include <micm/process/forcing_generator.h>
#include <micm/process/inter_phase_process.h>
#include <micm/process/jacobian_generator.h>
#include <micm/process/rate_constant.h>
#include <micm/system/species.h>

#include <vector>

class HenrysLawInterPhaseProcess : public InterPhaseProcess,
                                                    ForcingGenerator, JacobianGenerator {
  private:
    const RateConstant& forward_rate_constant_;
    const RateConstant& equilibrium_rate_constant_;

  public:
};

#endif