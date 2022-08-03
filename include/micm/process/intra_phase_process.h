/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_INTRA_PHASE_PROCESS_BUILDER_H
#define MICM_INTRA_PHASE_PROCESS_BUILDER_H

#include <micm/process/forcing_generator.h>
#include <micm/process/jacobian_generator.h>
#include <micm/process/rate_constant.h>
#include <micm/system/species.h>

#include <vector>

class IntraPhaseProcessBuilder : public ForcingGenerator, JacobianGenerator {
  private:
    const std::vector<const Species&> reactants_;
    const std::vector<const Species&> products_;
    const RateConstant& rate_constant_;

  public:
};

#endif