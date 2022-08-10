/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_INTRA_PHASE_PROCESS_H
#define MICM_INTRA_PHASE_PROCESS_H

#include <micm/process/forcing_generator.h>
#include <micm/process/jacobian_generator.h>
#include <micm/process/process.h>
#include <micm/process/rate_constant.h>
#include <micm/system/species.h>
#include <micm/jit.h>

#include <vector>

class IntraPhaseProcess : public Process, ForcingGenerator, JacobianGenerator {
  private:
    const std::vector<const Species&> reactants_;
    const std::vector<const Species&> products_;
    const RateConstant& rate_constant_;

  public:
    void GenerateForcing(JIT);
    void GenerateJacobian(JIT);
};

#endif