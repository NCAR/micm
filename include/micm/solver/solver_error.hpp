// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>
#include <micm/util/micm_exception.hpp>

enum class MicmSolverErrc
{
  UnusedSpecies = MICM_SOLVER_ERROR_CODE_UNUSED_SPECIES,
  MissingChemicalSystem = MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SYSTEM,
  MissingProcesses = MICM_SOLVER_ERROR_CODE_MISSING_PROCESSES,
  MissingChemicalSpecies = MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SPECIES,
  InvalidToleranceSize = MICM_SOLVER_ERROR_CODE_INVALID_TOLERANCE_SIZE,
};
