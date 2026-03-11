// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>
#include <micm/util/micm_exception.hpp>

enum class MicmConstraintErrc
{
  InvalidEquilibriumConstant = MICM_CONSTRAINT_ERROR_CODE_INVALID_EQUILIBRIUM_CONSTANT,
  UnknownSpecies = MICM_CONSTRAINT_ERROR_CODE_UNKNOWN_SPECIES,
  EmptyReactants = MICM_CONSTRAINT_ERROR_CODE_EMPTY_REACTANTS,
  EmptyProducts = MICM_CONSTRAINT_ERROR_CODE_EMPTY_PRODUCTS,
  InvalidStoichiometry = MICM_CONSTRAINT_ERROR_CODE_INVALID_STOICHIOMETRY,
};
