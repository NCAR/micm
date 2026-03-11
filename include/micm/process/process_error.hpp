// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>
#include <micm/util/micm_exception.hpp>

enum class MicmProcessErrc
{
  ReactantDoesNotExist = MICM_PROCESS_ERROR_CODE_REACTANT_DOES_NOT_EXIST,
  ProductDoesNotExist = MICM_PROCESS_ERROR_CODE_PRODUCT_DOES_NOT_EXIST,
  RateConstantIsNotSet = MICM_PROCESS_ERROR_CODE_RATE_CONSTANT_IS_NOT_SET,
  TransferCoefficientIsNotSet = MICM_PROCESS_ERROR_CODE_TRANSFER_COEFFICIENT_IS_NOT_SET,
  InvalidConfiguration = MICM_PROCESS_ERROR_CODE_INVALID_CONFIGURATION,
};
