// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>
#include <micm/util/micm_exception.hpp>

enum class MicmMatrixErrc
{
  RowSizeMismatch = MICM_MATRIX_ERROR_CODE_ROW_SIZE_MISMATCH,
  InvalidVector = MICM_MATRIX_ERROR_CODE_INVALID_VECTOR,
  ElementOutOfRange = MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE,
  MissingBlockIndex = MICM_MATRIX_ERROR_CODE_MISSING_BLOCK_INDEX,
  ZeroElementAccess = MICM_MATRIX_ERROR_CODE_ZERO_ELEMENT_ACCESS
};
