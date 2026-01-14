// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>

#include <string>
#include <system_error>

enum class MicmMatrixErrc
{
  RowSizeMismatch = MICM_MATRIX_ERROR_CODE_ROW_SIZE_MISMATCH,
  InvalidVector = MICM_MATRIX_ERROR_CODE_INVALID_VECTOR,
  ElementOutOfRange = MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE,
  MissingBlockIndex = MICM_MATRIX_ERROR_CODE_MISSING_BLOCK_INDEX,
  ZeroElementAccess = MICM_MATRIX_ERROR_CODE_ZERO_ELEMENT_ACCESS
};

namespace std
{
  template<>
  struct is_error_code_enum<MicmMatrixErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class MicmMatrixErrorCategory : public std::error_category
  {
   public:
    const char *name() const noexcept override
    {
      return MICM_ERROR_CATEGORY_MATRIX;
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmMatrixErrc>(ev))
      {
        case MicmMatrixErrc::RowSizeMismatch: return "Matrix row size mismatch in assignment from vector";
        case MicmMatrixErrc::InvalidVector: return "Invalid vector for matrix assignment";
        case MicmMatrixErrc::ElementOutOfRange: return "Element out of range";
        case MicmMatrixErrc::MissingBlockIndex: return "Missing block index";
        case MicmMatrixErrc::ZeroElementAccess: return "Zero element access";
        default: return "Unknown error";
      }
    }
  };

  const MicmMatrixErrorCategory MICM_MATRIX_ERROR{};
}  // namespace

inline std::error_code make_error_code(MicmMatrixErrc e)
{
  return { static_cast<int>(e), MICM_MATRIX_ERROR };
}