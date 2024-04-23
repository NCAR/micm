// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <string>
#include <system_error>
#include <micm/util/error.hpp>

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
  template <>
  struct is_error_condition_enum<MicmMatrixErrc> : true_type
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
        case MicmMatrixErrc::RowSizeMismatch:
          return "Matrix row size mismatch in assignment from vector";
        case MicmMatrixErrc::InvalidVector:
          return "Invalid vector for matrix assignment";
        case MicmMatrixErrc::ElementOutOfRange:
          return "Element out of range";
        case MicmMatrixErrc::MissingBlockIndex:
          return "Missing block index";
        case MicmMatrixErrc::ZeroElementAccess:
          return "Zero element access";
        default:
          return "Unknown error";
      }
    }
  };

  const MicmMatrixErrorCategory micmMatrixErrorCategory{};
}  // namespace

std::error_code make_error_code(MicmMatrixErrc e)
{
  return {static_cast<int>(e), micmMatrixErrorCategory};
}