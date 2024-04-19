// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <system_error>

enum class MicmErrc
{
  InvalidMechanism = 1, // Inconsistent or incomplete mechanism specification
  InvalidConfiguration = 2, // Invalid configuration data or format
  MissingSpecies = 3, // Missing species in mechanism
  ArraySizeMismatch = 4, // Array size mismatch
  InvalidArrayElementAccess = 5, // Invalid array element access
  InvalidDataType = 6, // Invalid data type for operation
  CUBLASInitializationFailure = 7, // cuBLAS initialization failed
  CUBLASOperationFailure = 8, // cuBLAS operation failed
};

std::error_condition make_error_condition(MicmErrc e);

namespace std
{
  template <>
  struct is_error_condition_enum<MicmErrc> : true_type
  {
  };
} // namespace std

namespace {

  class MicmErrorCategory : public std::error_category
  {
  public:
    const char* name() const noexcept override { return "MICM"; }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmErrc>(ev))
      {
      case MicmErrc::InvalidMechanism:
        return "Invalid mechanism";
      case MicmErrc::InvalidConfiguration:
        return "Invalid configuration";
      case MicmErrc::MissingSpecies:
        return "Missing species";
      case MicmErrc::ArraySizeMismatch:
        return "Array size mismatch";
      case MicmErrc::InvalidArrayElementAccess:
        return "Invalid array element access";
      case MicmErrc::InvalidDataType:
        return "Invalid data type";
      case MicmErrc::CUBLASInitializationFailure:
        return "cuBLAS initialization failure";
      case MicmErrc::CUBLASOperationFailure:
        return "cuBLAS operation failure";
      default:
        return "Unknown error";
      }
    }
  };

  const MicmErrorCategory micmErrorCategory{};

} // namespace

std::error_code make_error_code(MicmErrc e)
{
  return {static_cast<int>(e), micmErrorCategory};
}