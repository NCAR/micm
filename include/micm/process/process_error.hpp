// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>

#include <string>
#include <system_error>

enum class MicmProcessErrc
{
  ReactantDoesNotExist = MICM_PROCESS_ERROR_CODE_REACTANT_DOES_NOT_EXIST,
  ProductDoesNotExist  = MICM_PROCESS_ERROR_CODE_PRODUCT_DOES_NOT_EXIST,
  RateConstantIsNotSet = MICM_PROCESS_ERROR_CODE_RATE_CONSTANT_IS_NOT_SET,
  PhaseIsNotSet        = MICM_PROCESS_ERROR_CODE_PHASE_IS_NOT_SET,
};

namespace std
{
  template<>
  struct is_error_code_enum<MicmProcessErrc> : true_type
  {
  };
}  // namespace std


class MicmProcessErrorCategory : public std::error_category
{
  public:
  const char *name() const noexcept override
  {
    return MICM_ERROR_CATEGORY_PROCESS;
  }

  std::string message(int ev) const override
  {
    switch (static_cast<MicmProcessErrc>(ev))
    {
      case MicmProcessErrc::ReactantDoesNotExist: return "Reactant does not exist";
      case MicmProcessErrc::ProductDoesNotExist:  return "Product does not exist";
      case MicmProcessErrc::RateConstantIsNotSet: return "Rate constant is not set";
      case MicmProcessErrc::PhaseIsNotSet:        return "Phase is not set";
      default: return "Unknown error";
    }
  }
};

inline const MicmProcessErrorCategory& MicmProcessError()
{
  static const MicmProcessErrorCategory instance;
  return instance;
}

inline std::error_code make_error_code(MicmProcessErrc e)
{
  return { static_cast<int>(e), MicmProcessError() };
}