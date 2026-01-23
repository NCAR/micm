// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>

#include <string>
#include <system_error>

enum class MicmConstraintErrc
{
  InvalidEquilibriumConstant = MICM_CONSTRAINT_ERROR_CODE_INVALID_EQUILIBRIUM_CONSTANT,
  UnknownSpecies = MICM_CONSTRAINT_ERROR_CODE_UNKNOWN_SPECIES,
};

namespace std
{
  template<>
  struct is_error_code_enum<MicmConstraintErrc> : true_type
  {
  };
}  // namespace std

class MicmConstraintErrorCategory : public std::error_category
{
 public:
  const char* name() const noexcept override
  {
    return MICM_ERROR_CATEGORY_CONSTRAINT;
  }

  std::string message(int ev) const override
  {
    switch (static_cast<MicmConstraintErrc>(ev))
    {
      case MicmConstraintErrc::InvalidEquilibriumConstant: return "Equilibrium constant must be positive";
      case MicmConstraintErrc::UnknownSpecies: return "Unknown species in constraint";
      default: return "Unknown error";
    }
  }
};

inline const MicmConstraintErrorCategory& MicmConstraintError()
{
  static const MicmConstraintErrorCategory instance;
  return instance;
}

inline std::error_code make_error_code(MicmConstraintErrc e)
{
  return { static_cast<int>(e), MicmConstraintError() };
}
