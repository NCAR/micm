// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>

#include <string>
#include <system_error>

enum class MicmSolverErrc
{
  UnusedSpecies = MICM_SOLVER_ERROR_CODE_UNUSED_SPECIES,
  MissingChemicalSystem = MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SYSTEM,
  MissingProcesses = MICM_SOLVER_ERROR_CODE_MISSING_PROCESSES,
  MissingChemicalSpecies = MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SPECIES,
  InvalidToleranceSize = MICM_SOLVER_ERROR_CODE_INVALID_TOLERANCE_SIZE,
};

namespace std
{
  template<>
  struct is_error_code_enum<MicmSolverErrc> : true_type
  {
  };
}  // namespace std

class MicmSolverErrorCategory : public std::error_category
{
 public:
  const char* name() const noexcept override
  {
    return MICM_ERROR_CATEGORY_SOLVER;
  }

  std::string message(int ev) const override
  {
    switch (static_cast<MicmSolverErrc>(ev))
    {
      case MicmSolverErrc::UnusedSpecies:
        return "Unused species present in the chemical system. Use the ignore_unused_species_ parameter to allow unused "
                "species in the solve.";
      case MicmSolverErrc::MissingChemicalSystem:
        return "Missing chemical system. Use the SetSystem function to set the chemical system.";
      case MicmSolverErrc::MissingProcesses:
        return "Missing processes. Use the SetReactions function to set the processes.";
      case MicmSolverErrc::MissingChemicalSpecies: return "Provided chemical system contains no species.";
      case MicmSolverErrc::InvalidToleranceSize:
        return "Provided tolerances do not match the number of species in the chemical system. Either provide none and "
                "allow defaults to be set or pass in a number equal to the number of chemical species.";
      default: return "Unknown error";
    }
  }
};

inline const MicmSolverErrorCategory& MicmSolverError()
{
  static const MicmSolverErrorCategory instance;
  return instance;
}

inline std::error_code make_error_code(MicmSolverErrc e)
{
  return { static_cast<int>(e), MicmSolverError() };
}