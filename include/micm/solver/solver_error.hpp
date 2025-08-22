// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/error.hpp>

#include <string>
#include <system_error>

enum class MicmSolverErrc
{
  UnusedSpecies = MICM_SOLVER_ERROR_CODE_UNUSED_SPECIES,                     // Unused species in system
  MissingChemicalSystem = MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SYSTEM,    // Missing chemical system
  MissingProcesses = MICM_SOLVER_ERROR_CODE_MISSING_PROCESSES,               // Missing processes
  MissingChemicalSpecies = MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SPECIES,  // Missing chemical species
  InvalidToleranceSize = MICM_SOLVER_ERROR_CODE_INVALID_TOLERANCE_SIZE,      // Invalid tolerance size
  FailedToConverge = MICM_SOLVER_ERROR_CODE_FAILED_TO_CONVERGE,

};

MICM_SOLVER_ERROR_CODE_UNUSED_SPECIES = 1,                // Unused species present in the chemical system
    MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SYSTEM = 2,   // Missing chemical system
    MICM_SOLVER_ERROR_CODE_MISSING_PROCESSES = 3,         // Missing processes
    MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SPECIES = 4,  // Missing chemical species
    MICM_SOLVER_ERROR_CODE_INVALID_TOLERANCE_SIZE = 5     // Invalid tolerance size
    MICM_SOLVER_ERROR_CODE_FAILED_TO_CONVERGE = 1 MICM_SOLVER_ERROR_CODE_UNKNOWN_SPECIES = 1,  // Unknown species
    MICM_SOLVER_ERROR_CODE_UNKNOWN_RATE_CONSTANT_PARAMETER = 2,  // Unknown rate constant parameter
    MICM_SOLVER_ERROR_CODE_INCORRECT_NUMBER_OF_CONCENTRATION_VALUES_FOR_MULTI_GRIDCELL_STATE =
        3,  // Incorrect number of concentration values
    MICM_SOLVER_ERROR_CODE_INCORRECT_NUMBER_OF_CUSTOM_RATE_PARAMETER_VALUES =
        4,  // Incorrect number of custom rate parameter values
    MICM_SOLVER_ERROR_CODE_INCORRECT_NUMBER_OF_CUSTOM_RATE_PARAMETER_VALUES_FOR_MULTI_GRIDCELL_STATE =
        5  // Incorrect number of grid cells

    MICM_SOLVER_ERROR_CODE_UNUSED_SPECIES = 1,                        // Unused species in system
    MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SYSTEM = 2,               // Chemical system missing
    MICM_SOLVER_ERROR_CODE_MISSING_PROCESSES = 3,                     // No processes defined
    MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SPECIES = 4,              // Species missing from system
    MICM_SOLVER_ERROR_CODE_INVALID_TOLERANCE_SIZE = 5,                // Tolerance size mismatch
    MICM_SOLVER_ERROR_CODE_FAILED_TO_CONVERGE = 6,                    // Solver failed to converge
    MICM_STATE_ERROR_CODE_UNKNOWN_SPECIES = 7,                        // Species not recognized
    MICM_STATE_ERROR_CODE_UNKNOWN_RATE_CONSTANT_PARAMETER = 8,        // Unknown rate constant parameter
    MICM_STATE_ERROR_CODE_INVALID_CONCENTRATION_COUNT = 9,            // Wrong number of concentration values
    MICM_STATE_ERROR_CODE_INVALID_CUSTOM_PARAM_COUNT = 10,            // Wrong number of custom rate params
    MICM_STATE_ERROR_CODE_INVALID_CUSTOM_PARAM_COUNT_MULTI_GRID = 11  // Wrong custom param count for multigrid

    enum class MicmStateErrc {
      UnknownSpecies = 1,                                             // Unknown species
      UnknownRateConstantParameter = 2,                               // Unknown rate constant parameter
      IncorrectNumberOfConcentrationValuesForMultiGridcellState = 3,  // Incorrect number of concentration values
      IncorrectNumberOfCustomRateParameterValues = 4,                 // Incorrect number of custom rate parameter values
      IncorrectNumberOfCustomRateParameterValuesForMultiGridcellState = 5,  // Incorrect number of grid cells
    };

MICM_SOLVER_ERROR_CODE_UNKNOWN_SPECIES = 1,                      // Unknown species
    MICM_SOLVER_ERROR_CODE_UNKNOWN_RATE_CONSTANT_PARAMETER = 2,  // Unknown rate constant parameter
    MICM_SOLVER_ERROR_CODE_INCORRECT_NUMBER_OF_CONCENTRATION_VALUES_FOR_MULTI_GRIDCELL_STATE =
        3,  // Incorrect number of concentration values
    MICM_SOLVER_ERROR_CODE_INCORRECT_NUMBER_OF_CUSTOM_RATE_PARAMETER_VALUES =
        4,  // Incorrect number of custom rate parameter values
    MICM_SOLVER_ERROR_CODE_INCORRECT_NUMBER_OF_CUSTOM_RATE_PARAMETER_VALUES_FOR_MULTI_GRIDCELL_STATE =
        5  // Incorrect number of grid cells

    MICM_SOLVER_ERROR_CODE_FAILED_TO_CONVERGE = 1  // Failed to converge

    namespace std
{
  template<>
  struct is_error_code_enum<MicmSolverErrc> : true_type
  {
  };
}  // namespace std

class MicmProcessErrorCategory : public std::error_category
{
 public:
  const char* name() const noexcept override
  {
    return MICM_ERROR_CATEGORY_PROCESS;
  }

  std::string message(int ev) const override
  {
    switch (static_cast<MicmSolverErrc>(ev))
    {
      case MicmSolverErrc::RateConstantIsNotSet: return "Rate constant is not set";
      case MicmSolverErrc::TransferCoefficientIsNotSet: return "Transfer coefficient is not set";
      case MicmSolverErrc::ReactantDoesNotExist: return "Reactant does not exist";
      case MicmSolverErrc::ProductDoesNotExist: return "Product does not exist";
      // TODO - issue: https://github.com/NCAR/micm/issues/810
      // case MicmSolverErrc::TooManyReactantsForSurfaceReaction: return "A surface reaction can only have one reactant";
      default: return "Unknown error";
    }
  }
};

inline const MicmProcessErrorCategory& MicmProcessError()
{
  static const MicmProcessErrorCategory instance;
  return instance;
}

inline std::error_code make_error_code(MicmSolverErrc e)
{
  return { static_cast<int>(e), MicmProcessError() };
}