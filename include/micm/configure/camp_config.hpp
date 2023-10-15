// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <micm/process/process.hpp>
#include <micm/system/system.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/constants.hpp>

#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>

namespace micm
{
  enum class ConfigParseStatus
  {
    Success,
    None,
    InvalidKey,
    UnknownKey,
    CAMPFilesSectionNotFound,
    CAMPDataSectionNotFound,
    InvalidSpecies,
    InvalidMechanism,
    ObjectTypeNotFound,
    RequiredKeyNotFound
  };

  inline std::string configParseStatusToString(const ConfigParseStatus& status)
  {
    switch (status)
    {
      case ConfigParseStatus::Success: return "Success";
      case ConfigParseStatus::None: return "None";
      case ConfigParseStatus::InvalidKey: return "InvalidKey";
      case ConfigParseStatus::UnknownKey: return "UnknownKey";
      case ConfigParseStatus::CAMPFilesSectionNotFound: return "CAMPFilesSectionNotFound";
      case ConfigParseStatus::CAMPDataSectionNotFound: return "CAMPDataSectionNotFound";
      case ConfigParseStatus::InvalidSpecies: return "InvalidSpecies";
      case ConfigParseStatus::InvalidMechanism: return "InvalidMechanism";
      case ConfigParseStatus::ObjectTypeNotFound: return "ObjectTypeNotFound";
      case ConfigParseStatus::RequiredKeyNotFound: return "RequiredKeyNotFound";
      default: return "Unknown";
    }
  }

  // Solver parameters
  struct SolverParameters
  {
    System system_;
    std::vector<Process> processes_;

    SolverParameters(const System& system, std::vector<Process>&& processes)
        : system_(system),
          processes_(std::move(processes))
    {
    }

    SolverParameters(System&& system, std::vector<Process>&& processes)
        : system_(std::move(system)),
          processes_(std::move(processes))
    {
    }
  };

  class JsonReaderPolicy
  {
    using json = nlohmann::json;

    public:
      std::vector<Species> species_arr_;

      std::vector<UserDefinedRateConstant> user_defined_rate_arr_;
      std::vector<ArrheniusRateConstant> arrhenius_rate_arr_;
      std::vector<TroeRateConstant> troe_rate_arr_;
      std::vector<BranchedRateConstant> branched_rate_arr_;
      std::vector<SurfaceRateConstant> surface_rate_arr_;
      std::vector<TernaryChemicalActivationRateConstant> ternary_rate_arr_;
      std::vector<TunnelingRateConstant> tunneling_rate_arr_;

      // Specific for solver parameters
      Phase gas_phase_;
      std::unordered_map<std::string, Phase> phases_;
      std::vector<Process> processes_;

      // Constants
      // Configure files
      static const inline std::string CAMP_CONFIG = "config.json";
      static const inline std::string SPECIES_CONFIG = "species.json";
      static const inline std::string MECHANISM_CONFIG = "mechanism.json";
      static const inline std::string REACTIONS_CONFIG = "reactions.json";
      static const inline std::string TOLERANCE_CONFIG = "tolerance.json";

      // Common JSON
      static const inline std::string CAMP_DATA = "camp-data";
      static const inline std::string CAMP_FILES = "camp-files";
      static const inline std::string TYPE = "type";

  };
}  // namespace micm
