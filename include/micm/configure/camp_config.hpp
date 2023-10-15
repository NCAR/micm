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
    InvalidCAMPFilePath,
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
      case ConfigParseStatus::InvalidCAMPFilePath: return "InvalidCAMPFilePath";
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

      // Common JSON
      static const inline std::string CAMP_DATA = "camp-data";
      static const inline std::string CAMP_FILES = "camp-files";
      static const inline std::string TYPE = "type";

      // Functions

      /// @brief Parse configures
      /// @return True for successful parsing
      ConfigParseStatus Parse(const std::filesystem::path& config_dir)
      {
        std::vector<std::string> camp_files;

        // Look for CAMP config file
        std::filesystem::path camp_config(config_dir / CAMP_CONFIG);
        if (std::filesystem::exists(camp_config))
        {
          json camp_data = json::parse(std::ifstream(camp_config));
          if (!camp_data.contains(CAMP_FILES))
            return ConfigParseStatus::CAMPFilesSectionNotFound;

          for (const auto& element : camp_data[CAMP_FILES])
          {
            camp_files.push_back(element.get<std::string>());
          }
        }
        else
        {
          return ConfigParseStatus::InvalidCAMPFilePath;
        }

        // Merge config JSON from CAMP file list
        for (const auto& camp_file : camp_files)
        {
          json config_data = json::parse(std::ifstream(camp_file));
        }

        return ConfigParseStatus::Success;
      }
  };
}  // namespace micm
