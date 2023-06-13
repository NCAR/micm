/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 */

#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/property.hpp>
#include <micm/system/species.hpp>
#include <micm/system/system.hpp>
#include <variant>

#ifdef USE_JSON
#  include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace micm
{
  template<class Object>
  class ThrowPolicy
  {
    class Exception : public std::exception
    {
     public:
      const char* msg_;

     public:
      Exception(const char* msg)
          : msg_(msg)
      {
      }

      virtual const char* what()
      {
        return msg_;
      }
    };

   public:
    Object OnError(std::string message)
    {
      throw Exception(message.c_str());
    }
  };

  template<class Object>
  class NoThrowPolicy
  {
   public:
    Object OnError(std::string message)
    {
      std::cerr << message << std::endl;
      return Object();
    }
  };

  // Solover parameters
  struct SolverParameters
  {
    micm::System system_;
    std::vector<micm::Process> processes_;
  };

  // Error code
  enum class ConfigErrorCode
  {
    None = 0,
    NoAcess,
    FileNotFound,
    KeyNotFound,
  };

  // JSON Configure paser
  template<class ErrorPolicy>
  class JsonReaderPolicy : public ErrorPolicy
  {
   public:
    // Set to "public" for "NoThrowPolcy"
    std::vector<Species> species_;
    std::vector<Species> emissions_;
    std::vector<Species> first_order_loss_;
    std::vector<micm::Process> processes_;
    micm::Phase gas_phase_;
    std::unordered_map<std::string, micm::Phase> phases_;

    // Constants
    static const inline std::string SPECIES_CONFIG =
        "species.json";  // TODO:jiwon 6/6 - instead of searching, pass the configure path
    static const inline std::string REACTIONS_CONFIG = "mechanism.json";  // TODO:jiwon 6/6

    static const inline std::string CAMP_FILES = "camp-files";
    static const inline std::string CAMP_DATA = "camp-data";

    static const inline std::string TYPE = "type";

    // Functions

    /// @brief read and parse JSON objects
    /// @param
    /// @return SolverParameters if parsing is success, else returns ConfigErrorCode
    std::variant<micm::SolverParameters, micm::ConfigErrorCode> ReadAndParse(const std::filesystem::path& path)
    {
      // Check whether file exists
      if (!std::filesystem::exists(path))
      {
        std::string err_msg = "Configuration file at path " + path.string() + " does not exist\n";
        this->OnError(err_msg);

        return micm::ConfigErrorCode::FileNotFound;
      }

      // Read file to get the list of configure files
      json data = json::parse(std::ifstream(path));
      if (!ValidateJsonWithKey(data, CAMP_FILES))
      {
        return micm::ConfigErrorCode::KeyNotFound;
      }

      // Check whether the listed files exist and determine the sequence to read files.
      std::string species_file;
      std::vector<std::string> other_files;
      bool found_species_file = false;

      for (const auto& file : data[CAMP_FILES].get<std::vector<nlohmann::json::string_t>>())
      {
        if (!std::filesystem::exists(file))
        {
          std::string err_msg = "Configuration file at path " + file + " does not exist\n";
          this->OnError(err_msg);

          return micm::ConfigErrorCode::FileNotFound;
        }

        // Find species file to read first
        std::size_t found = file.find(SPECIES_CONFIG);
        if (found != std::string::npos)
        {
          species_file = file;
          found_species_file = true;
        }
        else
        {
          other_files.push_back(file);
        }
      }

      // Read species file to create Species and Phase that needs to be known to System and Process
      if (found_species_file)
      {
        if (!ConfigureSpecies(species_file))
        {
          return micm::ConfigErrorCode::KeyNotFound;
        }
      }
      else
      {
        std::string err_msg = "Species configure file does not exist\n";
        this->OnError(err_msg);

        return micm::ConfigErrorCode::FileNotFound;
      }

      // Read files, eg. reactions.json
      for (const auto& file : other_files)
      {
        json file_data = json::parse(std::ifstream(file));
        if (!ValidateJsonWithKey(file_data, CAMP_DATA))
        {
          return micm::ConfigErrorCode::KeyNotFound;
        }

        std::vector<json> objects;
        for (const auto& element : file_data[CAMP_DATA])
        {
          objects.push_back(element);
        }

        if (!ParseObjectArray(objects))
        {
          return micm::ConfigErrorCode::KeyNotFound;
        }
      }

      micm::SystemParameters sysParams = { gas_phase_, phases_ };

      return micm::SolverParameters{ micm::System(sysParams), processes_ };
    }

    /// @brief Create 'Species' and 'Phase'
    /// @param path to 'Species' file
    /// @return True at success
    bool ConfigureSpecies(const std::string& file)
    {
      json file_data = json::parse(std::ifstream(file));

      if (!ValidateJsonWithKey(file_data, CAMP_DATA))
        return false;

      std::vector<json> objects;
      for (const auto& element : file_data[CAMP_DATA])
        objects.push_back(element);

      for (const auto& object : objects)
      {
        if (!ValidateJsonWithKey(object, TYPE))
          return false;

        std::string type = object[TYPE].get<std::string>();

        if (type == "CHEM_SPEC")
        {
          if (!ParseChemicalSpecies(object))
            return false;
        }
        else if (type == "RELATIVE_TOLERANCE")
        {
          if (!ParseRelativeTolerance(object))
            return false;
        }
      }
      // After creating Species, create Phase
      gas_phase_.species_ = species_;

      return true;
    }

    bool ValidateJsonWithKey(const json& object, const std::string& key)
    {
      if (!object.contains(key))
      {
        this->OnError("Key " + key + " was not found in the config file");

        return false;
      }
      return true;
    }

    bool ParseObjectArray(const std::vector<json>& objects)
    {
      for (const auto& object : objects)
      {
        if (!ValidateJsonWithKey(object, TYPE))
          return false;

        std::string type = object[TYPE].get<std::string>();

        if (type == "MECHANISM")
        {
          if (!ParseMechanism(object))
            return false;
        }
        else if (type == "PHOTOLYSIS")
        {
          if (!ParsePhotolysis(object))
            return false;
        }
        else if (type == "ARRHENIUS")
        {
          if (!ParseArrhenius(object))
            return false;
        }
        else if (type == "EMISSION")
        {
          if (!ParseEmission(object))
            return false;
        }
        else if (type == "FIRST_ORDER_LOSS")
        {
          if (!ParseFirstOrderLoss(object))
            return false;
        }
        else
        {
          this->OnError("Unknown key in config file: " + type);
          return false;
        }
      }
      return true;
    }

    bool ParseChemicalSpecies(const json& object)
    {
      std::vector<std::string> required_keys = { "name" };
      std::vector<std::string> optional_keys = { "absolute tolerance" };

      for (const auto& key : required_keys)
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }

      std::string name = object["name"].get<std::string>();

      std::string key = "absolute tolerance";

      if (object.contains(key))
      {
        double abs_tol = object[key].get<double>();
        auto species = Species(name, Property(key, "", abs_tol));
        species_.push_back(species);
      }
      else
      {
        species_.push_back(Species(name));
      }

      return true;
    }

    bool ParseRelativeTolerance(const json& object)
    {
      // TODO: what is this?
      return true;
    }

    bool ParseMechanism(const json& object)
    {
      std::vector<std::string> required_keys = { "name", "reactions" };
      for (const auto& key : required_keys)
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }
      std::vector<json> objects;
      for (const auto& element : object["reactions"])
      {
        objects.push_back(element);
      }

      return ParseObjectArray(objects);
    }

    bool ParsePhotolysis(const json& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";
      const std::string MUSICA_NAME = "MUSICA name";
      const std::string YIELD = "yield";

      const double DEFAULT_YEILD = 1.0;

      for (const auto& key : { REACTANTS, PRODUCTS, MUSICA_NAME })
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }

      // Create process
      std::vector<micm::Species> reactants;
      for (auto& [key, value] : object[REACTANTS].items())
      {
        reactants.push_back(micm::Species(key));
      }

      std::vector<std::pair<micm::Species, double>> products;
      for (auto& [key, value] : object[PRODUCTS].items())
      {
        if (value.contains(YIELD))
        {
          products.push_back(std::make_pair(micm::Species(key), value[YIELD]));
        }
        else
        {
          products.push_back(std::make_pair(micm::Species(key), DEFAULT_YEILD));
        }
      }

      std::unique_ptr<micm::PhotolysisRateConstant> rate_ptr =
          std::make_unique<micm::PhotolysisRateConstant>(object[MUSICA_NAME].get<std::string>());

      processes_.push_back(micm::Process(reactants, products, std::move(rate_ptr), gas_phase_));

      return true;
    }

    bool ParseArrhenius(const json& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";
      const std::string YIELD = "yield";

      const double DEFAULT_YEILD = 1.0;

      // Check required json objects exist
      for (const auto& key : { REACTANTS, PRODUCTS })
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }

      // Create process
      std::vector<micm::Species> reactants;
      for (auto& [key, value] : object[REACTANTS].items())
      {
        reactants.push_back(micm::Species(key));
      }

      std::vector<std::pair<micm::Species, double>> products;
      for (auto& [key, value] : object[PRODUCTS].items())
      {
        if (value.contains(YIELD))
        {
          products.push_back(std::make_pair(micm::Species(key), value[YIELD]));
        }
        else
        {
          products.push_back(std::make_pair(micm::Species(key), DEFAULT_YEILD));
        }
      }

      micm::ArrheniusRateConstantParameters parameters;
      if (object.contains("A"))
      {
        parameters.A_ = object["A"].get<double>();
      }
      if (object.contains("B"))
      {
        parameters.B_ = object["B"].get<double>();
      }
      if (object.contains("C"))
      {
        parameters.C_ = object["C"].get<double>();
      }
      if (object.contains("D"))
      {
        parameters.D_ = object["D"].get<double>();
      }
      if (object.contains("E"))
      {
        parameters.E_ = object["E"].get<double>();
      }

      std::unique_ptr<micm::ArrheniusRateConstant> rate_ptr =
          std::make_unique<micm::ArrheniusRateConstant>(micm::ArrheniusRateConstantParameters(parameters));

      processes_.push_back(micm::Process(reactants, products, std::move(rate_ptr), gas_phase_));

      return true;
    }

    bool ParseEmission(const json& object)
    {
      std::vector<std::string> required_keys = { "species" };
      std::vector<std::string> optional_keys = { "MUSICA name" };
      for (const auto& key : required_keys)
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }

      std::string name = object["species"].get<std::string>();

      emissions_.push_back(Species(name));

      return true;
    }

    bool ParseFirstOrderLoss(const json& object)
    {
      std::vector<std::string> required_keys = { "species" };
      std::vector<std::string> optional_keys = { "MUSICA name" };
      for (const auto& key : required_keys)
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }

      std::string name = object["species"].get<std::string>();

      first_order_loss_.push_back(Species(name));

      return true;
    }
  };

  template<template<class> class ConfigTypePolicy = JsonReaderPolicy, template<class> class ErrorPolicy = NoThrowPolicy>
  class SolverConfig : public ConfigTypePolicy<ErrorPolicy<std::variant<micm::SolverParameters, micm::ConfigErrorCode>>>
  {
   public:
    std::variant<micm::SolverParameters, micm::ConfigErrorCode> Configure(const std::filesystem::path& path)
    {
      return this->ReadAndParse(path);
    }
  };

}  // namespace micm
#endif