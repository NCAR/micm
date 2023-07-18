/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 */

#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/property.hpp>
#include <micm/system/species.hpp>
#include <micm/system/system.hpp>
#include <micm/util/constants.hpp>
#include <nlohmann/json.hpp>
#include <variant>

namespace micm
{
  // JSON-Parsing error policies
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

  // JSON Configure paser
  template<class ErrorPolicy>
  class JsonReaderPolicy : public ErrorPolicy
  {
    using json = nlohmann::json;

   public:
    // Read from species configure
    std::vector<Species> species_arr_;

    // Read from reaction configure
    std::vector<PhotolysisRateConstant> photolysis_rate_arr_;
    std::vector<ArrheniusRateConstant> arrhenius_rate_arr_;
    std::vector<TroeRateConstant> troe_rate_arr_;
    std::vector<Species> emission_arr_;
    std::vector<Species> first_order_loss_arr_;

    // Specific for solver parameters
    Phase gas_phase_;
    std::unordered_map<std::string, Phase> phases_;
    std::vector<Process> processes_;

    // Status of parsing
    bool is_parse_success_ = false;

    // Constants
    // Configure files 
    static const inline std::string SPECIES_CONFIG = "species.json";
    static const inline std::string MECHANISM_CONFIG = "mechanism.json";
    static const inline std::string REACTIONS_CONFIG = "reactions.json";
    static const inline std::string TOLERANCE_CONFIG = "tolerance.json";

    // Common JSON 
    static const inline std::string CAMP_DATA = "camp-data";
    static const inline std::string TYPE = "type";

    // Functions

    /// @brief Parse configures
    /// @return True for successful parsing
    bool Parse(const std::filesystem::path& config_dir)
    {
      // Create configure paths
      std::filesystem::path species_config(config_dir / SPECIES_CONFIG);
      std::filesystem::path mechanism_config(config_dir / MECHANISM_CONFIG);
      std::filesystem::path reactions_config(config_dir / REACTIONS_CONFIG);
      std::filesystem::path tolerance_config(
          config_dir / TOLERANCE_CONFIG);  // TODO: jiwon 7/13 - we don't need to worry about this now

      // Current reaction configs should be either mechanism_config or reactions config
      std::filesystem::path cur_reactions_config;

      // Check if species config exists
      if (!std::filesystem::exists(species_config))
      {
        std::string err_msg = "Species configuration file at path " + species_config.string() + " does not exist\n";
        this->OnError(err_msg);

        return false;
      }

      // Check if a reaction configure exists and decide which one
      if (std::filesystem::exists(mechanism_config))
      {
        cur_reactions_config = mechanism_config;
      }
      else if (std::filesystem::exists(reactions_config))
      {
        cur_reactions_config = reactions_config;
      }
      else
      {
        std::string err_msg = "Reaction configuration file at path " + mechanism_config.string() + " or " +
                              reactions_config.string() + " does not exist\n";
        this->OnError(err_msg);

        return false;
      }

      // Read species file to create Species and Phase that needs to be known to System and Process
      if (!ConfigureSpecies(species_config))
        return false;

      // Assign the parsed 'Species' to 'Phase'
      gas_phase_ = Phase(species_arr_);

      // Read reactions file
      json reaction_data = json::parse(std::ifstream(cur_reactions_config));

      if (!ValidateJsonWithKey(reaction_data, CAMP_DATA))
        return false;

      std::vector<json> reaction_objects;
      for (const auto& element : reaction_data[CAMP_DATA])
      {
        reaction_objects.push_back(element);
      }

      if (!ParseObjectArray(reaction_objects))
        return false;

      is_parse_success_ = true;

      return is_parse_success_;
    }

   private:
    /// @brief Create 'Species' and 'Phase'
    /// @param path to 'Species' file
    /// @return True at success
    bool ConfigureSpecies(const std::filesystem::path& file)
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
        else if (type == "TROE")
        {
          if (!ParseTroe(object))
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
      // required keys
      const std::string NAME = "name";

      // optional keys
      const std::string ABS_TOL = "absolute tolerance";
      const std::string MOL_WEIGHT = "molecular weight [kg mol-1]";
      const std::string MOL_WEIGHT_UNIT = "kg mol-1";

      std::array<std::string, 1> required_keys = { NAME };

      // Check if it contains the required key(s)
      for (const auto& key : required_keys)
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }

      // Check if it contains optional key(s)
      std::string name = object[NAME].get<std::string>();

      if (object.contains(ABS_TOL))
      {
        auto species = Species(name, Property(ABS_TOL, "", object[ABS_TOL].get<double>()));
        species_arr_.push_back(species);
      }
      else if (object.contains(MOL_WEIGHT))
      {
        auto species = Species(name, Property(MOL_WEIGHT, MOL_WEIGHT_UNIT, object[MOL_WEIGHT].get<double>()));
        species_arr_.push_back(species);
      }
      else
      {
        species_arr_.push_back(Species(name));
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

      constexpr double DEFAULT_YEILD = 1.0;

      for (const auto& key : { REACTANTS, PRODUCTS, MUSICA_NAME })
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }

      std::vector<Species> reactants;
      for (auto& [key, value] : object[REACTANTS].items())
      {
        reactants.push_back(Species(key));
      }

      std::vector<std::pair<Species, double>> products;
      for (auto& [key, value] : object[PRODUCTS].items())
      {
        if (value.contains(YIELD))
        {
          products.push_back(std::make_pair(Species(key), value[YIELD]));
        }
        else
        {
          products.push_back(std::make_pair(Species(key), DEFAULT_YEILD));
        }
      }

      std::string name = object[MUSICA_NAME].get<std::string>();

      photolysis_rate_arr_.push_back(PhotolysisRateConstant(name));

      std::unique_ptr<PhotolysisRateConstant> rate_ptr = std::make_unique<PhotolysisRateConstant>(name);
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));

      return true;
    }

    bool ParseArrhenius(const json& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";
      const std::string YIELD = "yield";

      constexpr double DEFAULT_YEILD = 1.0;

      // Check required json objects exist
      for (const auto& key : { REACTANTS, PRODUCTS })
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }

      // Create process
      std::vector<Species> reactants;
      for (auto& [key, value] : object[REACTANTS].items())
      {
        reactants.push_back(Species(key));
      }

      std::vector<std::pair<Species, double>> products;
      for (auto& [key, value] : object[PRODUCTS].items())
      {
        if (value.contains(YIELD))
        {
          products.push_back(std::make_pair(Species(key), value[YIELD]));
        }
        else
        {
          products.push_back(std::make_pair(Species(key), DEFAULT_YEILD));
        }
      }

      ArrheniusRateConstantParameters parameters;
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
      if (object.contains("Ea"))
      {
        // Calculate 'C' using 'Ea'
        parameters.C_ = -1 * object["Ea"].get<double>() / BOLTZMANN_CONSTANT;
      }

      arrhenius_rate_arr_.push_back(ArrheniusRateConstant(parameters));

      std::unique_ptr<ArrheniusRateConstant> rate_ptr = std::make_unique<ArrheniusRateConstant>(parameters);

      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));

      return true;
    }

    bool ParseTroe(const json& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";
      const std::string YIELD = "yield";

      constexpr double DEFAULT_YEILD = 1.0;

      // Check required json objects exist
      for (const auto& key : { REACTANTS, PRODUCTS })
      {
        if (!ValidateJsonWithKey(object, key))
          return false;
      }

      // Create process
      std::vector<Species> reactants;
      for (auto& [key, value] : object[REACTANTS].items())
      {
        reactants.push_back(Species(key));
      }

      std::vector<std::pair<Species, double>> products;
      for (auto& [key, value] : object[PRODUCTS].items())
      {
        if (value.contains(YIELD))
        {
          products.push_back(std::make_pair(Species(key), value[YIELD]));
        }
        else
        {
          products.push_back(std::make_pair(Species(key), DEFAULT_YEILD));
        }
      }

      TroeRateConstantParameters parameters;
      if (object.contains("k0_A"))
      {
        parameters.k0_A_ = object["k0_A"].get<double>();
      }
      if (object.contains("k0_B"))
      {
        parameters.k0_B_ = object["k0_B"].get<double>();
      }
      if (object.contains("k0_C"))
      {
        parameters.k0_C_ = object["k0_C"].get<double>();
      }
      if (object.contains("kinf_A"))
      {
        parameters.kinf_A_ = object["kinf_A"].get<double>();
      }
      if (object.contains("kinf_B"))
      {
        parameters.kinf_B_ = object["kinf_B"].get<double>();
      }
      if (object.contains("kinf_C"))
      {
        parameters.kinf_C_ = object["kinf_C"].get<double>();
      }
      if (object.contains("Fc"))
      {
        parameters.Fc_ = object["Fc"].get<double>();
      }
      if (object.contains("N"))
      {
        parameters.N_ = object["N"].get<double>();
      }

      troe_rate_arr_.push_back(TroeRateConstant(parameters));

      std::unique_ptr<TroeRateConstant> rate_ptr = std::make_unique<TroeRateConstant>(parameters);

      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));

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

      emission_arr_.push_back(Species(name));

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

      first_order_loss_arr_.push_back(Species(name));

      return true;
    }
  };

  /// @brief Public interface to read and parse config
  template<template<class> class ConfigTypePolicy = JsonReaderPolicy, template<class> class ErrorPolicy = ThrowPolicy>
  class SolverConfig : public ConfigTypePolicy<ErrorPolicy<std::exception>>  // TODO jiwon 7/3 : what should be the template
                                                                             // speciailization for this?
  {
   public:
    /// @brief Reads and parses configures
    /// @return TRUE if parsing is success
    bool ReadAndParse(const std::filesystem::path& config_dir)
    {
      return this->Parse(config_dir);
    }

    /// @brief Creates and returns SolverParameters
    /// @return SolverParameters that contains 'System' and a collection of 'Process'
    SolverParameters GetSolverParams()
    {
      if (!this->is_parse_success_)
      {
        throw std::runtime_error("Parsing configure files hasn't been completed");
      }

      return SolverParameters(
          std::move(System(std::move(this->gas_phase_), std::move(this->phases_))), std::move(this->processes_));
    }

    /// @brief Get a collection of 'PhotolysisRateConstant'
    /// @return a collection of 'PhotolysisRateConstant'
    std::vector<PhotolysisRateConstant>& GetPhotolysisRateConstants()
    {
      if (!this->is_parse_success_)
      {
        throw std::runtime_error("Parsing configure files hasn't been completed");
      }

      return this->photolysis_rate_arr_;
    }
  };

}  // namespace micm