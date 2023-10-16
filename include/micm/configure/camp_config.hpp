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
      static const inline std::string SPECIES_CONFIG = "species.json";
      static const inline std::string MECHANISM_CONFIG = "mechanism.json";
      static const inline std::string REACTIONS_CONFIG = "reactions.json";
      static const inline std::string TOLERANCE_CONFIG = "tolerance.json";

      // Common JSON
      static const inline std::string CAMP_DATA = "camp-data";
      static const inline std::string CAMP_FILES = "camp-files";
      static const inline std::string TYPE = "type";

      // Functions

      /// @brief Parse configures
      /// @param config_dir Path to a the configuration directory
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
          // return ConfigParseStatus::InvalidCAMPFilePath;
          // TODO: add complete backward compatibility with original parser
          camp_files.push_back(config_dir / SPECIES_CONFIG);
          camp_files.push_back(config_dir / MECHANISM_CONFIG);
        }

        // Merge config JSON from CAMP file list
        json camp_data;
        for (const auto& camp_file : camp_files)
        {
          std::cout << "JsonReaderPolicy.Parse CAMP file" << camp_file << std::endl;
          json config_subset = json::parse(std::ifstream(config_dir / camp_file));

          if (!config_subset.contains(CAMP_DATA))
            return ConfigParseStatus::CAMPDataSectionNotFound;
          // json camp_subset = config_subset[CAMP_DATA];

          std::copy(config_subset.begin(), config_subset.end(),
            std::back_inserter(camp_data));
        }

        ConfigParseStatus status;

        status = Configure(camp_data);

        // Assign the parsed 'Species' to 'Phase'
        gas_phase_ = Phase(species_arr_);

        return status;
      }

    private:

      ConfigParseStatus Configure(const json& config_data)
      {
        ConfigParseStatus status = ConfigParseStatus::None;

        std::vector<json> objects;

        for (const auto& section : config_data)
        {
          for (const auto& object : section)
          {
            objects.push_back(object);
          }
        }

        status = ParseObjectArray(objects);

        return status;
      }

      ConfigParseStatus ParseObjectArray(const std::vector<json>& objects)
      {
        ConfigParseStatus status = ConfigParseStatus::None;

        for (const auto& object : objects)
        {
          if (!ValidateJsonWithKey(object, TYPE))
          {
            status = ConfigParseStatus::ObjectTypeNotFound;
            break;
          }
          std::string type = object[TYPE].get<std::string>();

          std::cout << type << std::endl;
          std::cout << object.dump(4) << std::endl;

          if (type == "CHEM_SPEC")
          {
            status = ParseChemicalSpecies(object);
          }
          else if (type == "RELATIVE_TOLERANCE")
          {
            status = ParseRelativeTolerance(object);
          }
          else if (type == "MECHANISM")
          {
            status = ParseMechanism(object);
          }
          else if (type == "PHOTOLYSIS")
          {
            status = ParsePhotolysis(object);
          }
          else if (type == "EMISSION")
          {
            status = ParseEmission(object);
          }
          else if (type == "FIRST_ORDER_LOSS")
          {
            status = ParseFirstOrderLoss(object);
          }
          else if (type == "ARRHENIUS")
          {
            status = ParseArrhenius(object);
          }
          else if (type == "TROE")
          {
            status = ParseTroe(object);
          }
          else if (type == "BRANCHED" || type == "WENNBERG_NO_RO2")
          {
            status = ParseBranched(object);
          }
          else if (type == "SURFACE")
          {
            status = ParseSurface(object);
          }
          else if (type == "TERNARY_CHEMICAL_ACTIVATION")
          {
            status = ParseTernaryChemicalActivation(object);
          }
          else if (type == "TUNNELING" || type == "WENNBERG_TUNNELING")
          {
            status = ParseTunneling(object);
          }
          else
          {
            status = ConfigParseStatus::UnknownKey;
          }

          if (status != ConfigParseStatus::Success)
            break;
        }

        return status;
      }

      bool ValidateJsonWithKey(const json& object, const std::string& key)
      {
        if (!object.contains(key))
        {
          std::string msg = "Key " + key + " was not found in the config file";
          std::cerr << msg << std::endl;
          return false;
        }
        return true;
      }

      ConfigParseStatus ParseChemicalSpecies(const json& object)
      {
        // required keys
        const std::string NAME = "name";
        std::array<std::string, 1> required_keys = { NAME };

        // Check if it contains the required key(s)
        for (const auto& key : required_keys)
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }
        std::string name = object[NAME].get<std::string>();

        // Load remaining keys as properties
        std::map<std::string, double> properties{};
        for (auto& [key, value] : object.items())
        {
          if (value.is_number_float())
            properties[key] = value;
        }
        species_arr_.push_back(Species(name, properties));

        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseRelativeTolerance(const json& object)
      {
        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseMechanism(const json& object)
      {
        std::vector<std::string> required_keys = { "name", "reactions" };
        for (const auto& key : required_keys)
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }
        std::vector<json> objects;
        for (const auto& element : object["reactions"])
        {
          objects.push_back(element);
        }

        return ParseObjectArray(objects);
      }

      std::vector<Species> ParseReactants(const json& object)
      {
        const std::string QTY = "qty";
        std::vector<Species> reactants;
        for (auto& [key, value] : object.items())
        {
          std::size_t qty = 1;
          if (value.contains(QTY))
            qty = value[QTY];
          for (std::size_t i = 0; i < qty; ++i)
            reactants.push_back(Species(key));
        }
        return reactants;
      }

      std::vector<std::pair<Species, double>> ParseProducts(const json& object)
      {
        const std::string YIELD = "yield";
        constexpr double DEFAULT_YEILD = 1.0;
        std::vector<std::pair<Species, double>> products;
        for (auto& [key, value] : object.items())
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
        return products;
      }

      ConfigParseStatus ParsePhotolysis(const json& object)
      {
        const std::string REACTANTS = "reactants";
        const std::string PRODUCTS = "products";
        const std::string MUSICA_NAME = "MUSICA name";

        for (const auto& key : { REACTANTS, PRODUCTS, MUSICA_NAME })
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }

        auto reactants = ParseReactants(object[REACTANTS]);
        auto products = ParseProducts(object[PRODUCTS]);

        std::string name = "PHOTO." + object[MUSICA_NAME].get<std::string>();

        user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name }));

        std::unique_ptr<UserDefinedRateConstant> rate_ptr =
          std::make_unique<UserDefinedRateConstant>(UserDefinedRateConstantParameters{ .label_ = name });
        processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));

        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseEmission(const json& object)
      {
        const std::string SPECIES = "species";
        const std::string MUSICA_NAME = "MUSICA name";
        for (const auto& key : { SPECIES, MUSICA_NAME })
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }

        std::string species = object["species"].get<std::string>();
        json reactants_object{};
        json products_object{};
        products_object[species] = { { "YIELD", 1.0 } };
        auto reactants = ParseReactants(reactants_object);
        auto products = ParseProducts(products_object);

        std::string name = "EMIS." + object[MUSICA_NAME].get<std::string>();

        user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name }));

        std::unique_ptr<UserDefinedRateConstant> rate_ptr =
          std::make_unique<UserDefinedRateConstant>(UserDefinedRateConstantParameters{ .label_ = name });
        processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));

        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseFirstOrderLoss(const json& object)
      {
        const std::string SPECIES = "species";
        const std::string MUSICA_NAME = "MUSICA name";
        for (const auto& key : { SPECIES, MUSICA_NAME })
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }

        std::string species = object["species"].get<std::string>();
        json reactants_object{};
        json products_object{};
        reactants_object[species] = { {} };
        auto reactants = ParseReactants(reactants_object);
        auto products = ParseProducts(products_object);

        std::string name = "LOSS." + object[MUSICA_NAME].get<std::string>();

        user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name }));

        std::unique_ptr<UserDefinedRateConstant> rate_ptr =
          std::make_unique<UserDefinedRateConstant>(UserDefinedRateConstantParameters{ .label_ = name });
        processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));

        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseArrhenius(const json& object)
      {
        const std::string REACTANTS = "reactants";
        const std::string PRODUCTS = "products";

        // Check required json objects exist
        for (const auto& key : { REACTANTS, PRODUCTS })
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }

        auto reactants = ParseReactants(object[REACTANTS]);
        auto products = ParseProducts(object[PRODUCTS]);

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

        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseTroe(const json& object)
      {
        const std::string REACTANTS = "reactants";
        const std::string PRODUCTS = "products";

        // Check required json objects exist
        for (const auto& key : { REACTANTS, PRODUCTS })
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }

        auto reactants = ParseReactants(object[REACTANTS]);
        auto products = ParseProducts(object[PRODUCTS]);

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

        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseBranched(const json& object)
      {
        const std::string REACTANTS = "reactants";
        const std::string ALKOXY_PRODUCTS = "alkoxy products";
        const std::string NITRATE_PRODUCTS = "nitrate products";
        const std::string X = "X";
        const std::string Y = "Y";
        const std::string A0 = "a0";
        const std::string N = "n";

        // Check required json objects exist
        for (const auto& key : { REACTANTS, ALKOXY_PRODUCTS, NITRATE_PRODUCTS, X, Y, A0, N })
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }

        auto reactants = ParseReactants(object[REACTANTS]);
        auto alkoxy_products = ParseProducts(object[ALKOXY_PRODUCTS]);
        auto nitrate_products = ParseProducts(object[NITRATE_PRODUCTS]);

        BranchedRateConstantParameters parameters;
        parameters.X_ = object[X].get<double>();
        parameters.Y_ = object[Y].get<double>();
        parameters.a0_ = object[A0].get<double>();
        parameters.n_ = object[N].get<int>();

        // Alkoxy branch
        parameters.branch_ = BranchedRateConstantParameters::Branch::Alkoxy;
        branched_rate_arr_.push_back(BranchedRateConstant(parameters));
        std::unique_ptr<BranchedRateConstant> rate_ptr = std::make_unique<BranchedRateConstant>(parameters);
        processes_.push_back(Process(reactants, alkoxy_products, std::move(rate_ptr), gas_phase_));

        // Nitrate branch
        parameters.branch_ = BranchedRateConstantParameters::Branch::Nitrate;
        branched_rate_arr_.push_back(BranchedRateConstant(parameters));
        rate_ptr = std::make_unique<BranchedRateConstant>(parameters);
        processes_.push_back(Process(reactants, nitrate_products, std::move(rate_ptr), gas_phase_));

        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseSurface(const json& object)
      {
        const std::string REACTANTS = "gas-phase reactant";
        const std::string PRODUCTS = "gas-phase products";
        const std::string MUSICA_NAME = "MUSICA name";
        const std::string PROBABILITY = "reaction probability";
        for (const auto& key : { REACTANTS, PRODUCTS, MUSICA_NAME })
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }

        std::string species_name = object[REACTANTS].get<std::string>();
        json reactants_object{};
        reactants_object[species_name] = { {} };
        auto reactants = ParseReactants(reactants_object);
        auto products = ParseProducts(object[PRODUCTS]);

        Species reactant_species = Species("");
        for (auto& species : species_arr_)
        {
          if (species.name_ == species_name)
          {
            reactant_species = species;
            break;
          }
        }
        SurfaceRateConstantParameters parameters{ .label_ = "SURF." + object[MUSICA_NAME].get<std::string>(),
                                                  .species_ = reactant_species };

        if (object.contains(PROBABILITY))
        {
          parameters.reaction_probability_ = object[PROBABILITY].get<double>();
        }

        surface_rate_arr_.push_back(SurfaceRateConstant(parameters));

        std::unique_ptr<SurfaceRateConstant> rate_ptr = std::make_unique<SurfaceRateConstant>(parameters);
        processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));

        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseTernaryChemicalActivation(const json& object)
      {
        const std::string REACTANTS = "reactants";
        const std::string PRODUCTS = "products";

        // Check required json objects exist
        for (const auto& key : { REACTANTS, PRODUCTS })
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }

        auto reactants = ParseReactants(object[REACTANTS]);
        auto products = ParseProducts(object[PRODUCTS]);

        TernaryChemicalActivationRateConstantParameters parameters;
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

        ternary_rate_arr_.push_back(TernaryChemicalActivationRateConstant(parameters));

        std::unique_ptr<TernaryChemicalActivationRateConstant> rate_ptr =
          std::make_unique<TernaryChemicalActivationRateConstant>(parameters);

        processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));

        return ConfigParseStatus::Success;
      }

      ConfigParseStatus ParseTunneling(const json& object)
      {
        const std::string REACTANTS = "reactants";
        const std::string PRODUCTS = "products";

        // Check required json objects exist
        for (const auto& key : { REACTANTS, PRODUCTS })
        {
          if (!ValidateJsonWithKey(object, key))
            return ConfigParseStatus::RequiredKeyNotFound;
        }

        auto reactants = ParseReactants(object[REACTANTS]);
        auto products = ParseProducts(object[PRODUCTS]);

        TunnelingRateConstantParameters parameters;
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

        tunneling_rate_arr_.push_back(TunnelingRateConstant(parameters));

        std::unique_ptr<TunnelingRateConstant> rate_ptr = std::make_unique<TunnelingRateConstant>(parameters);

        processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));

        return ConfigParseStatus::Success;
      }
  };

  /// @brief Public interface to read and parse config
  template<class ConfigTypePolicy = JsonReaderPolicy>
  class SolverConfig : public ConfigTypePolicy
  {
    private:
      ConfigParseStatus last_parse_status_ = ConfigParseStatus::None;

    public:
      /// @brief Reads and parses configures
      /// @param config_dir Path to a the configuration directory
      /// @return an enum indicating the success or failure of the parse
      [[nodiscard]] ConfigParseStatus ReadAndParse(const std::filesystem::path& config_dir)
      {
        last_parse_status_ = this->Parse(config_dir);
        return last_parse_status_;
      }

      /// @brief Creates and returns SolverParameters
      /// @return SolverParameters that contains 'System' and a collection of 'Process'
      SolverParameters GetSolverParams()
      {
        if (last_parse_status_ != ConfigParseStatus::Success)
        {
          std::string msg = "Parsing configuration files failed. The parsing failed with error: " +
                            configParseStatusToString(last_parse_status_);
          throw std::runtime_error(msg);
        }

        return SolverParameters(
          std::move(System(std::move(this->gas_phase_), std::move(this->phases_))), std::move(this->processes_));
      }
  };
}  // namespace micm
