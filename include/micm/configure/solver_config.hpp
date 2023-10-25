// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/system.hpp>
#include <micm/util/constants.hpp>
#include <nlohmann/json.hpp>
#include <variant>

namespace micm
{
  enum class ConfigParseStatus
  {
    Success,
    None,
    InvalidSpeciesFilePath,
    InvalidReactionsFilePath,
    InvalidKey,
    UnknownKey,
    InvalidSpecies,
    CAMPFilesSectionNotFound,
    InvalidCAMPFileCount,
    CAMPDataSectionNotFound,
    InvalidMechanism,
    ObjectTypeNotFound,
    RequiredKeyNotFound,
    ContainsNonStandardKey,
    MutuallyExclusiveOption
  };

  inline std::string configParseStatusToString(const ConfigParseStatus& status)
  {
    switch (status)
    {
      case ConfigParseStatus::Success: return "Success";
      case ConfigParseStatus::None: return "None";
      case ConfigParseStatus::InvalidSpeciesFilePath: return "InvalidSpeciesFilePath";
      case ConfigParseStatus::InvalidReactionsFilePath: return "InvalidReactionsFilePath";
      case ConfigParseStatus::InvalidKey: return "InvalidKey";
      case ConfigParseStatus::UnknownKey: return "UnknownKey";
      case ConfigParseStatus::InvalidSpecies: return "InvalidSpecies";
      case ConfigParseStatus::CAMPFilesSectionNotFound: return "CAMPFilesSectionNotFound";
      case ConfigParseStatus::InvalidCAMPFileCount: return "InvalidCAMPFileCount";
      case ConfigParseStatus::CAMPDataSectionNotFound: return "CAMPDataSectionNotFound";
      case ConfigParseStatus::InvalidMechanism: return "InvalidMechanism";
      case ConfigParseStatus::ObjectTypeNotFound: return "ObjectTypeNotFound";
      case ConfigParseStatus::RequiredKeyNotFound: return "RequiredKeyNotFound";
      case ConfigParseStatus::ContainsNonStandardKey: return "ContainsNonStandardKey";
      case ConfigParseStatus::MutuallyExclusiveOption: return "MutuallyExclusiveOption";
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

  // JSON Configure paser
  class JsonReaderPolicy
  {
    using json = nlohmann::json;

   public:
    // Read from species configure
    std::vector<Species> species_arr_;

    // Read from reaction configure
    std::vector<UserDefinedRateConstant> user_defined_rate_arr_;
    std::vector<ArrheniusRateConstant> arrhenius_rate_arr_;
    std::vector<BranchedRateConstant> branched_rate_arr_;
    std::vector<SurfaceRateConstant> surface_rate_arr_;
    std::vector<TroeRateConstant> troe_rate_arr_;
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
    /// @return True for successful parsing
    ConfigParseStatus Parse(const std::filesystem::path& config_dir)
    {
      // Create configure paths
      std::filesystem::path species_config(config_dir / SPECIES_CONFIG);
      std::filesystem::path mechanism_config(config_dir / MECHANISM_CONFIG);
      std::filesystem::path reactions_config(config_dir / REACTIONS_CONFIG);
      // Note tolerance_config is defined here but not used
      std::filesystem::path tolerance_config(config_dir / TOLERANCE_CONFIG);

      // Look for CAMP config file
      std::filesystem::path camp_config(config_dir / CAMP_CONFIG);
      if (std::filesystem::exists(camp_config))
      {
        json camp_data = json::parse(std::ifstream(camp_config));
        if (!camp_data.contains(CAMP_FILES))
          return ConfigParseStatus::CAMPFilesSectionNotFound;

        std::vector<std::string> camp_files;
        for (const auto& element : camp_data[CAMP_FILES])
        {
          camp_files.push_back(element.get<std::string>());
        }
        if (camp_files.size() != 2)
        {
          std::string err_msg = "CAMP file list should contain two files [species.json, mechanism.json]";
          std::cerr << err_msg << std::endl;
          return ConfigParseStatus::InvalidCAMPFileCount;
        }
        // As a temporary implementation, assume camp files are ordered
        species_config = config_dir / camp_files[0];
        mechanism_config = config_dir / camp_files[1];
      }

      // Current reaction configs should be either mechanism_config or reactions_config
      std::filesystem::path cur_reactions_config;

      // Check if species config exists
      if (!std::filesystem::exists(species_config))
      {
        std::string err_msg = "Species configuration file at path " + species_config.string() + " does not exist\n";
        std::cerr << err_msg << std::endl;
        return ConfigParseStatus::InvalidSpeciesFilePath;
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
        std::cerr << err_msg << std::endl;
        return ConfigParseStatus::InvalidReactionsFilePath;
      }

      // Read species file to create Species and Phase that needs to be known to System and Process

      auto species_status = ConfigureSpecies(species_config);
      if (species_status != ConfigParseStatus::Success)
        return species_status;

      // Assign the parsed 'Species' to 'Phase'
      gas_phase_ = Phase(species_arr_);

      // Read reactions file
      json reaction_data = json::parse(std::ifstream(cur_reactions_config));

      if (!reaction_data.contains(CAMP_DATA))
        return ConfigParseStatus::CAMPDataSectionNotFound;

      std::vector<json> reaction_objects;
      for (const auto& element : reaction_data[CAMP_DATA])
      {
        reaction_objects.push_back(element);
      }

      return ParseObjectArray(reaction_objects);
    }

   private:
    /// @brief Create 'Species' and 'Phase'
    /// @param path to 'Species' file
    /// @return True at success
    ConfigParseStatus ConfigureSpecies(const std::filesystem::path& file)
    {
      ConfigParseStatus status = ConfigParseStatus::None;
      json file_data = json::parse(std::ifstream(file));

      if (!file_data.contains(CAMP_DATA))
        return ConfigParseStatus::CAMPDataSectionNotFound;

      std::vector<json> objects;
      for (const auto& element : file_data[CAMP_DATA])
        objects.push_back(element);

      for (const auto& object : objects)
      {
        if (!ValidateJsonWithKey(object, TYPE))
        {
          status = ConfigParseStatus::ObjectTypeNotFound;
          break;
        }

        std::string type = object[TYPE].get<std::string>();

        if (type == "CHEM_SPEC")
        {
          status = ParseChemicalSpecies(object);
        }
        else if (type == "RELATIVE_TOLERANCE")
        {
          status = ParseRelativeTolerance(object);
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

        if (type == "MECHANISM")
        {
          status = ParseMechanism(object);
        }
        else if (type == "PHOTOLYSIS")
        {
          status = ParsePhotolysis(object);
        }
        else if (type == "ARRHENIUS")
        {
          status = ParseArrhenius(object);
        }
        else if (type == "BRANCHED" || type == "WENNBERG_NO_RO2")
        {
          status = ParseBranched(object);
        }
        else if (type == "TERNARY_CHEMICAL_ACTIVATION")
        {
          status = ParseTernaryChemicalActivation(object);
        }
        else if (type == "TROE")
        {
          status = ParseTroe(object);
        }
        else if (type == "TUNNELING" || type == "WENNBERG_TUNNELING")
        {
          status = ParseTunneling(object);
        }
        else if (type == "EMISSION")
        {
          status = ParseEmission(object);
        }
        else if (type == "FIRST_ORDER_LOSS")
        {
          status = ParseFirstOrderLoss(object);
        }
        else if (type == "SURFACE")
        {
          status = ParseSurface(object);
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

    ConfigParseStatus ParseChemicalSpecies(const json& object)
    {
      // required keys
      const std::string NAME = "name";

      auto status = ValidateSchema(object, { NAME, "type" }, { "tracer type", "absolute tolerance", "diffusion coefficient [m2 s-1]", "molecular weight [kg mol-1]" });
      if (status != ConfigParseStatus::Success)
      {
        return status;
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
      auto status = ValidateSchema(object, { "name", "reactions", "type" }, {});
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      std::vector<json> objects;
      for (const auto& element : object["reactions"])
      {
        objects.push_back(element);
      }

      return ParseObjectArray(objects);
    }

    std::pair<ConfigParseStatus, std::vector<Species>> ParseReactants(const json& object)
    {
      const std::string QTY = "qty";
      std::vector<Species> reactants;

      ConfigParseStatus status = ConfigParseStatus::Success;

      for (auto& [key, value] : object.items())
      {
        std::size_t qty = 1;
        auto new_status = ValidateSchema(value, {}, { "qty" });
        if (new_status != ConfigParseStatus::Success) {
          status = new_status;
        }
        if (value.contains(QTY))
          qty = value[QTY];
        for (std::size_t i = 0; i < qty; ++i)
          reactants.push_back(Species(key));
      }
      return std::make_pair(status, reactants);
    }

    std::pair<ConfigParseStatus, std::vector<std::pair<Species, double>>> ParseProducts(const json& object)
    {
      const std::string YIELD = "yield";
        
      ConfigParseStatus status = ConfigParseStatus::Success;

      constexpr double DEFAULT_YEILD = 1.0;
      std::vector<std::pair<Species, double>> products;
      for (auto& [key, value] : object.items())
      {
        auto new_status = ValidateSchema(value, {}, { "yield" });
        if (new_status != ConfigParseStatus::Success) {
          status = new_status;
        }
        if (value.contains(YIELD))
        {
          products.push_back(std::make_pair(Species(key), value[YIELD]));
        }
        else
        {
          products.push_back(std::make_pair(Species(key), DEFAULT_YEILD));
        }
      }
      return std::make_pair(status, products);
    }

    ConfigParseStatus ParsePhotolysis(const json& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";
      const std::string MUSICA_NAME = "MUSICA name";

      auto status = ValidateSchema(object, { "type", REACTANTS, PRODUCTS, MUSICA_NAME }, {});
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);

      if (reactants.first != ConfigParseStatus::Success)
      {
        return reactants.first;
      }

      if (products.first != ConfigParseStatus::Success)
      {
        return products.first;
      }

      std::string name = "PHOTO." + object[MUSICA_NAME].get<std::string>();

      user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name }));

      std::unique_ptr<UserDefinedRateConstant> rate_ptr =
          std::make_unique<UserDefinedRateConstant>(UserDefinedRateConstantParameters{ .label_ = name });
      processes_.push_back(Process(reactants.second, products.second, std::move(rate_ptr), gas_phase_));

      return ConfigParseStatus::Success;
    }

    ConfigParseStatus ParseArrhenius(const json& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";

      auto status = ValidateSchema(object, { "type", REACTANTS, PRODUCTS }, { "A", "B", "C", "D", "E", "Ea", "MUSICA name" });
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);

      if (reactants.first != ConfigParseStatus::Success)
      {
        return reactants.first;
      }

      if (products.first != ConfigParseStatus::Success)
      {
        return products.first;
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
        if (parameters.C_ != 0)
        {
          std::cerr << "Ea is specified when C is also specified for an Arrhenius reaction. Pick one." << std::endl;
          return ConfigParseStatus::MutuallyExclusiveOption;
        }
        // Calculate 'C' using 'Ea'
        parameters.C_ = -1 * object["Ea"].get<double>() / BOLTZMANN_CONSTANT;
      }

      arrhenius_rate_arr_.push_back(ArrheniusRateConstant(parameters));

      std::unique_ptr<ArrheniusRateConstant> rate_ptr = std::make_unique<ArrheniusRateConstant>(parameters);

      processes_.push_back(Process(reactants.second, products.second, std::move(rate_ptr), gas_phase_));

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

      auto status = ValidateSchema(object, { "type", REACTANTS, ALKOXY_PRODUCTS, NITRATE_PRODUCTS, X, Y, A0, N }, {});
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      auto reactants = ParseReactants(object[REACTANTS]);
      auto alkoxy_products = ParseProducts(object[ALKOXY_PRODUCTS]);
      auto nitrate_products = ParseProducts(object[NITRATE_PRODUCTS]);

      if (reactants.first != ConfigParseStatus::Success)
      {
        return reactants.first;
      }
      if (alkoxy_products.first != ConfigParseStatus::Success)
      {
        return alkoxy_products.first;
      }
      if (nitrate_products.first != ConfigParseStatus::Success)
      {
        return nitrate_products.first;
      }

      BranchedRateConstantParameters parameters;
      parameters.X_ = object[X].get<double>();
      parameters.Y_ = object[Y].get<double>();
      parameters.a0_ = object[A0].get<double>();
      parameters.n_ = object[N].get<int>();

      // Alkoxy branch
      parameters.branch_ = BranchedRateConstantParameters::Branch::Alkoxy;
      branched_rate_arr_.push_back(BranchedRateConstant(parameters));
      std::unique_ptr<BranchedRateConstant> rate_ptr = std::make_unique<BranchedRateConstant>(parameters);
      processes_.push_back(Process(reactants.second, alkoxy_products.second, std::move(rate_ptr), gas_phase_));

      // Nitrate branch
      parameters.branch_ = BranchedRateConstantParameters::Branch::Nitrate;
      branched_rate_arr_.push_back(BranchedRateConstant(parameters));
      rate_ptr = std::make_unique<BranchedRateConstant>(parameters);
      processes_.push_back(Process(reactants.second, nitrate_products.second, std::move(rate_ptr), gas_phase_));

      return ConfigParseStatus::Success;
    }

    ConfigParseStatus ParseTroe(const json& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";

      auto status = ValidateSchema(
          object, { "type", REACTANTS, PRODUCTS }, { "k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N" });
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);
      if (reactants.first != ConfigParseStatus::Success)
      {
        return reactants.first;
      }
      if (products.first != ConfigParseStatus::Success)
      {
        return products.first;
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

      processes_.push_back(Process(reactants.second, products.second, std::move(rate_ptr), gas_phase_));

      return ConfigParseStatus::Success;
    }

    ConfigParseStatus ParseTernaryChemicalActivation(const json& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";

      auto status = ValidateSchema(
          object, { "type", REACTANTS, PRODUCTS }, { "k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N" });
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);
      if (reactants.first != ConfigParseStatus::Success)
      {
        return reactants.first;
      }
      if (products.first != ConfigParseStatus::Success)
      {
        return products.first;
      }

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

      processes_.push_back(Process(reactants.second, products.second, std::move(rate_ptr), gas_phase_));

      return ConfigParseStatus::Success;
    }

    ConfigParseStatus ParseTunneling(const json& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";

      auto status = ValidateSchema(object, { "type", REACTANTS, PRODUCTS }, { "A", "B", "C" });
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);
      if (reactants.first != ConfigParseStatus::Success)
      {
        return reactants.first;
      }
      if (products.first != ConfigParseStatus::Success)
      {
        return products.first;
      }

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

      processes_.push_back(Process(reactants.second, products.second, std::move(rate_ptr), gas_phase_));

      return ConfigParseStatus::Success;
    }

    ConfigParseStatus ParseEmission(const json& object)
    {
      const std::string SPECIES = "species";
      const std::string MUSICA_NAME = "MUSICA name";

      auto status = ValidateSchema(object, { "type", SPECIES, MUSICA_NAME }, {});
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      std::string species = object["species"].get<std::string>();
      json reactants_object{};
      json products_object{};
      products_object[species] = { { "yield", 1.0 } };
      auto reactants = ParseReactants(reactants_object);
      auto products = ParseProducts(products_object);
      if (reactants.first != ConfigParseStatus::Success)
      {
        return reactants.first;
      }
      if (products.first != ConfigParseStatus::Success)
      {
        return products.first;
      }

      std::string name = "EMIS." + object[MUSICA_NAME].get<std::string>();

      user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name }));

      std::unique_ptr<UserDefinedRateConstant> rate_ptr =
          std::make_unique<UserDefinedRateConstant>(UserDefinedRateConstantParameters{ .label_ = name });
      processes_.push_back(Process(reactants.second, products.second, std::move(rate_ptr), gas_phase_));

      return ConfigParseStatus::Success;
    }

    ConfigParseStatus ParseFirstOrderLoss(const json& object)
    {
      const std::string SPECIES = "species";
      const std::string MUSICA_NAME = "MUSICA name";

      auto status = ValidateSchema(object, { "type", SPECIES, MUSICA_NAME }, {});
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      std::string species = object["species"].get<std::string>();
      json reactants_object{};
      json products_object{};
      reactants_object[species] = { };
      auto reactants = ParseReactants(reactants_object);
      auto products = ParseProducts(products_object);
      if (reactants.first != ConfigParseStatus::Success)
      {
        return reactants.first;
      }
      if (products.first != ConfigParseStatus::Success)
      {
        return products.first;
      }

      std::string name = "LOSS." + object[MUSICA_NAME].get<std::string>();

      user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name }));

      std::unique_ptr<UserDefinedRateConstant> rate_ptr =
          std::make_unique<UserDefinedRateConstant>(UserDefinedRateConstantParameters{ .label_ = name });
      processes_.push_back(Process(reactants.second, products.second, std::move(rate_ptr), gas_phase_));

      return ConfigParseStatus::Success;
    }

    ConfigParseStatus ParseSurface(const json& object)
    {
      const std::string REACTANTS = "gas-phase reactant";
      const std::string PRODUCTS = "gas-phase products";
      const std::string MUSICA_NAME = "MUSICA name";
      const std::string PROBABILITY = "reaction probability";

      auto status = ValidateSchema(object, { "type", REACTANTS, PRODUCTS, MUSICA_NAME }, { PROBABILITY });
      if (status != ConfigParseStatus::Success)
      {
        return status;
      }

      std::string species_name = object[REACTANTS].get<std::string>();
      json reactants_object{};
      reactants_object[species_name] = { };
      auto reactants = ParseReactants(reactants_object);
      auto products = ParseProducts(object[PRODUCTS]);
      if (reactants.first != ConfigParseStatus::Success)
      {
        return reactants.first;
      }
      if (products.first != ConfigParseStatus::Success)
      {
        return products.first;
      }

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
      processes_.push_back(Process(reactants.second, products.second, std::move(rate_ptr), gas_phase_));

      return ConfigParseStatus::Success;
    }

    /// @brief Search for nonstandard keys. Only nonstandard keys starting with __ are allowed. Others are considered typos
    /// @param object the object whose keys need to be validated
    /// @param required_keys The required keys
    /// @param optional_keys The optional keys
    /// @return true if only standard keys are found
    ConfigParseStatus ValidateSchema(
        const json& object,
        const std::vector<std::string>& required_keys,
        const std::vector<std::string>& optional_keys)
    {
      // standard keys are:
      // those in required keys
      // those in optional keys
      // starting with __
      // anything else is reported as an error so that typos are caught, specifically for optional keys

      std::vector<std::string> sorted_object_keys;
      for (auto& [key, value] : object.items())
        sorted_object_keys.push_back(key);

      auto sorted_required_keys = required_keys;
      auto sorted_optional_keys = optional_keys;
      std::sort(sorted_object_keys.begin(), sorted_object_keys.end());
      std::sort(sorted_required_keys.begin(), sorted_required_keys.end());
      std::sort(sorted_optional_keys.begin(), sorted_optional_keys.end());

      // get the difference between the object keys and those required
      // what's left should be the optional keys and valid comments
      std::vector<std::string> difference;
      std::set_difference(
          sorted_object_keys.begin(),
          sorted_object_keys.end(),
          sorted_required_keys.begin(),
          sorted_required_keys.end(),
          std::back_inserter(difference));

      // check that the number of keys remaining is exactly equal to the expected number of required keys
      if (difference.size() != (sorted_object_keys.size() - required_keys.size()))
      {
        return ConfigParseStatus::RequiredKeyNotFound;
      }

      std::vector<std::string> remaining;
      std::set_difference(
          difference.begin(),
          difference.end(),
          sorted_optional_keys.begin(),
          sorted_optional_keys.end(),
          std::back_inserter(remaining));

      // now, anything left must be standard comment starting with __
      for (auto& key : remaining)
      {
        if (!key.starts_with("__"))
        {
          return ConfigParseStatus::ContainsNonStandardKey;
        }
      }
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
    /// @param config_dir A path to a configuration file
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
