// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/system.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/error.hpp>

#include <yaml-cpp/yaml.h>

#include <array>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <system_error>

enum class MicmConfigErrc
{
  InvalidKey = MICM_CONFIGURATION_ERROR_CODE_INVALID_KEY,
  UnknownKey = MICM_CONFIGURATION_ERROR_CODE_UNKNOWN_KEY,
  InvalidFilePath = MICM_CONFIGURATION_ERROR_CODE_INVALID_FILE_PATH,
  NoConfigFilesFound = MICM_CONFIGURATION_ERROR_CODE_NO_CONFIG_FILES_FOUND,
  CAMPFilesNotFound = MICM_CONFIGURATION_ERROR_CODE_CAMP_FILES_NOT_FOUND,
  CAMPDataNotFound = MICM_CONFIGURATION_ERROR_CODE_CAMP_DATA_NOT_FOUND,
  InvalidSpecies = MICM_CONFIGURATION_ERROR_CODE_INVALID_SPECIES,
  InvalidMechanism = MICM_CONFIGURATION_ERROR_CODE_INVALID_MECHANISM,
  InvalidType = MICM_CONFIGURATION_ERROR_CODE_INVALID_TYPE,
  ObjectTypeNotFound = MICM_CONFIGURATION_ERROR_CODE_OBJECT_TYPE_NOT_FOUND,
  RequiredKeyNotFound = MICM_CONFIGURATION_ERROR_CODE_REQUIRED_KEY_NOT_FOUND,
  ContainsNonStandardKey = MICM_CONFIGURATION_ERROR_CODE_CONTAINS_NON_STANDARD_KEY,
  MutuallyExclusiveOption = MICM_CONFIGURATION_ERROR_CODE_MUTUALLY_EXCLUSIVE_OPTION
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmConfigErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class MicmConfigErrorCategory : public std::error_category
  {
   public:
    const char* name() const noexcept override
    {
      return MICM_ERROR_CATEGORY_CONFIGURATION;
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmConfigErrc>(ev))
      {
        case MicmConfigErrc::InvalidKey: return "Invalid key";
        case MicmConfigErrc::UnknownKey: return "Unknown key";
        case MicmConfigErrc::InvalidFilePath: return "Invalid file path";
        case MicmConfigErrc::NoConfigFilesFound: return "No config files found";
        case MicmConfigErrc::CAMPFilesNotFound: return "CAMP files not found";
        case MicmConfigErrc::CAMPDataNotFound: return "CAMP data not found";
        case MicmConfigErrc::InvalidSpecies: return "Invalid species";
        case MicmConfigErrc::InvalidMechanism: return "Invalid mechanism";
        case MicmConfigErrc::InvalidType: return "Invalid type";
        case MicmConfigErrc::ObjectTypeNotFound: return "Object type not found";
        case MicmConfigErrc::RequiredKeyNotFound: return "Required key not found";
        case MicmConfigErrc::ContainsNonStandardKey: return "Contains non-standard key";
        case MicmConfigErrc::MutuallyExclusiveOption: return "Mutually exclusive option";
        default: return "Unknown error";
      }
    }
  };

  const MicmConfigErrorCategory MICM_CONFIG_ERROR{};
}  // namespace

inline std::error_code make_error_code(MicmConfigErrc e)
{
  return { static_cast<int>(e), MICM_CONFIG_ERROR };
}

namespace micm
{

  constexpr double MOLES_M3_TO_MOLECULES_CM3 = 1.0e-6 * constants::AVOGADRO_CONSTANT;

  // Solver parameters
  struct SolverParameters
  {
    System system_;
    std::vector<Process> processes_;
    RosenbrockSolverParameters parameters_;
    double relative_tolerance_;

    SolverParameters(const System& system, std::vector<Process>&& processes, const RosenbrockSolverParameters&& parameters)
        : system_(system),
          processes_(processes),
          parameters_(parameters),
          relative_tolerance_(0)
    {
    }

    SolverParameters(System&& system, std::vector<Process>&& processes, RosenbrockSolverParameters&& parameters)
        : system_(system),
          processes_(processes),
          parameters_(parameters),
          relative_tolerance_(0)
    {
    }
    SolverParameters(
        System&& system,
        std::vector<Process>&& processes,
        RosenbrockSolverParameters&& parameters,
        double relative_tolerance)
        : system_(system),
          processes_(processes),
          parameters_(parameters),
          relative_tolerance_(relative_tolerance)
    {
    }
  };

  class ConfigReaderPolicy
  {
    using objectType = YAML::Node;

   public:
    std::map<std::string, Species> species_;

    std::vector<UserDefinedRateConstant> user_defined_rate_arr_;
    std::vector<ArrheniusRateConstant> arrhenius_rate_arr_;
    std::vector<TroeRateConstant> troe_rate_arr_;
    std::vector<TernaryChemicalActivationRateConstant> ternary_rate_arr_;
    std::vector<BranchedRateConstant> branched_rate_arr_;
    std::vector<TunnelingRateConstant> tunneling_rate_arr_;
    std::vector<SurfaceRateConstant> surface_rate_arr_;

    // Specific for solver parameters
    Phase gas_phase_;
    std::unordered_map<std::string, Phase> phases_;
    std::vector<Process> processes_;
    RosenbrockSolverParameters parameters_;
    double relative_tolerance_;

    // Common YAML
    inline static const std::string DEFAULT_CONFIG_FILE_JSON = "config.json";
    inline static const std::string DEFAULT_CONFIG_FILE_YAML = "config.yaml";
    inline static const std::string CAMP_FILES = "camp-files";
    inline static const std::string CAMP_DATA = "camp-data";
    inline static const std::string TYPE = "type";

    // Error string
    std::stringstream last_json_object_;

    // Constructor

    ConfigReaderPolicy(const RosenbrockSolverParameters& parameters)
        : parameters_(parameters)
    {
    }

    // Functions

    /// @brief Parse configures
    /// @param config_path Path to a the CAMP configuration directory or file
    void Parse(const std::filesystem::path& config_path)
    {
      // Look for CAMP config path
      if (!std::filesystem::exists(config_path))
      {
        throw std::system_error{ make_error_code(MicmConfigErrc::InvalidFilePath), config_path.string() };
      }

      std::filesystem::path config_dir;
      std::filesystem::path config_file;

      if (std::filesystem::is_directory(config_path))
      {
        // If config path is a directory, use default config file name
        config_dir = config_path;
        if (std::filesystem::exists(config_dir / DEFAULT_CONFIG_FILE_YAML))
        {
          config_file = config_dir / DEFAULT_CONFIG_FILE_YAML;
        }
        else
        {
          config_file = config_dir / DEFAULT_CONFIG_FILE_JSON;
        }
      }
      else
      {
        // Extract configuration dir from configuration file path
        config_dir = config_path.parent_path();
        config_file = config_path;
      }

      // Load the CAMP file list YAML
      objectType camp_data = YAML::LoadFile(config_file.string());
      if (!camp_data[CAMP_FILES])
      {
        throw std::system_error{ make_error_code(MicmConfigErrc::CAMPFilesNotFound), config_file.string() };
      }

      // Build a list of individual CAMP config files
      std::vector<std::filesystem::path> camp_files;
      for (const auto& element : camp_data[CAMP_FILES])
      {
        std::filesystem::path camp_file = config_dir / element.as<std::string>();
        if (!std::filesystem::exists(camp_file))
        {
          throw std::system_error{ make_error_code(MicmConfigErrc::CAMPFilesNotFound), camp_file.string() };
        }
        camp_files.push_back(camp_file);
      }

      // No config files found
      if (camp_files.size() < 1)
      {
        throw std::system_error{ make_error_code(MicmConfigErrc::NoConfigFilesFound), config_file.string() };
      }

      std::vector<objectType> species_objects;
      std::vector<objectType> mechanism_objects;

      // Iterate CAMP file list and form CAMP data object arrays
      for (const auto& camp_file : camp_files)
      {
        objectType config_subset = YAML::LoadFile(camp_file.string());

        if (!config_subset[CAMP_DATA])
        {
          YAML::Emitter out;
          out << config_subset;
          throw std::system_error{ make_error_code(MicmConfigErrc::CAMPDataNotFound), out.c_str() };
        }
        // Iterate YAML objects from CAMP data entry
        for (const auto& object : config_subset[CAMP_DATA])
        {
          if (object)
          {
            // Require object to have a type entry
            if (!object[TYPE])
            {
              YAML::Emitter out;
              out << object;
              std::string msg = "object: " + std::string(out.c_str()) + "; type: " + TYPE;

              throw std::system_error{ make_error_code(MicmConfigErrc::ObjectTypeNotFound), msg };
            }
            // Sort into object arrays by type
            std::string type = object[TYPE].as<std::string>();
            // CHEM_SPEC and RELATIVE_TOLERANCE parsed first by ParseSpeciesArray
            if ((type == "CHEM_SPEC") || (type == "RELATIVE_TOLERANCE"))
            {
              species_objects.push_back(object);
            }
            // All other objects will be parsed by ParseMechanismArray
            else
            {
              mechanism_objects.push_back(object);
            }
          }
        }
      }

      // Clear vectors and maps
      species_.clear();
      user_defined_rate_arr_.clear();
      arrhenius_rate_arr_.clear();
      troe_rate_arr_.clear();
      ternary_rate_arr_.clear();
      branched_rate_arr_.clear();
      tunneling_rate_arr_.clear();
      surface_rate_arr_.clear();
      phases_.clear();
      processes_.clear();

      // Parse species object array
      ParseSpeciesArray(species_objects);

      // Assign the parsed 'Species' to 'Phase'
      std::vector<Species> species_arr;
      for (const auto& [name, species] : species_)
      {
        species_arr.push_back(species);
      }
      gas_phase_ = Phase(species_arr);

      // Parse mechanism object array
      ParseMechanismArray(mechanism_objects);
    }

   private:
    void ParseSpeciesArray(const std::vector<objectType>& objects)
    {
      for (const auto& object : objects)
      {
        std::string type = object[TYPE].as<std::string>();

        if (type == "CHEM_SPEC")
        {
          ParseChemicalSpecies(object);
        }
        else if (type == "RELATIVE_TOLERANCE")
        {
          ParseRelativeTolerance(object);
        }
      }
    }

    void ParseMechanismArray(const std::vector<objectType>& objects)
    {
      for (const auto& object : objects)
      {
        std::string type = object[TYPE].as<std::string>();

        if (type == "MECHANISM")
        {
          ParseMechanism(object);
        }
        else if (type == "PHOTOLYSIS")
        {
          ParsePhotolysis(object);
        }
        else if (type == "EMISSION")
        {
          ParseEmission(object);
        }
        else if (type == "FIRST_ORDER_LOSS")
        {
          ParseFirstOrderLoss(object);
        }
        else if (type == "ARRHENIUS")
        {
          ParseArrhenius(object);
        }
        else if (type == "TROE")
        {
          ParseTroe(object);
        }
        else if (type == "TERNARY_CHEMICAL_ACTIVATION")
        {
          ParseTernaryChemicalActivation(object);
        }
        else if (type == "BRANCHED" || type == "WENNBERG_NO_RO2")
        {
          ParseBranched(object);
        }
        else if (type == "TUNNELING" || type == "WENNBERG_TUNNELING")
        {
          ParseTunneling(object);
        }
        else if (type == "SURFACE")
        {
          ParseSurface(object);
        }
        else if (type == "USER_DEFINED")
        {
          ParseUserDefined(object);
        }
        else
        {
          throw std::system_error{ make_error_code(MicmConfigErrc::UnknownKey), type };
        }
      }
    }

    void ParseChemicalSpecies(const objectType& object)
    {
      // required keys
      const std::string NAME = "name";
      const std::string TYPE = "type";

      const std::string TRACER_TYPE = "tracer type";
      const std::string ABS_TOLERANCE = "absolute tolerance";
      const std::string DIFFUSION_COEFF = "diffusion coefficient [m2 s-1]";
      const std::string MOL_WEIGHT = "molecular weight [kg mol-1]";
      const std::string THIRD_BODY = "THIRD_BODY";

      ValidateSchema(object, { NAME, TYPE }, { TRACER_TYPE, ABS_TOLERANCE, DIFFUSION_COEFF, MOL_WEIGHT });

      std::string name = object[NAME].as<std::string>();
      Species species{ name };

      // Load remaining keys as properties
      for (auto it = object.begin(); it != object.end(); ++it)
      {
        auto key = it->first.as<std::string>();
        auto value = it->second;

        if (key.empty())
        {
          throw std::system_error{ make_error_code(MicmConfigErrc::InvalidType), key };
        }

        if (key != NAME && key != TYPE)
        {
          std::string stringValue = value.as<std::string>();

          if (stringValue.empty())
          {
            species.SetProperty<std::string>(key, stringValue);
          }
          else if (IsInt(stringValue))
          {
            species.SetProperty<int>(key, value.as<int>());
          }
          else if (IsFloat(stringValue))
          {
            species.SetProperty<double>(key, value.as<double>());
          }
          else if (IsBool(stringValue))
          {
            species.SetProperty<bool>(key, value.as<bool>());
          }
          else if (key == TRACER_TYPE && stringValue == THIRD_BODY)
          {
            species.SetThirdBody();
          }
          else
          {
            species.SetProperty<std::string>(key, stringValue);
          }
        }
      }
      species_[name] = species;
    }

    void ParseRelativeTolerance(const objectType& object)
    {
      ValidateSchema(object, { "value", "type" }, {});
      this->relative_tolerance_ = object["value"].as<double>();
    }

    void ParseMechanism(const objectType& object)
    {
      ValidateSchema(object, { "name", "reactions", "type" }, {});
      std::vector<objectType> objects;
      for (const auto& element : object["reactions"])
      {
        objects.push_back(element);
      }
      ParseMechanismArray(objects);
    }

    std::vector<Species> ParseReactants(const objectType& object)
    {
      const std::string QTY = "qty";
      std::vector<Species> reactants;

      for (auto it = object.begin(); it != object.end(); ++it)
      {
        auto key = it->first.as<std::string>();
        auto value = it->second;

        std::size_t qty = 1;
        ValidateSchema(value, {}, { "qty" });
        if (value[QTY])
          qty = value[QTY].as<std::size_t>();
        for (std::size_t i = 0; i < qty; ++i)
          reactants.push_back(species_[key]);
      }
      return reactants;
    }

    std::vector<std::pair<Species, double>> ParseProducts(const objectType& object)
    {
      const std::string YIELD = "yield";

      constexpr double DEFAULT_YIELD = 1.0;
      std::vector<std::pair<Species, double>> products;
      for (auto it = object.begin(); it != object.end(); ++it)
      {
        auto key = it->first.as<std::string>();
        auto value = it->second;

        ValidateSchema(value, {}, { "yield" });
        auto species = species_[key];
        if (value[YIELD])
        {
          double yield = value[YIELD].as<double>();
          products.push_back(std::make_pair(species, yield));
        }
        else
        {
          products.push_back(std::make_pair(species, DEFAULT_YIELD));
        }
      }
      return products;
    }

    void ParsePhotolysis(const objectType& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";
      const std::string MUSICA_NAME = "MUSICA name";
      const std::string SCALING_FACTOR = "scaling factor";

      ValidateSchema(object, { "type", REACTANTS, PRODUCTS, MUSICA_NAME }, { SCALING_FACTOR });

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);
      double scaling_factor = object[SCALING_FACTOR] ? object[SCALING_FACTOR].as<double>() : 1.0;

      std::string name = "PHOTO." + object[MUSICA_NAME].as<std::string>();

      user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name, .scaling_factor_ = scaling_factor }));

      std::unique_ptr<UserDefinedRateConstant> rate_ptr = std::make_unique<UserDefinedRateConstant>(
          UserDefinedRateConstantParameters{ .label_ = name, .scaling_factor_ = scaling_factor });
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));
    }

    void ParseArrhenius(const objectType& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";

      ValidateSchema(object, { "type", REACTANTS, PRODUCTS }, { "A", "B", "C", "D", "E", "Ea", "MUSICA name" });

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);

      ArrheniusRateConstantParameters parameters;
      if (object["A"])
      {
        parameters.A_ = object["A"].as<double>();
      }
      int n_reactants = reactants.size();
      parameters.A_ *= std::pow(MOLES_M3_TO_MOLECULES_CM3, n_reactants - 1);
      if (object["B"])
      {
        parameters.B_ = object["B"].as<double>();
      }
      if (object["C"])
      {
        parameters.C_ = object["C"].as<double>();
      }
      if (object["D"])
      {
        parameters.D_ = object["D"].as<double>();
      }
      if (object["E"])
      {
        parameters.E_ = object["E"].as<double>();
      }
      if (object["Ea"])
      {
        if (parameters.C_ != 0)
        {
          throw std::system_error{ make_error_code(MicmConfigErrc::MutuallyExclusiveOption),
                                   "Ea is specified when C is also specified for an Arrhenius reaction. Pick one." };
        }
        // Calculate 'C' using 'Ea'
        parameters.C_ = -1 * object["Ea"].as<double>() / constants::BOLTZMANN_CONSTANT;
      }
      arrhenius_rate_arr_.push_back(ArrheniusRateConstant(parameters));
      std::unique_ptr<ArrheniusRateConstant> rate_ptr = std::make_unique<ArrheniusRateConstant>(parameters);
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));
    }

    void ParseTroe(const objectType& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";

      ValidateSchema(
          object, { "type", REACTANTS, PRODUCTS }, { "k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N" });

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);

      TroeRateConstantParameters parameters;
      if (object["k0_A"])
      {
        parameters.k0_A_ = object["k0_A"].as<double>();
      }
      // Account for the conversion of reactant concentrations (including M) to molecules cm-3
      int n_reactants = reactants.size();
      parameters.k0_A_ *= std::pow(MOLES_M3_TO_MOLECULES_CM3, n_reactants);
      if (object["k0_B"])
      {
        parameters.k0_B_ = object["k0_B"].as<double>();
      }
      if (object["k0_C"])
      {
        parameters.k0_C_ = object["k0_C"].as<double>();
      }
      if (object["kinf_A"])
      {
        parameters.kinf_A_ = object["kinf_A"].as<double>();
      }
      // Account for terms in denominator and exponent that include [M] but not other reactants
      parameters.kinf_A_ *= std::pow(MOLES_M3_TO_MOLECULES_CM3, n_reactants - 1);
      if (object["kinf_B"])
      {
        parameters.kinf_B_ = object["kinf_B"].as<double>();
      }
      if (object["kinf_C"])
      {
        parameters.kinf_C_ = object["kinf_C"].as<double>();
      }
      if (object["Fc"])
      {
        parameters.Fc_ = object["Fc"].as<double>();
      }
      if (object["N"])
      {
        parameters.N_ = object["N"].as<double>();
      }
      troe_rate_arr_.push_back(TroeRateConstant(parameters));
      std::unique_ptr<TroeRateConstant> rate_ptr = std::make_unique<TroeRateConstant>(parameters);
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));
    }

    void ParseTernaryChemicalActivation(const objectType& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";

      ValidateSchema(
          object, { "type", REACTANTS, PRODUCTS }, { "k0_A", "k0_B", "k0_C", "kinf_A", "kinf_B", "kinf_C", "Fc", "N" });

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);

      TernaryChemicalActivationRateConstantParameters parameters;
      if (object["k0_A"])
      {
        parameters.k0_A_ = object["k0_A"].as<double>();
      }
      // Account for the conversion of reactant concentrations (including M) to molecules cm-3
      int n_reactants = reactants.size();
      parameters.k0_A_ *= std::pow(MOLES_M3_TO_MOLECULES_CM3, n_reactants - 1);
      if (object["k0_B"])
      {
        parameters.k0_B_ = object["k0_B"].as<double>();
      }
      if (object["k0_C"])
      {
        parameters.k0_C_ = object["k0_C"].as<double>();
      }
      if (object["kinf_A"])
      {
        parameters.kinf_A_ = object["kinf_A"].as<double>();
      }
      // Account for terms in denominator and exponent that include [M] but not other reactants
      parameters.kinf_A_ *= std::pow(MOLES_M3_TO_MOLECULES_CM3, n_reactants - 2);
      if (object["kinf_B"])
      {
        parameters.kinf_B_ = object["kinf_B"].as<double>();
      }
      if (object["kinf_C"])
      {
        parameters.kinf_C_ = object["kinf_C"].as<double>();
      }
      if (object["Fc"])
      {
        parameters.Fc_ = object["Fc"].as<double>();
      }
      if (object["N"])
      {
        parameters.N_ = object["N"].as<double>();
      }

      ternary_rate_arr_.push_back(TernaryChemicalActivationRateConstant(parameters));
      std::unique_ptr<TernaryChemicalActivationRateConstant> rate_ptr =
          std::make_unique<TernaryChemicalActivationRateConstant>(parameters);
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));
    }

    void ParseBranched(const objectType& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string ALKOXY_PRODUCTS = "alkoxy products";
      const std::string NITRATE_PRODUCTS = "nitrate products";
      const std::string X = "X";
      const std::string Y = "Y";
      const std::string A0 = "a0";
      const std::string N = "n";

      ValidateSchema(object, { "type", REACTANTS, ALKOXY_PRODUCTS, NITRATE_PRODUCTS, X, Y, A0, N }, {});

      auto reactants = ParseReactants(object[REACTANTS]);
      auto alkoxy_products = ParseProducts(object[ALKOXY_PRODUCTS]);
      auto nitrate_products = ParseProducts(object[NITRATE_PRODUCTS]);

      BranchedRateConstantParameters parameters;
      parameters.X_ = object[X].as<double>();
      // Account for the conversion of reactant concentrations to molecules cm-3
      int n_reactants = reactants.size();
      parameters.X_ *= std::pow(MOLES_M3_TO_MOLECULES_CM3, n_reactants - 1);
      parameters.Y_ = object[Y].as<double>();
      parameters.a0_ = object[A0].as<double>();
      parameters.n_ = object[N].as<int>();

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
    }

    void ParseTunneling(const objectType& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";

      ValidateSchema(object, { "type", REACTANTS, PRODUCTS }, { "A", "B", "C" });

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);

      TunnelingRateConstantParameters parameters;
      if (object["A"])
      {
        parameters.A_ = object["A"].as<double>();
      }
      // Account for the conversion of reactant concentrations to molecules cm-3
      int n_reactants = reactants.size();
      parameters.A_ *= std::pow(MOLES_M3_TO_MOLECULES_CM3, n_reactants - 1);
      if (object["B"])
      {
        parameters.B_ = object["B"].as<double>();
      }
      if (object["C"])
      {
        parameters.C_ = object["C"].as<double>();
      }

      tunneling_rate_arr_.push_back(TunnelingRateConstant(parameters));
      std::unique_ptr<TunnelingRateConstant> rate_ptr = std::make_unique<TunnelingRateConstant>(parameters);
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));
    }

    void ParseEmission(const objectType& object)
    {
      const std::string SPECIES = "species";
      const std::string MUSICA_NAME = "MUSICA name";
      const std::string PRODUCTS = "products";
      const std::string SCALING_FACTOR = "scaling factor";

      ValidateSchema(object, { "type", SPECIES, MUSICA_NAME }, { SCALING_FACTOR, PRODUCTS });

      std::string species = object["species"].as<std::string>();
      objectType reactants_object{};
      objectType products_object{};
      products_object[species]["yield"] = 1.0;
      auto reactants = ParseReactants(reactants_object);
      auto products = ParseProducts(products_object);
      double scaling_factor = object[SCALING_FACTOR] ? object[SCALING_FACTOR].as<double>() : 1.0;

      std::string name = "EMIS." + object[MUSICA_NAME].as<std::string>();
      user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name, .scaling_factor_ = scaling_factor }));
      std::unique_ptr<UserDefinedRateConstant> rate_ptr = std::make_unique<UserDefinedRateConstant>(
          UserDefinedRateConstantParameters{ .label_ = name, .scaling_factor_ = scaling_factor });
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));
    }

    void ParseFirstOrderLoss(const objectType& object)
    {
      const std::string SPECIES = "species";
      const std::string MUSICA_NAME = "MUSICA name";
      const std::string SCALING_FACTOR = "scaling factor";

      ValidateSchema(object, { "type", SPECIES, MUSICA_NAME }, { SCALING_FACTOR });

      std::string species = object["species"].as<std::string>();
      objectType reactants_object{};
      objectType products_object{};
      reactants_object[species] = {};
      auto reactants = ParseReactants(reactants_object);
      auto products = ParseProducts(products_object);
      double scaling_factor = object[SCALING_FACTOR] ? object[SCALING_FACTOR].as<double>() : 1.0;

      std::string name = "LOSS." + object[MUSICA_NAME].as<std::string>();
      user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name, .scaling_factor_ = scaling_factor }));
      std::unique_ptr<UserDefinedRateConstant> rate_ptr = std::make_unique<UserDefinedRateConstant>(
          UserDefinedRateConstantParameters{ .label_ = name, .scaling_factor_ = scaling_factor });
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));
    }

    void ParseUserDefined(const objectType& object)
    {
      const std::string REACTANTS = "reactants";
      const std::string PRODUCTS = "products";
      const std::string MUSICA_NAME = "MUSICA name";
      const std::string SCALING_FACTOR = "scaling factor";

      ValidateSchema(object, { "type", REACTANTS, PRODUCTS, MUSICA_NAME }, { SCALING_FACTOR });

      auto reactants = ParseReactants(object[REACTANTS]);
      auto products = ParseProducts(object[PRODUCTS]);
      double scaling_factor = object[SCALING_FACTOR] ? object[SCALING_FACTOR].as<double>() : 1.0;

      std::string name = "USER." + object[MUSICA_NAME].as<std::string>();
      user_defined_rate_arr_.push_back(UserDefinedRateConstant({ .label_ = name, .scaling_factor_ = scaling_factor }));
      std::unique_ptr<UserDefinedRateConstant> rate_ptr = std::make_unique<UserDefinedRateConstant>(
          UserDefinedRateConstantParameters{ .label_ = name, .scaling_factor_ = scaling_factor });
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));
    }

    void ParseSurface(const objectType& object)
    {
      const std::string REACTANTS = "gas-phase reactant";
      const std::string PRODUCTS = "gas-phase products";
      const std::string MUSICA_NAME = "MUSICA name";
      const std::string PROBABILITY = "reaction probability";

      ValidateSchema(object, { "type", REACTANTS, PRODUCTS, MUSICA_NAME }, { PROBABILITY });

      std::string species_name = object[REACTANTS].as<std::string>();
      objectType reactants_object{};
      reactants_object[species_name] = {};

      auto reactants = ParseReactants(reactants_object);
      auto products = ParseProducts(object[PRODUCTS]);

      Species reactant_species = Species("");
      reactant_species = species_[species_name];
      SurfaceRateConstantParameters parameters{ .label_ = "SURF." + object[MUSICA_NAME].as<std::string>(),
                                                .species_ = reactant_species };

      if (object[PROBABILITY])
      {
        parameters.reaction_probability_ = object[PROBABILITY].as<double>();
      }

      surface_rate_arr_.push_back(SurfaceRateConstant(parameters));
      std::unique_ptr<SurfaceRateConstant> rate_ptr = std::make_unique<SurfaceRateConstant>(parameters);
      processes_.push_back(Process(reactants, products, std::move(rate_ptr), gas_phase_));
    }

    // Utility functions to check types and perform conversions
    bool IsBool(const std::string& value)
    {
      return (value == "true" || value == "false");
    }

    bool IsInt(const std::string& value)
    {
      std::istringstream iss(value);
      int result;
      return (iss >> result >> std::ws).eof();  // Check if the entire string is an integer
    }

    bool IsFloat(const std::string& value)
    {
      std::istringstream iss(value);
      float result;
      return (iss >> result >> std::ws).eof();  // Check if the entire string is a float
    }

    /// @brief Search for nonstandard keys. Only nonstandard keys starting with __ are allowed. Others are considered typos
    /// @param object the object whose keys need to be validated
    /// @param required_keys The required keys
    /// @param optional_keys The optional keys
    void ValidateSchema(
        const objectType& object,
        const std::vector<std::string>& required_keys,
        const std::vector<std::string>& optional_keys)
    {
      YAML::Emitter out;
      out << object;  // Serialize the object to a string

      // standard keys are:
      // those in required keys
      // those in optional keys
      // starting with __
      // anything else is reported as an error so that typos are caught, specifically for optional keys
      if (object && object.size() > 0 && object.begin()->second.IsNull())
      {
        return;
      }

      std::vector<std::string> sorted_object_keys;
      for (auto it = object.begin(); it != object.end(); ++it)
      {
        std::string key = it->first.as<std::string>();
        sorted_object_keys.push_back(key);
      }

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
        std::vector<std::string> missing_keys;
        std::set_difference(
            sorted_required_keys.begin(),
            sorted_required_keys.end(),
            sorted_object_keys.begin(),
            sorted_object_keys.end(),
            std::back_inserter(missing_keys));
        std::string msg;
        for (auto& key : missing_keys)
        {
          msg += "Missing required key '" + key + "' in object: " + out.c_str() + "\n";
        }

        throw std::system_error{ make_error_code(MicmConfigErrc::RequiredKeyNotFound), msg };
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
          std::string msg = "Non-standard key '" + key + "' found in object" + out.c_str();
          throw std::system_error{ make_error_code(MicmConfigErrc::ContainsNonStandardKey), msg };
        }
      }
    }
  };

  /// @brief Public interface to read and parse config
  template<class ConfigTypePolicy = ConfigReaderPolicy>
  class SolverConfig : public ConfigTypePolicy
  {
   public:
    SolverConfig()
        : ConfigTypePolicy(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
    {
    }
    SolverConfig(const RosenbrockSolverParameters& parameters)
        : ConfigTypePolicy(parameters)
    {
    }

    /// @brief Reads and parses configures
    /// @param config_dir Path to a the configuration directory
    void ReadAndParse(const std::filesystem::path& config_dir)
    {
      this->Parse(config_dir);
    }

    /// @brief Creates and returns SolverParameters
    /// @return SolverParameters that contains 'System' and a collection of 'Process'
    SolverParameters GetSolverParams()
    {
      return SolverParameters(
          std::move(System(this->gas_phase_, this->phases_)),
          std::move(this->processes_),
          std::move(this->parameters_),
          this->relative_tolerance_);
    }
  };
}  // namespace micm
