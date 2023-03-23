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
#include <micm/process/intraphase_process.hpp>
#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/property.hpp>
#include <micm/system/species.hpp>
#include <micm/system/system.hpp>

#ifdef USE_JSON
#  include <nlohmann/json.hpp>
#endif

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

  template<class ErrorPolicy>
  class JsonReaderPolicy : public ErrorPolicy
  {
    using json = nlohmann::json;
    std::vector<Species> species_;
    std::vector<Species> emissions_;
    std::vector<Species> first_order_loss_;
    std::vector<IntraPhaseProcess<PhotolysisRateConstant>> photolysis_reactions_;
    std::vector<IntraPhaseProcess<ArrheniusRateConstant>> arrhenius_reactions_;

   public:
    std::unique_ptr<micm::System> ReadAndParse(std::filesystem::path path)
    {
      species_.clear();
      emissions_.clear();
      first_order_loss_.clear();
      photolysis_reactions_.clear();
      arrhenius_reactions_.clear();

      if (!std::filesystem::exists(path))
      {
        std::string err_msg = "Configuration file at path " + path.string() + " does not exist\n";
        return this->OnError(err_msg);
      }

      json data = json::parse(std::ifstream(path));

      std::string key = "camp-files";
      ValidateJsonWithKey(data, key);

      std::vector<std::string> files = data[key].get<std::vector<nlohmann::json::string_t>>();

      key = "camp-data";
      for (const auto& file : files)
      {
        if (!std::filesystem::exists(path))
        {
          std::string err_msg = "Configuration file at path " + path.string() + " does not exist\n";
          return this->OnError(err_msg);
        }
        json file_data = json::parse(std::ifstream(file));
        ValidateJsonWithKey(file_data, key);
        std::vector<json> objects;
        for (const auto& element : file_data[key])
        {
          objects.push_back(element);
        }

        ParseObjectArray(objects);
      }

      return std::make_unique<micm::System>(micm::System(micm::Phase(species_)));
    }

    void ValidateJsonWithKey(const json& object, std::string key)
    {
      if (!object.contains(key))
      {
        this->OnError("Key " + key + " was not found in the config file");
      }
    }

    void ParseObjectArray(const std::vector<json>& objects)
    {
      for (const auto& object : objects)
      {
        std::string key = "type";
        ValidateJsonWithKey(object, key);
        std::string type = object[key].get<std::string>();

        if (type == "CHEM_SPEC")
        {
          ParseChemicalSpecies(object);
        }
        else if (type == "RELATIVE_TOLERANCE")
        {
          ParseRelativeTolerance(object);
        }
        else if (type == "MECHANISM")
        {
          ParseMechanism(object);
        }
        else if (type == "PHOTOLYSIS")
        {
          ParsePhotolysis(object);
        }
        else if (type == "ARRHENIUS")
        {
          ParseArrhenius(object);
        }
        else if (type == "EMISSION")
        {
          ParseEmission(object);
        }
        else if (type == "FIRST_ORDER_LOSS")
        {
          ParseFirstOrderLoss(object);
        }
        else
        {
          this->OnError("Unknown key in config file: " + type);
        }
      }
    }

    void ParseChemicalSpecies(const json& object)
    {
      std::vector<std::string> required_keys = { "name" };
      std::vector<std::string> optional_keys = { "absolute tolerance" };
      for (const auto& key : required_keys)
        ValidateJsonWithKey(object, key);

      std::string name = object["name"].get<std::string>();

      if (object.contains("absolute tolerance"))
      {
        double abs_tol = object["absolute tolerance"].get<double>();
        species_.push_back(Species(name, Property("absolute tolerance", "", abs_tol)));
      }
      else
      {
        species_.push_back(Species(name));
      }
    }

    void ParseRelativeTolerance(const json& object)
    {
      // TODO: what is this?
    }

    void ParseMechanism(const json& object)
    {
      std::vector<std::string> required_keys = { "name", "reactions" };
      for (const auto& key : required_keys)
        ValidateJsonWithKey(object, key);

      std::vector<json> objects;
      for (const auto& element : object["reactions"])
      {
        objects.push_back(element);
      }

      ParseObjectArray(objects);
    }

    void ParsePhotolysis(const json& object)
    {
      std::vector<std::string> required_keys = { "reactants", "products", "MUSICA name" };
      for (const auto& key : required_keys)
        ValidateJsonWithKey(object, key);

      std::vector<Species> reactants;
      for (auto& [key, value] : object["reactants"].items()) {
          reactants.push_back(Species(key));
      }
      std::vector<Species> products;
      for (auto& [key, value] : object["products"].items()) {
          products.push_back(Species(key));
      }
      std::string name = object["MUSICA name"].get<std::string>();

      using reaction = IntraPhaseProcess<PhotolysisRateConstant>;
      photolysis_reactions_.push_back(
        reaction(reactants, products, PhotolysisRateConstant(0, name))
      );
    }

    void ParseArrhenius(const json& object)
    {
      std::vector<std::string> required_keys = { "reactants", "products" };
      for (const auto& key : required_keys)
        ValidateJsonWithKey(object, key);

      std::vector<Species> reactants;
      for (auto& [key, value] : object["reactants"].items()) {
          reactants.push_back(Species(key));
      }
      std::vector<Species> products;
      for (auto& [key, value] : object["products"].items()) {
          products.push_back(Species(key));
      }

      ArrheniusRateConstantParameters parameters;

      if (object.contains("A")){
        parameters.A_ = object["A"].get<double>();
      }
      if (object.contains("B")){
        parameters.B_ = object["B"].get<double>();
      }
      if (object.contains("C")){
        parameters.C_ = object["C"].get<double>();
      }
      if (object.contains("D")){
        parameters.D_ = object["D"].get<double>();
      }
      if (object.contains("E")){
        parameters.E_ = object["E"].get<double>();
      }

      using reaction = IntraPhaseProcess<ArrheniusRateConstant>;
      arrhenius_reactions_.push_back(
        reaction(reactants, products, ArrheniusRateConstant(parameters))
      );

    }

    void ParseEmission(const json& object)
    {
      std::vector<std::string> required_keys = { "species" };
      std::vector<std::string> optional_keys = { "MUSICA name" };
      for (const auto& key : required_keys)
        ValidateJsonWithKey(object, key);

      std::string name = object["species"].get<std::string>();

      emissions_.push_back(Species(name));
    }

    void ParseFirstOrderLoss(const json& object)
    {
      std::vector<std::string> required_keys = { "species" };
      std::vector<std::string> optional_keys = { "MUSICA name" };
      for (const auto& key : required_keys)
        ValidateJsonWithKey(object, key);

      std::string name = object["species"].get<std::string>();

      first_order_loss_.push_back(Species(name));
    }
  };

  template<
    template<class> class ConfigTypePolicy = JsonReaderPolicy, 
    template<class> class ErrorPolicy = ThrowPolicy
  >
  class SystemBuilder : public ConfigTypePolicy<ErrorPolicy<std::unique_ptr<micm::System>>>
  {
   public:
    std::unique_ptr<micm::System> Build(std::filesystem::path path)
    {
      return this->ReadAndParse(path);
    }
  };

}  // namespace micm