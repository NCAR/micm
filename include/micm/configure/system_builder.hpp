/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 */

#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <micm/system/system.hpp>
#include <micm/system/species.hpp>
#include <micm/system/property.hpp>

#ifdef USE_JSON
#include <nlohmann/json.hpp>
#endif

namespace micm
{

  class ThrowPolicy {
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

    void OnError(std::string message) {
      throw Exception(message.c_str());
    }
  };

  template<class Object>
  class NoThrowPolicy {
    public:

      void OnError(std::string message) {
        std::cerr << message << std::endl;
      }
  };

  template<
    class ErrorPolicy
  >
  class JsonReaderPolicy : public ErrorPolicy
  {
    using json = nlohmann::json;
    std::vector<std::unique_ptr<Species>> species_;

   public:
    std::unique_ptr<micm::System> ReadAndParse(std::filesystem::path path)
    {
      species_.clear();

      ValidateFilePath(path);

      json data = json::parse(std::ifstream(path));

      std::string key = "camp-files";
      ValidateJsonWithKey(data, key);

      std::vector<std::string> files = data[key].get<std::vector<nlohmann::json::string_t>>();

      key = "camp-data";
      for(const auto& file : files){
        ValidateFilePath(file);
        json file_data = json::parse(std::ifstream(file));
        std::vector<json> objects;
        for(const auto& element : file_data[key]){
          objects.push_back(element);
        }

        ParseObjectArray(objects);
      }

      return std::make_unique<micm::System>(micm::System());
    }

    void ValidateFilePath(std::filesystem::path path) {
      if (!std::filesystem::exists(path))
      {
        std::string err_msg = "Configuration file at path " + path.string() + " does not exist\n";
        this->OnError(err_msg);
      }
    }

    void ValidateJsonWithKey(const json& object, std::string key) {
      if (!object.contains(key)) {
        this->OnError("Key " + key + " was not found in the config file");
      }
    }

    void ParseObjectArray(const std::vector<json>& objects){
      for(const auto& object : objects){
        std::string key = "type";
        ValidateJsonWithKey(object, key);
        std::string type = object[key].get<std::string>();
        if (type == "CHEM_SPEC"){
          ParseChemicalSpecies(object);
        }
        else if (type == "RELATIVE_TOLERANCE") {
          ParseRelativeTolerance(object);
        }
        else if (type == "MECHANISM") {
          ParseMechanism(object);
        }
        else {
          this->OnError("Unknown key in config file: " + type);
        }
      }
    }

    void ParseChemicalSpecies(const json& object){
      std::vector<std::string> required_keys = {"name"};
      std::vector<std::string> optional_keys = {"absolute tolerance"};

      for(const auto& key : required_keys) ValidateJsonWithKey(object, key);

      std::string name = object["name"].get<std::string>();

      std::unique_ptr<Species> species;

      if (object.contains("absolute tolerance")){
        double abs_tol = object["absolute tolerance"].get<double>();
        species_.push_back(std::make_unique<Species>(Species(name, Property("absolute tolerance", "", abs_tol))));
      }
      else{
        species_.push_back(std::make_unique<Species>(Species(name)));
      }
    }

    void ParseRelativeTolerance(const json& object){
    }

    void ParseMechanism(const json& object){
    }
  };

  template<
    template<class> class ConfigTypePolicy = JsonReaderPolicy,
    class ErrorPolicy = ThrowPolicy
  >
  class SystemBuilder : public ConfigTypePolicy<ErrorPolicy>
  {
   public:
    std::unique_ptr<micm::System> Build(std::filesystem::path path)
    {
      return this->ReadAndParse(path);
    }
  };

}  // namespace micm