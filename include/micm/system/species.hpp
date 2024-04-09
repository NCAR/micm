// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "conditions.hpp"

namespace micm
{

  /// @brief A representation of a chemcial species
  class Species
  {
   public:
    /// @brief The name of this species
    std::string name_;

    /// @brief A list of properties of this species
    std::map<std::string, std::string> properties_string_;
    std::map<std::string, double> properties_double_;
    std::map<std::string, bool> properties_bool_;
    std::map<std::string, int> properties_int_;

    /// @brief A function that if provided will be used to parameterize
    ///        the concentration of this species during solving.
    ///        Species with this function defined will be excluded from
    ///        the solver state.
    std::function<double(const Conditions)> parameterize_{ nullptr };

    /// @brief Copy assignment
    /// @param other species to copy
    Species& operator=(const Species& other);

    /// @brief Copy Constructor
    /// @param other
    Species(const Species& other);

    /// @brief Construct a species by name only
    /// @param name The name of the species
    Species(const std::string& name);

    /// @brief Construct a species by name and properties
    /// @param name The name of the species
    /// @param properties The properties of the species
    Species(const std::string& name, const std::map<std::string, double>& properties);

    /// @brief Returns whether a species is parameterized
    bool IsParameterized() const;

    /// @brief Return the value of a species property
    template<class T>
    T GetProperty(const std::string& key) const;

    /// @brief Set the value of a species property
    template<class T>
    void SetProperty(const std::string& key, T value);

    /// @brief Return a Species instance parameterized on air density
    static Species ThirdBody();
  };

  inline Species& Species::operator=(const Species& other)
  {
    if (this == &other)
    {
      return *this;  // Handle self-assignment
    }

    name_ = other.name_;
    properties_string_ = other.properties_string_;
    properties_double_ = other.properties_double_;
    properties_int_ = other.properties_int_;
    properties_bool_ = other.properties_bool_;
    parameterize_ = other.parameterize_;

    return *this;
  }

  inline Species::Species(const Species& other)
      : name_(other.name_),
        properties_string_(other.properties_string_),
        properties_double_(other.properties_double_),
        properties_int_(other.properties_int_),
        properties_bool_(other.properties_bool_),
        parameterize_(other.parameterize_){};

  inline Species::Species(const std::string& name)
      : name_(name){};

  inline Species::Species(const std::string& name, const std::map<std::string, double>& properties)
      : name_(name),
        properties_double_(properties){};

  inline bool Species::IsParameterized() const
  {
    return parameterize_ != nullptr;
  }

  template<class T>
  inline T Species::GetProperty(const std::string& key) const
  {
    if constexpr (std::is_same<T, std::string>::value)
    {
      try {
        return properties_string_.at(key);
      } catch (const std::out_of_range& e) {
        throw std::runtime_error("Species property not found");
      }
    }
    else if constexpr (std::is_same<T, double>::value)
    {
      try {
        return properties_double_.at(key);
      } catch (const std::out_of_range& e) {
        throw std::runtime_error("Species property not found");
      }
    }
    else if constexpr (std::is_same<T, bool>::value)
    {
      try {
        return properties_bool_.at(key);
      } catch (const std::out_of_range& e) {
        throw std::runtime_error("Species property not found");
      }
    }
    else if constexpr (std::is_same<T, int>::value)
    {
      try {
        return properties_int_.at(key);
      } catch (const std::out_of_range& e) {
        throw std::runtime_error("Species property not found");
      }
    }
    else
    {
      throw std::runtime_error("Invalid type for species property");
    }
  }

  template<class T>
  inline void Species::SetProperty(const std::string& key, T value)
  {
    if constexpr (std::is_same<T, std::string>::value || std::is_same<T, const char*>::value)
    {
      properties_string_[key] = value;
    }
    else if constexpr (std::is_same<T, double>::value)
    {
      properties_double_[key] = value;
    }
    else if constexpr (std::is_same<T, bool>::value)
    {
      properties_bool_[key] = value;
    }
    else if constexpr (std::is_same<T, int>::value)
    {
      properties_int_[key] = value;
    }
    else
    {
      throw std::runtime_error("Invalid type for property");
    }
  }

  inline Species Species::ThirdBody()
  {
    Species third_body{ "M" };
    third_body.parameterize_ = [](const Conditions& c) { return c.air_density_; };
    return third_body;
  }
}  // namespace micm
