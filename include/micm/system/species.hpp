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
    std::map<std::string, double> properties_;

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
    properties_ = other.properties_;  // This performs a shallow copy
    parameterize_ = other.parameterize_;

    return *this;
  }

  inline Species::Species(const Species& other)
      : name_(other.name_),
        properties_(other.properties_),
        parameterize_(other.parameterize_){};

  inline Species::Species(const std::string& name)
      : name_(name){};

  inline Species::Species(const std::string& name, const std::map<std::string, double>& properties)
      : name_(name),
        properties_(properties){};

  inline bool Species::IsParameterized() const
  {
    return parameterize_ != nullptr;
  }

  inline Species Species::ThirdBody()
  {
    Species third_body{ "M" };
    third_body.parameterize_ = [](const Conditions& c) { return c.air_density_; };
    return third_body;
  }
}  // namespace micm
