// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <string>
#include <vector>

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
  };

  inline Species& Species::operator=(const Species& other)
  {
      if (this == &other) {
       return *this; // Handle self-assignment
     }

     name_ = other.name_;
     properties_ = other.properties_; // This performs a shallow copy

     return *this;
  }

  inline Species::Species(const Species& other)
      : name_(other.name_),
        properties_(other.properties_){};

  inline Species::Species(const std::string& name)
      : name_(name){};

  inline Species::Species(const std::string& name, const std::map<std::string, double>& properties)
      : name_(name),
        properties_(properties){};

}  // namespace micm
