/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/property.hpp>
#include <string>
#include <vector>

namespace micm
{

  /**
   * @brief A representation of a chemcial species
   *
   */
  class Species
  {
   public:
    /// @brief The name of this species
    const std::string name_;
    /// @brief A list of properties of this species
    const std::vector<Property> properties_;

    /// @brief Default constructor is not allowed
    Species() = delete;

    /// @brief Copy Constructor
    /// @param other
    Species(const Species& other);

    /// @brief Construct a species by name only
    /// @param name The name of the species
    Species(const std::string& name);

    /// @brief Construct a species by name and properties
    /// @param name The name of the species
    /// @param properties The properties of the species
    Species(const std::string& name, std::vector<Property> properties);

    /// @brief Construct a species by name and property
    /// @param name The name of the species
    /// @param property A property of the species
    Species(const std::string& name, Property property);
  };

  inline Species::Species(const Species& other)
      : name_(other.name_),
        properties_(other.properties_){};

  inline Species::Species(const std::string& name)
      : name_(name){};

  inline Species::Species(const std::string& name, std::vector<Property> properties)
      : name_(name),
        properties_(properties){};

  inline Species::Species(const std::string& name, Property property)
      : name_(name),
        properties_({ property }){};

}  // namespace micm
