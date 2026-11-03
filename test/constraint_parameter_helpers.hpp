// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Equilibrium constraint parameters are named "<constraint name>.<uuid>.K_eq", so tests
// cannot spell them out and must look them up instead.

/// @brief Maps each constraint's parameter names to consecutive state parameter columns, in constraint order
template<class Constraint>
std::unordered_map<std::string, micm::Index> KEqParameterIndices(const std::vector<Constraint>& constraints)
{
  std::unordered_map<std::string, micm::Index> indices;
  for (const auto& constraint : constraints)
  {
    for (const auto& name : constraint.GetParameterNames())
    {
      indices.emplace(name, indices.size());
    }
  }
  return indices;
}

/// @brief Returns the single parameter name in the map that belongs to the named constraint
template<class Map>
std::string FindConstraintParameter(const Map& parameter_map, const std::string& constraint_name)
{
  const std::string prefix = constraint_name + ".";
  std::string found;
  for (const auto& [name, index] : parameter_map)
  {
    if (name.rfind(prefix, 0) == 0)
    {
      if (!found.empty())
      {
        throw std::runtime_error("Multiple parameters found for constraint '" + constraint_name + "'");
      }
      found = name;
    }
  }
  if (found.empty())
  {
    throw std::runtime_error("No parameter found for constraint '" + constraint_name + "'");
  }
  return found;
}
