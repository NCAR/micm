// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/phase.hpp>

#include <functional>
#include <string>
#include <vector>

namespace micm
{

  /// @brief Defines the gas-phase species available in the chemical system
  class System
  {
   public:
    /// @brief The gas phase, defining a set of species present in the system
    Phase gas_phase_;

    System() = default;

    System(const Phase& gas_phase)
        : gas_phase_(gas_phase)
    {
    }

    /// @brief Returns the number of gas-phase state variables
    size_t StateSize() const
    {
      return gas_phase_.StateSize();
    }

    /// @brief Returns the unique gas-phase species names
    /// @return vector of unique state variable names
    std::vector<std::string> UniqueNames() const
    {
      return UniqueNames(nullptr);
    }

    /// @brief Returns the unique gas-phase species names, optionally reordered
    /// @param f Function used to apply a specific order to unique names
    /// @return vector of unique state variable names
    std::vector<std::string> UniqueNames(
        const std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> f) const
    {
      auto names = gas_phase_.SpeciesNames();

      if (f)
      {
        std::vector<std::string> reordered;
        reordered.reserve(names.size());
        for (std::size_t i = 0; i < names.size(); ++i)
          reordered.push_back(f(names, i));
        return reordered;
      }

      return names;
    }
  };

}  // namespace micm
