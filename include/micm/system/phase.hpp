// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/species.hpp>

#include <algorithm>
#include <vector>

namespace micm
{

  /**
   * @brief A class which represnts a certain chemical phase (gaseous, aquous)
   *
   * Each phase represents a set of species that participate in chemical reactions in that phase.
   */
  class Phase
  {
   public:
    /// @brief The list of species
    std::vector<Species> species_;
    std::string name_;

    /// @brief Defaulted constructors and assignment operators
    Phase() = default;
    Phase(const Phase&) = default;
    Phase(Phase&&) noexcept = default;
    Phase& operator=(const Phase&) = default;
    Phase& operator=(Phase&&) noexcept = default;

    /// @brief Create a phase with a set of species
     Phase(const std::vector<Species>& species)
        : name_(),
          species_(species)
    {
    }

    /// @brief Create a phase with a name and a set of species
    Phase(const std::string& name, const std::vector<Species>& species)
        : name_(name),
          species_(species)
    {
    }

    /// @brief Returns the number of non-parameterized species
    std::size_t StateSize() const
    {
      return std::count_if(species_.begin(), species_.end(), [](const Species& s) { return !s.IsParameterized(); });
    }

    /// @brief Returns a set of unique names for each non-parameterized species
    std::vector<std::string> UniqueNames() const
    {
      std::vector<std::string> names{};
      for (const auto& species : species_)
        if (!species.IsParameterized())
          names.push_back(species.name_);
      return names;
    }
  };
}  // namespace micm