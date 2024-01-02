/* Copyright (C) 2023-2024 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <algorithm>
#include <micm/system/species.hpp>
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

   public:
    /// @brief Default constructor
    Phase()
        : species_()
    {
    }

    /// @brief Create a phase with a set of species
    /// @param species A unique list of species
    Phase(const std::vector<Species>& species)
        : species_(species)
    {
    }
    Phase(const Phase& other)
        : species_(other.species_)
    {
    }

    Phase(Phase&& phase) noexcept
        : species_(std::move(phase.species_))
    {
    }

    Phase& operator=(const Phase& other)
    {
      species_ = other.species_;
      return *this;
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