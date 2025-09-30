// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/species.hpp>

#include <algorithm>
#include <vector>

namespace micm
{

  /// @brief Represents a chemical species within a specific phase, 
  ///        storing the species information and its optional diffusion coefficient
  class PhaseSpecies
  {
  public:
    Species species_;
    std::optional<double> diffusion_coefficient_;

    PhaseSpecies(const Species& species)
      : species_(species)
    {}

    PhaseSpecies(const Species& species, double diffusion_coefficient)
      : species_(species), diffusion_coefficient_(diffusion_coefficient) 
    {}

    void SetDiffusionCoefficient(double diffusion_coefficient)
    {
      diffusion_coefficient_ = diffusion_coefficient;
    }
  };

  /// @brief Represents a chemical phase (e.g., gaseous, aqueous)
  ///        Each phase defines a set of species that participate in chemical
  ///        reactions within that phase.
  class Phase
  {
   public:
    std::string name_;
    /// @brief The list of phase-specific species
    std::vector<PhaseSpecies> phase_species_;

    /// @brief Defaulted constructors and assignment operators
    Phase() = default;
    Phase(const Phase&) = default;
    Phase(Phase&&) noexcept = default;
    Phase& operator=(const Phase&) = default;
    Phase& operator=(Phase&&) noexcept = default;

    /// @brief Create a phase with a name and a set of species
    Phase(const std::string& name, const std::vector<PhaseSpecies>& phase_species)
        : name_(name),
          phase_species_(phase_species)
    {
    }

    /// @brief Returns the number of non-parameterized species
    std::size_t StateSize() const
    {
      return std::count_if(phase_species_.begin(), phase_species_.end(), [](const PhaseSpecies& ps) { return !ps.species_.IsParameterized(); });
    }

    /// @brief Returns a set of unique names for each non-parameterized species
    std::vector<std::string> UniqueNames() const
    {
      std::vector<std::string> names{};
      for (const auto& phase_species : phase_species_)
        if (!phase_species.species_.IsParameterized())
          names.push_back(phase_species.species_.name_);
      return names;
    }
  };

}  // namespace micm