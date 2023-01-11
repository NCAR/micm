/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

namespace micm
{

  /**
   * @brief Creates any type of intraphase process builder
   *
   */
  class IntraPhaseProcessBuilder
  {
   private:
   public:
    /// @brief Adds an additional phase to the state of the process
    /// @param phase A phase
    /// @return A reference to this object
    IntraPhaseProcessBuilder& For(const Phase& phase);
    /// @brief Adds an additional species to the process
    /// @param species A species
    /// @return A reference to this object
    IntraPhaseProcessBuilder& With(const Species& species);
    /// @brief Adds a reacting species
    /// @param reactant A Species
    /// @return A reference to this object
    IntraPhaseProcessBuilder& Reacting(const Species& reactant);
    /// @brief Adds a species that will be produced
    /// @param product A species
    /// @return A reference to this object
    IntraPhaseProcessBuilder& Producing(const Species& product);
    /// @brief Provides a yield in the amount of TODO UNITS
    /// @param yield A value
    /// @return A reference to this object
    IntraPhaseProcessBuilder& WithYield(const double yield);
    /// @brief Add a rate contant
    /// @param rate_constant A rate constant
    /// @return A reference to this object
    IntraPhaseProcessBuilder& WithRateConstant(const RateConstant& rate_constant);
    /// @brief Create the final process
    /// @return A reference to this object
    IntraPhaseProcessBuilder& Build();
  };

  inline IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::For(const Phase& phase)
  {
    return *this;
  }

  inline IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::With(const Species& species)
  {
    return *this;
  }

  inline IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::Reacting(const Species& reactant)
  {
    return *this;
  }

  inline IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::Producing(const Species& product)
  {
    return *this;
  }

  inline IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::WithYield(const double yield)
  {
    return *this;
  }

  inline IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::WithRateConstant(const RateConstant& rate_constant)
  {
    return *this;
  }

  inline IntraPhaseProcessBuilder& IntraPhaseProcessBuilder::Build()
  {
    return *this;
  }

}  // namespace micm