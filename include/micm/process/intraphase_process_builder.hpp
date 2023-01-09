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
   * @tparam T The underlying datatype of the species and phases
   */
  template <typename T>
  class IntraPhaseProcessBuilder
  {
   private:
   public:
    /// @brief Adds an additional phase to the state of the process
    /// @param phase A phase
    /// @return A reference to this object
    IntraPhaseProcessBuilder& For(const Phase<T>& phase);
    /// @brief Adds an additional species to the process
    /// @param species A species 
    /// @return A reference to this object
    IntraPhaseProcessBuilder& With(const Species<T>& species);
    /// @brief Adds a reacting species
    /// @param reactant A Species
    /// @return A reference to this object
    IntraPhaseProcessBuilder& Reacting(const Species<T>& reactant);
    /// @brief Adds a species that will be produced
    /// @param product A species
    /// @return A reference to this object
    IntraPhaseProcessBuilder& Producing(const Species<T>& product);
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

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::For(const Phase<T>& phase)
  {
    return *this;
  }

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::With(const Species<T>& species)
  {
    return *this;
  }

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::Reacting(const Species<T>& reactant)
  {
    return *this;
  }

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::Producing(const Species<T>& product)
  {
    return *this;
  }

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::WithYield(const double yield)
  {
    return *this;
  }

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::WithRateConstant(const RateConstant& rate_constant)
  {
    return *this;
  }

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::Build()
  {
    return *this;
  }

}  // namespace micm