/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

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
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::Build()
  {
    return *this;
  }

}  // namespace micm