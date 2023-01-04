/* Copyright (C) 2022 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

namespace micm
{

  template <typename T>
  class IntraPhaseProcessBuilder
  {
   private:
   public:
    IntraPhaseProcessBuilder& For(const Phase<T>& phase);
    IntraPhaseProcessBuilder& With(const Species<T>& phase);
    IntraPhaseProcessBuilder& Build();
  };

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::For(const Phase<T>& phase)
  {
    return *this;
  }

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::With(const Species<T>& phase)
  {
    return *this;
  }

  template<typename T>
  inline IntraPhaseProcessBuilder<T>& IntraPhaseProcessBuilder<T>::Build()
  {
    return *this;
  }

}  // namespace micm