/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/species.hpp>
#include <vector>

namespace micm
{

  /**
   * @brief An intraphase process
   *
   */
  template<class Rate>
  class IntraPhaseProcess
  {
   public:
    const std::vector<Species> reactants_;
    const std::vector<Species> products_;
    const Rate rate_;

   public:
    /// @brief 
    /// @param reactants 
    /// @param products 
    /// @param rate 
    IntraPhaseProcess(const std::vector<Species>& reactants, const std::vector<Species>& products, const Rate& rate);
  };

  template<class Rate>
  inline IntraPhaseProcess<Rate>::IntraPhaseProcess(const std::vector<Species>& reactants, const std::vector<Species>& products, const Rate& rate)
    : reactants_(reactants), 
      products_(products), 
      rate_(rate)
    {}

}  // namespace micm