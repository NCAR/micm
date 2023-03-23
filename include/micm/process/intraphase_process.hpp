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
    std::vector<Species> reactants_;
    std::vector<Species> products_;
    Rate rate_;

   public:
    /// @brief Default constructor
    IntraPhaseProcess();

    /// @brief 
    /// @param reactants 
    /// @param products 
    /// @param rate 
    IntraPhaseProcess(std::vector<Species> reactants, std::vector<Species> products, Rate rate);
  };

  template<class Rate>
  inline IntraPhaseProcess<Rate>::IntraPhaseProcess()
      : reactants_(),
        products_(),
        rate_()
  {
  }

  template<class Rate>
  inline IntraPhaseProcess<Rate>::IntraPhaseProcess(std::vector<Species> reactants, std::vector<Species> products, Rate rate)
      : reactants_(reactants),
        products_(products),
        rate_(rate)
  {

  }

}  // namespace micm