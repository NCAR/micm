/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/condition.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/process/intraphase_process.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/photolysis_rate_constant.hpp>
#include <string>
#include <vector>

namespace micm
{

  /**
   * @brief A `System` holds all physical information that represents a grid cell.
   *
   */
  class System
  {
   public:
    /// @brief The gas phase is a micm::Phase and determines what species are present.
    const Phase gas_phase_;
    /// @brief This is a catchall for anything that is not the gas phase.
    const std::vector<Phase> phases_;
    /// @brief This represents any physical measurement of a grid cell.
    const std::vector<Condition> conditions_;
    /// @brief Photolysis reactions in this system
    const IntraPhaseProcess<PhotolysisRateConstant> photolysis_reactions_;
    /// @brief Arrhenius reactions in this system
    const IntraPhaseProcess<ArrheniusRateConstant> arrhenius_reactions_;

   public:
    /// @brief Default constructor
    System();

    /// @brief 
    /// @param gas_phase A Gas phase
    System(Phase gas_phase);

    /// @brief 
    /// @param gas_phase 
    /// @param condition 
    System(Phase gas_phase, Condition condition);

    /// @brief 
    /// @param gas_phase 
    /// @param conditions 
    System(Phase gas_phase, std::vector<Condition> conditions);
  };

  inline micm::System::System()
      : gas_phase_(),
        phases_(),
        conditions_(),
        photolysis_reactions_(),
        arrhenius_reactions_()
  {
  }

  inline micm::System::System(Phase gas_phase)
      : gas_phase_(gas_phase),
        phases_(),
        conditions_(),
        photolysis_reactions_(),
        arrhenius_reactions_()
  {
  }

  inline System::System(Phase gas_phase, Condition condition)
      : gas_phase_(gas_phase),
        phases_(),
        conditions_({condition}),
        photolysis_reactions_(),
        arrhenius_reactions_()
  {
  }

  inline System::System(Phase gas_phase, std::vector<Condition> conditions)
      : gas_phase_(gas_phase),
        phases_(),
        conditions_(conditions),
        photolysis_reactions_(),
        arrhenius_reactions_()
  {
  }

}  // namespace micm
