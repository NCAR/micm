/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/intraphase_process.hpp>
#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/system/condition.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <string>
#include <vector>

namespace micm
{
  struct SystemParameters
  {
    Phase gas_phase_{};
    std::vector<Phase> phases_{};
    std::vector<Condition> conditions_{};
    std::vector<IntraPhaseProcess<PhotolysisRateConstant>> photolysis_reactions_{};
    std::vector<IntraPhaseProcess<ArrheniusRateConstant>> arrhenius_reactions_{};
  };

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
    const std::vector<IntraPhaseProcess<PhotolysisRateConstant>> photolysis_reactions_;
    /// @brief Arrhenius reactions in this system
    const std::vector<IntraPhaseProcess<ArrheniusRateConstant>> arrhenius_reactions_;

   public:
    /// @brief Default constructor
    System();

    /// @brief
    /// @param gas_phase
    System(SystemParameters parameters);
  };

  inline micm::System::System()
      : gas_phase_(),
        phases_(),
        conditions_(),
        photolysis_reactions_(),
        arrhenius_reactions_()
  {
  }

  inline System::System(SystemParameters parameters)
      : gas_phase_(parameters.gas_phase_),
        phases_(parameters.phases_),
        conditions_(parameters.conditions_),
        photolysis_reactions_(parameters.photolysis_reactions_),
        arrhenius_reactions_(parameters.arrhenius_reactions_)
  {
  }

}  // namespace micm
