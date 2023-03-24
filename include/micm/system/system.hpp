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
    IntraPhaseProcess<PhotolysisRateConstant> photolysis_reactions_{};
    IntraPhaseProcess<ArrheniusRateConstant> arrhenius_reactions_{};
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
    const IntraPhaseProcess<PhotolysisRateConstant> photolysis_reactions_;
    /// @brief Arrhenius reactions in this system
    const IntraPhaseProcess<ArrheniusRateConstant> arrhenius_reactions_;

   public:
    /// @brief Default constructor
    System();

    /// @brief Copy Constructor
    /// @param other
    System(System& other);

    /// @brief Move Constructor
    /// @param other
    System(System&& other);

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

  inline System::System(System& other)
      : gas_phase_(other.gas_phase_),
        phases_(other.phases_),
        conditions_(other.conditions_),
        photolysis_reactions_(other.photolysis_reactions_),
        arrhenius_reactions_(other.arrhenius_reactions_)
  {
  }

  inline System::System(System&& other)
      : gas_phase_(std::move(other.gas_phase_)),
        phases_(std::move(other.phases_)),
        conditions_(std::move(other.conditions_)),
        photolysis_reactions_(std::move(other.photolysis_reactions_)),
        arrhenius_reactions_(std::move(other.arrhenius_reactions_))
  {
  }

}  // namespace micm
