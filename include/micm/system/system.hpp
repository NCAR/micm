/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <string>
#include <vector>

namespace micm
{
  // TODO: jiwon 5/18 - commented out because there is no default constructor for Phase class
  // struct SystemParameters
  // {
  //   /// @brief The gas phase is a micm::Phase and determines what species are present.
  //   Phase gas_phase_{};
  //   /// @brief This is a catchall for anything that is not the gas phase.
  //   std::vector<Phase> phases_{};
  // };

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

   public:
    /// @brief Default constructor is not allowed
    System() = delete;

    /// @brief
    /// @param gas_phase
    System(const Phase& gas_phase, const std::vector<Phase>& phases)
      : gas_phase_(gas_phase),
        phases_(phases)
      {}
  };

}  // namespace micm
