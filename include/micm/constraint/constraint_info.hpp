// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <cstddef>
#include <vector>

namespace micm
{
  /// @brief Information for each constraint (built during ConstraintSet construction)
  struct ConstraintInfo
  {
    Index index_;                             // Index in constraints_ vector
    Index row_index_;                         // Row in the forcing/Jacobian
    Index number_of_dependencies_;            // Number of species this constraint depends on
    Index dependency_offset_;                 // Starting offset in dependency_ids_
    Index jacobian_flat_offset_;              // Starting offset in jacobian_flat_ids_
    std::vector<Index> state_indices_;        // Dependency indices in state_variables_
    std::vector<Index> state_param_indices_;  // Dependency indices in custom_rate_parameters_
  };
}  // namespace micm