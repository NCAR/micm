// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>
#include <cstddef>

namespace micm
{
/// @brief Information for each constraint (built during ConstraintSet construction)
  struct ConstraintInfo
  {
    std::size_t index_;                       // Index in constraints_ vector
    std::size_t row_index_;                   // Row in the forcing/Jacobian
    std::size_t number_of_dependencies_;      // Number of species this constraint depends on
    std::size_t dependency_offset_;           // Starting offset in dependency_ids_
    std::size_t jacobian_flat_offset_;        // Starting offset in jacobian_flat_ids_
    std::vector<std::size_t> state_indices_;  // Dependency indices in state_variables_
  };
}