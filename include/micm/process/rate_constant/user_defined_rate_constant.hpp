// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <string>

namespace micm
{
  struct UserDefinedRateConstantParameters
  {
    /// @brief Label for the reaction used to identify user-defined parameters
    std::string label_;
    /// @brief Scaling factor to apply to user-provided rate constants
    double scaling_factor_{ 1.0 };
  };

  /// @brief GPU-safe calculation data for a user-defined rate constant.
  ///        Populated by ReactionRateConstantStore::BuildFrom from UserDefinedRateConstantParameters;
  ///        do not construct directly.
  struct UserDefinedRateConstantData
  {
    /// @brief Scaling factor applied to the user-provided rate constant value
    double scaling_factor_{ 1.0 };
    /// @brief Index into custom_rate_parameters_[cell] holding the rate constant value
    std::size_t custom_param_index_{ 0 };
  };
}  // namespace micm
