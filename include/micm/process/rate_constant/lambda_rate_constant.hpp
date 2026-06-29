// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/system/conditions.hpp>

#include <functional>
#include <string>

namespace micm
{
  struct LambdaRateConstantParameters
  {
    /// @brief Label for the reaction used to identify user-defined parameters
    std::string label_;
    /// @brief Lambda function for calculating the rate constant
    std::function<double(const Conditions&)> lambda_function_;
  };
}  // namespace micm
