// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/lambda_rate_constant.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

constexpr micm::Real TOLERANCE = 1e-13;

TEST(LambdaRateConstant, CalculateWithSystem)
{
  micm::LambdaRateConstantParameters params{ .label_ = "test lambda",
                                             .lambda_function_ = [](const micm::Conditions& conditions)
                                             { return 2.0 * conditions.temperature_; } };
  micm::Conditions conditions = {
    .temperature_ = 300.0  // [K]
  };
  micm::Real k = params.lambda_function_(conditions);
  micm::Real expected = 600.0;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}
