// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

constexpr micm::Real TOLERANCE = 1e-13;

TEST(TunnelingRateConstant, CalculateWithMinimalArguments)
{
  micm::Conditions conditions{
    .temperature_ = 301.24,  // [K]
  };
  micm::TunnelingRateConstantParameters tunneling_params;
  micm::Real k = micm::CalculateTunneling(tunneling_params, conditions.temperature_);
  micm::Real expected = 1.0;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(TunnelingRateConstant, CalculateWithAllArguments)
{
  micm::Real temperature = 301.24;
  micm::Conditions conditions{
    .temperature_ = temperature,  // [K]
  };
  micm::TunnelingRateConstantParameters tunneling_params{ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 };
  micm::Real k = micm::CalculateTunneling(tunneling_params, conditions.temperature_);
  micm::Real expected = 1.2 * std::exp(-2.3 / temperature) * std::exp(302.3 / std::pow(temperature, 3));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}
