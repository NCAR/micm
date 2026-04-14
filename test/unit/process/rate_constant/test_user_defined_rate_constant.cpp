// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>

#include <gtest/gtest.h>

TEST(UserDefinedRateConstant, CalculateWithSystem)
{
  double custom_params[1] = { 0.5 };
  micm::UserDefinedRateConstantData data{ .scaling_factor_ = 1.0, .custom_param_index_ = 0 };
  double k;
  micm::CalculateUserDefined(&data, 1, custom_params, &k);
  EXPECT_EQ(k, 0.5);
}

TEST(UserDefinedRateConstant, ConstructorWithRate)
{
  double custom_params[1] = { 1.1 };
  micm::UserDefinedRateConstantData data{ .scaling_factor_ = 1.0, .custom_param_index_ = 0 };
  double k;
  micm::CalculateUserDefined(&data, 1, custom_params, &k);
  EXPECT_EQ(k, 1.1);
}

TEST(UserDefinedRateConstant, ConstructorWithRateAndName)
{
  double custom_params[1] = { 1.1 };
  micm::UserDefinedRateConstantData data{ .scaling_factor_ = 1.0, .custom_param_index_ = 0 };
  double k;
  micm::CalculateUserDefined(&data, 1, custom_params, &k);
  EXPECT_EQ(k, 1.1);
}

TEST(UserDefinedRateConstant, CustomScalingFactor)
{
  double custom_params[1] = { 1.2 };
  micm::UserDefinedRateConstantData data{ .scaling_factor_ = 2.0, .custom_param_index_ = 0 };
  double k;
  micm::CalculateUserDefined(&data, 1, custom_params, &k);
  EXPECT_EQ(k, 1.2 * 2.0);
}
