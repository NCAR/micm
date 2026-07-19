// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

TEST(UserDefinedRateConstant, CalculateWithSystem)
{
  micm::Real custom_params[1] = { 0.5 };
  micm::UserDefinedRateConstantData data{ .scaling_factor_ = 1.0, .custom_param_index_ = 0 };
  micm::Real k;
  k = micm::CalculateUserDefined(data, custom_params[data.custom_param_index_]);
  EXPECT_EQ(k, 0.5);
}

TEST(UserDefinedRateConstant, ConstructorWithRate)
{
  micm::Real custom_params[1] = { 1.1 };
  micm::UserDefinedRateConstantData data{ .scaling_factor_ = 1.0, .custom_param_index_ = 0 };
  micm::Real k;
  k = micm::CalculateUserDefined(data, custom_params[data.custom_param_index_]);
  EXPECT_EQ(k, 1.1);
}

TEST(UserDefinedRateConstant, ConstructorWithRateAndName)
{
  micm::Real custom_params[1] = { 1.1 };
  micm::UserDefinedRateConstantData data{ .scaling_factor_ = 1.0, .custom_param_index_ = 0 };
  micm::Real k;
  k = micm::CalculateUserDefined(data, custom_params[data.custom_param_index_]);
  EXPECT_EQ(k, 1.1);
}

TEST(UserDefinedRateConstant, CustomScalingFactor)
{
  micm::Real custom_params[1] = { 1.2 };
  micm::UserDefinedRateConstantData data{ .scaling_factor_ = 2.0, .custom_param_index_ = 0 };
  micm::Real k;
  k = micm::CalculateUserDefined(data, custom_params[data.custom_param_index_]);
  EXPECT_EQ(k, 1.2 * 2.0);
}
