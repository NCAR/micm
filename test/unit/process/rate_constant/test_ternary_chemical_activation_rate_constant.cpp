// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>

#include <gtest/gtest.h>

TEST(TernaryChemicalActivationRateConstant, CalculateWithMinimalArguments)
{
  micm::Conditions conditions{
    .temperature_ = 301.24,  // [K]
    .air_density_ = 42.2,    // [mol mol-1]
  };
  micm::TernaryChemicalActivationRateConstantParameters ternary_params;
  ternary_params.k0_A_ = 1.0;
  ternary_params.kinf_A_ = 1.0;
  double k;
  micm::CalculateTernaryChemicalActivation(&ternary_params, 1, conditions.temperature_, conditions.air_density_, &k);
  double k0 = 1.0;
  double kinf = 1.0;
  EXPECT_EQ(k, k0 / (1.0 + k0 * 42.2 / kinf) * std::pow(0.6, 1.0 / (1 + std::pow(std::log10(k0 * 42.2 / kinf), 2))));
}

TEST(TernaryChemicalActivationRateConstant, CalculateWithAllArguments)
{
  double temperature = 301.24;  // [K]
  micm::Conditions conditions{
    .temperature_ = temperature,
    .air_density_ = 42.2,  // [mol mol-1]
  };
  micm::TernaryChemicalActivationRateConstantParameters params{
    .k0_A_ = 1.2,
    .k0_B_ = 2.3,
    .k0_C_ = 302.3,
    .kinf_A_ = 2.6,
    .kinf_B_ = -3.1,
    .kinf_C_ = 402.1,
    .Fc_ = 0.9,
    .N_ = 1.2
  };
  double k;
  micm::CalculateTernaryChemicalActivation(&params, 1, conditions.temperature_, conditions.air_density_, &k);
  double k0 = 1.2 * std::exp(302.3 / temperature) * std::pow(temperature / 300.0, 2.3);
  double kinf = 2.6 * std::exp(402.1 / temperature) * std::pow(temperature / 300.0, -3.1);
  EXPECT_EQ(
      k, k0 / (1.0 + k0 * 42.2 / kinf) * std::pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * std::pow(std::log10(k0 * 42.2 / kinf), 2))));
}
