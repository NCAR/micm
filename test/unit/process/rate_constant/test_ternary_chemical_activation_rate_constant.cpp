// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <type_traits>

// double mode keeps the original exact-equality check; float mode allows a few ULPs
constexpr micm::Real TOLERANCE = std::is_same_v<micm::Real, double> ? 0.0 : 1e-6;

TEST(TernaryChemicalActivationRateConstant, CalculateWithMinimalArguments)
{
  micm::Conditions conditions{
    .temperature_ = 301.24,  // [K]
    .air_density_ = 42.2,    // [mol mol-1]
  };
  micm::TernaryChemicalActivationRateConstantParameters ternary_params;
  ternary_params.k0_A_ = 1.0;
  ternary_params.kinf_A_ = 1.0;
  micm::Real k = micm::CalculateTernaryChemicalActivation(ternary_params, conditions.temperature_, conditions.air_density_);
  micm::Real k0 = 1.0;
  micm::Real kinf = 1.0;
  micm::Real expected =
      k0 / (1.0 + k0 * 42.2 / kinf) * std::pow(0.6, 1.0 / (1 + std::pow(std::log10(k0 * 42.2 / kinf), 2)));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(TernaryChemicalActivationRateConstant, CalculateWithAllArguments)
{
  micm::Real temperature = 301.24;  // [K]
  micm::Conditions conditions{
    .temperature_ = temperature,
    .air_density_ = 42.2,  // [mol mol-1]
  };
  micm::TernaryChemicalActivationRateConstantParameters params{
    .k0_A_ = 1.2, .k0_B_ = 2.3, .k0_C_ = 302.3, .kinf_A_ = 2.6, .kinf_B_ = -3.1, .kinf_C_ = 402.1, .Fc_ = 0.9, .N_ = 1.2
  };
  micm::Real k = micm::CalculateTernaryChemicalActivation(params, conditions.temperature_, conditions.air_density_);
  micm::Real k0 = 1.2 * std::exp(302.3 / temperature) * std::pow(temperature / 300.0, 2.3);
  micm::Real kinf = 2.6 * std::exp(402.1 / temperature) * std::pow(temperature / 300.0, -3.1);
  micm::Real expected =
      k0 / (1.0 + k0 * 42.2 / kinf) * std::pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * std::pow(std::log10(k0 * 42.2 / kinf), 2)));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}
