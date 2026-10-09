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
  micm::Real expected = k0 / (1.0 + k0 * 42.2 / kinf) * std::pow(0.6, 1.0 / (1 + std::pow(std::log10(k0 * 42.2 / kinf), 2)));
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

TEST(TernaryChemicalActivationJPL19RateConstant, CalculateWithAllArguments)
{
  micm::Real temperature = 301.24;  // [K]
  micm::Real air_density = 42.2;
  micm::TernaryChemicalActivationJPL19Parameters params{ .k0_A_ = 1.2,
                                                          .k0_B_ = 2.3,
                                                          .k0_C_ = 302.3,
                                                          .k0_D_ = 298.0,
                                                          .kinf_A_ = 2.6,
                                                          .kinf_B_ = -3.1,
                                                          .kinf_C_ = 402.1,
                                                          .kinf_D_ = 310.0,
                                                          .kint_A_ = 0.7,
                                                          .kint_B_ = 0.5,
                                                          .kint_C_ = 120.0,
                                                          .kint_D_ = 290.0,
                                                          .Fc_ = 0.9,
                                                          .N_ = 1.2 };
  micm::Real k = micm::CalculateTernaryChemicalActivationJPL19(params, temperature, air_density);
  micm::Real k0 = 1.2 * std::exp(302.3 / temperature) * std::pow(temperature / 298.0, 2.3);
  micm::Real kinf = 2.6 * std::exp(402.1 / temperature) * std::pow(temperature / 310.0, -3.1);
  micm::Real kint = 0.7 * std::exp(120.0 / temperature) * std::pow(temperature / 290.0, 0.5);
  micm::Real ratio = k0 * air_density / kinf;
  micm::Real k_f = k0 * air_density / (1.0 + ratio) * std::pow(0.9, 1.2 / (1.2 + std::pow(std::log10(ratio), 2)));
  micm::Real expected = k_f + kint * (1.0 - k_f / kinf);
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(TernaryChemicalActivationJPL19RateConstant, ZeroKintReducesToTroe)
{
  micm::Real temperature = 250.0;
  micm::Real air_density = 30.0;
  micm::TernaryChemicalActivationJPL19Parameters params;
  params.k0_A_ = 3.0e-3;
  params.k0_B_ = -2.0;
  params.kinf_A_ = 1.5;
  params.kinf_B_ = -1.0;
  params.kint_A_ = 0.0;
  micm::Real k = micm::CalculateTernaryChemicalActivationJPL19(params, temperature, air_density);
  // With k_int = 0 the result is the Troe falloff expression using the D = 298 reference temperature
  micm::Real k0 = 3.0e-3 * std::pow(temperature / 298.0, -2.0);
  micm::Real kinf = 1.5 * std::pow(temperature / 298.0, -1.0);
  micm::Real ratio = k0 * air_density / kinf;
  micm::Real expected = k0 * air_density / (1.0 + ratio) * std::pow(0.6, 1.0 / (1.0 + std::pow(std::log10(ratio), 2)));
  EXPECT_NEAR(k, expected, TOLERANCE * expected + 1e-30);
}
