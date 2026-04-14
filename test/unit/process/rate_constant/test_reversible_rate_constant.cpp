// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>

#include <gtest/gtest.h>

#include <cmath>

constexpr double TOLERANCE = 1e-13;

TEST(ReversibleRateConstant, DefaultConstructor)
{
  micm::ReversibleRateConstantParameters params{};
  micm::Conditions conditions = { .temperature_ = 298.15 };

  double k;
  micm::CalculateReversible(&params, 1, conditions.temperature_, &k);
  // Default: A_ = 1, C_ = 0, k_r_ = 0 → k = 1 * exp(0/T) * 0 = 0
  double expected = 0.0;
  EXPECT_NEAR(k, expected, TOLERANCE);
}

TEST(ReversibleRateConstant, CalculateWithSystem)
{
  micm::Conditions conditions = { .temperature_ = 301.24 };

  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.14e-2;
  parameters.C_ = 2300.0;
  parameters.k_r_ = 0.32;

  double k;
  micm::CalculateReversible(&parameters, 1, conditions.temperature_, &k);

  double K_eq = 1.14e-2 * std::exp(2300.0 / 301.24);
  double expected = K_eq * 0.32;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(ReversibleRateConstant, CalculateWithPrescribedArguments)
{
  double temperature = 298.15;

  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.14e-2;
  parameters.C_ = 2300.0;
  parameters.k_r_ = 0.32;

  double k;
  micm::CalculateReversible(&parameters, 1, temperature, &k);

  double K_eq = 1.14e-2 * std::exp(2300.0 / temperature);
  double expected = K_eq * 0.32;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(ReversibleRateConstant, TemperatureDependence)
{
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.5;
  parameters.C_ = 2000.0;
  parameters.k_r_ = 0.5;

  double T1 = 273.15;
  double T2 = 298.15;
  double T3 = 323.15;

  double k1, k2, k3;
  micm::CalculateReversible(&parameters, 1, T1, &k1);
  micm::CalculateReversible(&parameters, 1, T2, &k2);
  micm::CalculateReversible(&parameters, 1, T3, &k3);

  double expected1 = 1.5 * std::exp(2000.0 / T1) * 0.5;
  double expected2 = 1.5 * std::exp(2000.0 / T2) * 0.5;
  double expected3 = 1.5 * std::exp(2000.0 / T3) * 0.5;

  EXPECT_NEAR(k1, expected1, TOLERANCE * expected1);
  EXPECT_NEAR(k2, expected2, TOLERANCE * expected2);
  EXPECT_NEAR(k3, expected3, TOLERANCE * expected3);
}

TEST(ReversibleRateConstant, ZeroReverseRate)
{
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 5.0;
  parameters.C_ = 1000.0;
  parameters.k_r_ = 0.0;

  micm::Conditions conditions = { .temperature_ = 298.15 };

  double k;
  micm::CalculateReversible(&parameters, 1, conditions.temperature_, &k);
  EXPECT_NEAR(k, 0.0, TOLERANCE);
}

TEST(ReversibleRateConstant, EquilibriumConstantOnly)
{
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 3.5;
  parameters.C_ = 1800.0;
  parameters.k_r_ = 1.0;

  double temperature = 298.15;
  double k;
  micm::CalculateReversible(&parameters, 1, temperature, &k);

  double K_eq = 3.5 * std::exp(1800.0 / temperature);
  EXPECT_NEAR(k, K_eq, TOLERANCE * K_eq);
}
