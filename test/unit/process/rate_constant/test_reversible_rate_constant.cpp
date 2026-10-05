// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <type_traits>

constexpr micm::Real TOLERANCE = std::is_same_v<micm::Real, double> ? 1e-13 : 1e-6;

TEST(ReversibleRateConstant, DefaultConstructor)
{
  micm::ReversibleRateConstantParameters params{};
  micm::Conditions conditions = { .temperature_ = 298.15 };

  micm::Real k = micm::CalculateReversible(params, conditions.temperature_);
  // Default: A_ = 1, C_ = 0, k_r_ = 0 → k = 1 * exp(0/T) * 0 = 0
  micm::Real expected = 0.0;
  EXPECT_NEAR(k, expected, TOLERANCE);
}

TEST(ReversibleRateConstant, CalculateWithSystem)
{
  micm::Conditions conditions = { .temperature_ = 301.24 };

  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.14e-2;
  parameters.C_ = 2300.0;
  parameters.k_r_ = 0.32;

  micm::Real k = micm::CalculateReversible(parameters, conditions.temperature_);

  micm::Real K_eq = 1.14e-2 * std::exp(2300.0 / 301.24);
  micm::Real expected = K_eq * 0.32;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(ReversibleRateConstant, CalculateWithPrescribedArguments)
{
  micm::Real temperature = 298.15;

  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.14e-2;
  parameters.C_ = 2300.0;
  parameters.k_r_ = 0.32;

  micm::Real k = micm::CalculateReversible(parameters, temperature);

  micm::Real K_eq = 1.14e-2 * std::exp(2300.0 / temperature);
  micm::Real expected = K_eq * 0.32;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(ReversibleRateConstant, TemperatureDependence)
{
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.5;
  parameters.C_ = 2000.0;
  parameters.k_r_ = 0.5;

  micm::Real T1 = 273.15;
  micm::Real T2 = 298.15;
  micm::Real T3 = 323.15;

  micm::Real k1 = micm::CalculateReversible(parameters, T1);
  micm::Real k2 = micm::CalculateReversible(parameters, T2);
  micm::Real k3 = micm::CalculateReversible(parameters, T3);

  micm::Real expected1 = 1.5 * std::exp(2000.0 / T1) * 0.5;
  micm::Real expected2 = 1.5 * std::exp(2000.0 / T2) * 0.5;
  micm::Real expected3 = 1.5 * std::exp(2000.0 / T3) * 0.5;

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

  micm::Real k = micm::CalculateReversible(parameters, conditions.temperature_);
  EXPECT_NEAR(k, 0.0, TOLERANCE);
}

TEST(ReversibleRateConstant, EquilibriumConstantOnly)
{
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 3.5;
  parameters.C_ = 1800.0;
  parameters.k_r_ = 1.0;

  micm::Real temperature = 298.15;
  micm::Real k = micm::CalculateReversible(parameters, temperature);

  micm::Real K_eq = 3.5 * std::exp(1800.0 / temperature);
  EXPECT_NEAR(k, K_eq, TOLERANCE * K_eq);
}
