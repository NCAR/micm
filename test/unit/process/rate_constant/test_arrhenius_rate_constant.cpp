// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>

#include <gtest/gtest.h>

constexpr double TOLERANCE = 1e-13;

TEST(ArrheniusRateConstant, CalculateWithSystem)
{
  micm::Conditions conditions = {
    .temperature_ = 301.24  // [K]
  };

  micm::ArrheniusRateConstantParameters zero{};
  double k = micm::CalculateArrhenius(zero, conditions.temperature_, conditions.pressure_);
  double expected = 1;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = expected;
  k = micm::CalculateArrhenius(parameters, conditions.temperature_, conditions.pressure_);
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  k = micm::CalculateArrhenius(parameters, conditions.temperature_, conditions.pressure_);
  expected = 2.2e-10;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  k = micm::CalculateArrhenius(parameters, conditions.temperature_, conditions.pressure_);
  expected = 3e-11 * std::exp(-200 / 301.24);
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  k = micm::CalculateArrhenius(parameters, conditions.temperature_, conditions.pressure_);
  expected = 2.6e-12 * std::exp(-350 / 301.24);
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(ArrheniusRateConstant, CalculateWithPrescribedArguments)
{
  micm::Conditions conditions = {
    .temperature_ = 301.24  // [K]
  };

  micm::ArrheniusRateConstantParameters zero{};
  double k = micm::CalculateArrhenius(zero, conditions.temperature_, conditions.pressure_);
  double expected = 1;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;
  k = micm::CalculateArrhenius(parameters, conditions.temperature_, conditions.pressure_);
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  k = micm::CalculateArrhenius(parameters, conditions.temperature_, conditions.pressure_);
  expected = 2.2e-10;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  k = micm::CalculateArrhenius(parameters, conditions.temperature_, conditions.pressure_);
  expected = 3e-11 * std::exp(-200 / 301.24);
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  k = micm::CalculateArrhenius(parameters, conditions.temperature_, conditions.pressure_);
  expected = 2.6e-12 * std::exp(-350 / 301.24);
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}
