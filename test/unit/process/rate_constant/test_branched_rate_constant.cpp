// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/constants.hpp>

#include <gtest/gtest.h>

constexpr double TOLERANCE = 1e-13;

/// @brief Compute k0_ and z_ from BranchedRateConstantParameters as BuildFrom does.
static void ComputeDerivedFields(micm::BranchedRateConstantParameters& params)
{
  params.k0_ = 2.0e-22 * micm::constants::AVOGADRO_CONSTANT * 1.0e-6 *
               std::exp(static_cast<double>(params.n_));
  double air_ref = 2.45e19 / micm::constants::AVOGADRO_CONSTANT * 1.0e6;
  double a       = params.k0_ * air_ref;
  double b       = 0.43 * std::pow(293.0 / 298.0, -8.0);
  double A_val   = a / (1.0 + a / b) *
                   std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2.0)));
  params.z_ = A_val * (1.0 - params.a0_) / params.a0_;
}

TEST(BranchedRateConstant, CalculateAlkoxyBranchWithAllArguments)
{
  double temperature = 301.24;
  micm::Conditions conditions = {
    .temperature_ = temperature,  // [K]
    .air_density_ = 42.2          // [mol mol-1]
  };

  micm::BranchedRateConstantParameters params{
    .branch_ = micm::BranchedRateConstantParameters::Branch::Alkoxy,
    .X_ = 1.2,
    .Y_ = 204.3,
    .a0_ = 1.0e-3,
    .n_ = 2
  };
  ComputeDerivedFields(params);

  double k = micm::CalculateBranched(params, conditions.temperature_, conditions.air_density_);

  double air_dens_n_cm3 = 42.2 * micm::constants::AVOGADRO_CONSTANT * 1.0e-6;
  double a = 2.0e-22 * std::exp(2) * 2.45e19;
  double b = 0.43 * std::pow((293.0 / 298.0), -8.0);
  double z = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
  b = 0.43 * std::pow((temperature / 298.0), -8.0);
  double A = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2)));
  double expected = 1.2 * std::exp(-204.3 / temperature) * (z / (z + A));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(BranchedRateConstant, CalculateNitrateBranchWithAllArguments)
{
  double temperature = 301.24;
  micm::Conditions conditions = {
    .temperature_ = temperature,  // [K]
    .air_density_ = 42.2          // [mol mol-1]
  };

  micm::BranchedRateConstantParameters params{
    .branch_ = micm::BranchedRateConstantParameters::Branch::Nitrate,
    .X_ = 1.2,
    .Y_ = 204.3,
    .a0_ = 1.0e-3,
    .n_ = 2
  };
  ComputeDerivedFields(params);

  double k = micm::CalculateBranched(params, conditions.temperature_, conditions.air_density_);

  double air_dens_n_cm3 = 42.2 * micm::constants::AVOGADRO_CONSTANT * 1.0e-6;
  double a = 2.0e-22 * std::exp(2) * 2.45e19;
  double b = 0.43 * std::pow((293.0 / 298.0), -8.0);
  double z = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
  b = 0.43 * std::pow((temperature / 298.0), -8.0);
  double A = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2)));
  double expected = 1.2 * std::exp(-204.3 / temperature) * (A / (z + A));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}
