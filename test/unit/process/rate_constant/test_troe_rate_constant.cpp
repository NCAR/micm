// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>

#include <gtest/gtest.h>

#include <iostream>
#include <stdexcept>

TEST(TroeRateConstant, CalculateWithMinimalArguments)
{
  micm::Conditions conditions{
    .temperature_ = 301.24,  // [K]
    .air_density_ = 42.2,    // [mol mol-1]
  };
  micm::TroeRateConstantParameters troe_params;
  troe_params.k0_A_ = 1.0;
  troe_params.kinf_A_ = 1.0;
  double k = micm::CalculateTroe(troe_params, conditions.temperature_, conditions.air_density_);
  double k0 = 1.0;
  double kinf = 1.0;
  EXPECT_EQ(k, 42.2 * k0 / (1.0 + 42.2 * k0 / kinf) * std::pow(0.6, 1.0 / (1 + std::pow(std::log10(42.2 * k0 / kinf), 2))));
}

TEST(TroeRateConstant, CalculateWithAllArguments)
{
  double temperature = 301.24;  // [K]
  micm::Conditions conditions{
    .temperature_ = temperature,
    .air_density_ = 42.2,  // [mol mol-1]
  };
  micm::TroeRateConstantParameters params{ .k0_A_ = 1.2,
                                           .k0_B_ = 2.3,
                                           .k0_C_ = 302.3,
                                           .kinf_A_ = 2.6,
                                           .kinf_B_ = -3.1,
                                           .kinf_C_ = 402.1,
                                           .Fc_ = 0.9,
                                           .N_ = 1.2 };
  double k = micm::CalculateTroe(params, conditions.temperature_, conditions.air_density_);
  double k0 = 1.2 * std::exp(302.3 / temperature) * std::pow(temperature / 300.0, 2.3);
  double kinf = 2.6 * std::exp(402.1 / temperature) * std::pow(temperature / 300.0, -3.1);
  EXPECT_EQ(
      k,
      42.2 * k0 / (1.0 + 42.2 * k0 / kinf) *
          std::pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * std::pow(std::log10(42.2 * k0 / kinf), 2))));
}

TEST(TroeRateConstant, AnalyticalTroeExampleAB)
{
  // based off of the troe rate constants in the analytical integration test:
  micm::Conditions conditions{
    .temperature_ = 301.24,  // [K]
    .air_density_ = 42.2,    // [mol mol-1]
  };

  micm::TroeRateConstantParameters troe_params;
  troe_params.k0_A_ = 4.0e-18;
  double k = micm::CalculateTroe(troe_params, conditions.temperature_, conditions.air_density_);

  double k_0 = 4.0e-18;
  double k_inf = 1;
  double k1 = k_0 * 42.2 / (1.0 + k_0 * 42.2 / k_inf) *
              std::pow(0.6, 1.0 / (1.0 + 1.0 / 1.0 * std::pow(std::log10(k_0 * 42.2 / k_inf), 2)));

  EXPECT_EQ(k, k1);
}

TEST(TroeRateConstant, AnalyticalTroeExampleBC)
{
  // based off of the troe rate constants in the analytical integration test:
  micm::Conditions conditions{
    .temperature_ = 301.24,  // [K]
    .air_density_ = 42.2,    // [mol mol-1]
  };

  micm::TroeRateConstantParameters troe_params{
    .k0_A_ = 1.2e-12, .k0_B_ = 167.0, .k0_C_ = 3.0, .kinf_A_ = 136.0, .kinf_B_ = 5.0, .kinf_C_ = 24.0, .Fc_ = 0.9, .N_ = 0.8
  };
  double k = micm::CalculateTroe(troe_params, conditions.temperature_, conditions.air_density_);

  double k_0 = 1.2e-12 * std::exp(3.0 / 301.24) * std::pow(301.24 / 300.0, 167.0);
  double k_inf = 136.0 * std::exp(24.0 / 301.24) * std::pow(301.24 / 300.0, 5.0);
  double k1 = k_0 * 42.2 / (1.0 + k_0 * 42.2 / k_inf) *
              std::pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * std::pow(std::log10(k_0 * 42.2 / k_inf), 2)));

  auto relative_error = std::abs(k - k1) / std::max(std::abs(k), std::abs(k1));
  if (relative_error > 1.e-14)
  {
    std::cout << "k: " << std::setprecision(15) << k << std::endl;
    std::cout << "k1: " << std::setprecision(15) << k1 << std::endl;
    std::cout << "relative_error: " << std::setprecision(15) << relative_error << std::endl;
    throw std::runtime_error("Fail to match k and k1.\n");
  }
}
