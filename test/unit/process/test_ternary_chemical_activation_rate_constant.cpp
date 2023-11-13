#include <gtest/gtest.h>

#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

TEST(TernaryChemicalActivationRateConstant, CalculateWithMinimalArugments)
{
  micm::Conditions conditions{
    .temperature_ = 301.24,  // [K]
    .air_density_ = 42.2,    // [mol mol-1]
  };
  micm::TernaryChemicalActivationRateConstantParameters ternary_params;
  ternary_params.k0_A_ = 1.0;
  ternary_params.kinf_A_ = 1.0;
  micm::TernaryChemicalActivationRateConstant ternary{ ternary_params };
  auto k = ternary.calculate(conditions);
  double k0 = 1.0;
  double kinf = 1.0;
  EXPECT_EQ(k, k0 / (1.0 + k0 * 42.2 / kinf) * std::pow(0.6, 1.0 / (1 + std::pow(std::log10(k0 * 42.2 / kinf), 2))));
}

TEST(TernaryChemicalActivationRateConstant, CalculateWithAllArugments)
{
  double temperature = 301.24;  // [K]
  micm::Conditions conditions{
    .temperature_ = temperature,
    .air_density_ = 42.2,  // [mol mol-1]
  };
  micm::TernaryChemicalActivationRateConstant ternary{ micm::TernaryChemicalActivationRateConstantParameters{
      .k0_A_ = 1.2,
      .k0_B_ = 2.3,
      .k0_C_ = 302.3,
      .kinf_A_ = 2.6,
      .kinf_B_ = -3.1,
      .kinf_C_ = 402.1,
      .Fc_ = 0.9,
      .N_ = 1.2 } };
  auto k = ternary.calculate(conditions);
  double k0 = 1.2 * std::exp(302.3 / temperature) * std::pow(temperature / 300.0, 2.3);
  double kinf = 2.6 * std::exp(402.1 / temperature) * std::pow(temperature / 300.0, -3.1);
  EXPECT_EQ(
      k, k0 / (1.0 + k0 * 42.2 / kinf) * std::pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * std::pow(std::log10(k0 * 42.2 / kinf), 2))));
}
