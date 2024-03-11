#include <gtest/gtest.h>

#include <micm/process/branched_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
#include <micm/util/constants.hpp>

TEST(BranchedRateConstant, CalculateAlkoxyBranchWithAllArugments)
{
  double temperature = 301.24;
  micm::Conditions conditions = {
    .temperature_ = temperature,  // [K]
    .air_density_ = 42.2          // [mol mol-1]
  };

  micm::BranchedRateConstant branched{ micm::BranchedRateConstantParameters{
      .branch_ = micm::BranchedRateConstantParameters::Branch::Alkoxy, .X_ = 1.2, .Y_ = 204.3, .a0_ = 1.0e-3, .n_ = 2 } };
  auto k = branched.calculate(conditions);
  double air_dens_n_cm3 = 42.2 * AVOGADRO_CONSTANT * 1.0e-6;
  double a = 2.0e-22 * std::exp(2) * 2.45e19;
  double b = 0.43 * std::pow((293.0 / 298.0), -8.0);
  double z = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
  b = 0.43 * std::pow((temperature / 298.0), -8.0);
  double A = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2)));

  EXPECT_NEAR(k, 1.2 * std::exp(-204.3 / temperature) * (z / (z + A)), 1.0e-8);
}

TEST(BranchedRateConstant, CalculateNitrateBranchWithAllArugments)
{
  double temperature = 301.24;
  micm::Conditions conditions = {
    .temperature_ = temperature,  // [K]
    .air_density_ = 42.2          // [mol mol-1]
  };

  micm::BranchedRateConstant branched{ micm::BranchedRateConstantParameters{
      .branch_ = micm::BranchedRateConstantParameters::Branch::Nitrate, .X_ = 1.2, .Y_ = 204.3, .a0_ = 1.0e-3, .n_ = 2 } };
  auto k = branched.calculate(conditions);
  double air_dens_n_cm3 = 42.2 * AVOGADRO_CONSTANT * 1.0e-6;
  double a = 2.0e-22 * std::exp(2) * 2.45e19;
  double b = 0.43 * std::pow((293.0 / 298.0), -8.0);
  double z = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2))) * (1.0 - 1.0e-3) / 1.0e-3;
  a = 2.0e-22 * std::exp(2) * air_dens_n_cm3;
  b = 0.43 * std::pow((temperature / 298.0), -8.0);
  double A = a / (1.0 + a / b) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a / b), 2)));

  EXPECT_NEAR(k, 1.2 * std::exp(-204.3 / temperature) * (A / (z + A)), 1.0e-8);
}
