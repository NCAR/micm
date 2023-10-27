#include <gtest/gtest.h>

#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

TEST(TunnelingRateConstant, CalculateWithMinimalArugments)
{
  micm::Conditions conditions{
    .temperature_ = 301.24,  // [K]
  };
  micm::TunnelingRateConstantParameters tunneling_params;
  micm::TunnelingRateConstant tunneling{ tunneling_params };
  auto k = tunneling.calculate(conditions);
  EXPECT_NEAR(k, 1.0, 1.0e-8);
}

TEST(TunnelingRateConstant, CalculateWithAllArugments)
{
  double temperature = 301.24;
  micm::Conditions conditions{
    .temperature_ = temperature,  // [K]
  };
  micm::TunnelingRateConstant tunneling{ micm::TunnelingRateConstantParameters{ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 } };
  auto k = tunneling.calculate(conditions);
  EXPECT_NEAR(k, 1.2 * std::exp(-2.3 / temperature) * std::exp(302.3 / std::pow(temperature, 3)), 1.0e-8);
}
