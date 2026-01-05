#include <micm/process/rate_constant/tunneling_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

constexpr double TOLERANCE = 1e-13;

TEST(TunnelingRateConstant, CalculateWithMinimalArguments)
{
  micm::Conditions conditions{
    .temperature_ = 301.24,  // [K]
  };
  micm::TunnelingRateConstantParameters tunneling_params;
  micm::TunnelingRateConstant tunneling{ tunneling_params };
  auto k = tunneling.Calculate(conditions);
  double expected = 1.0;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(TunnelingRateConstant, CalculateWithAllArguments)
{
  double temperature = 301.24;
  micm::Conditions conditions{
    .temperature_ = temperature,  // [K]
  };
  micm::TunnelingRateConstant tunneling{ micm::TunnelingRateConstantParameters{ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 } };
  auto k = tunneling.Calculate(conditions);
  double expected = 1.2 * std::exp(-2.3 / temperature) * std::exp(302.3 / std::pow(temperature, 3));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}
