#include <gtest/gtest.h>

#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

TEST(TunnelingRateConstant, CalculateWithMinimalArugments)
{
  micm::State<micm::Matrix> state{ 0, 0, 1 };
  state.conditions_[0].temperature_ = 301.24;  // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::TunnelingRateConstantParameters tunneling_params;
  micm::TunnelingRateConstant tunneling{ tunneling_params };
  auto k = tunneling.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 1.0, 1.0e-8);
}

TEST(TunnelingRateConstant, CalculateWithAllArugments)
{
  micm::State<micm::Matrix> state{ 0, 0, 1 };
  double temperature = 301.24;
  state.conditions_[0].temperature_ = temperature;  // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::TunnelingRateConstant tunneling{ micm::TunnelingRateConstantParameters{ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 } };
  auto k = tunneling.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 1.2 * std::exp(-2.3 / temperature) * std::exp(302.3 / std::pow(temperature, 3)), 1.0e-8);
}
