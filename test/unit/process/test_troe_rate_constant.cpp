#include <gtest/gtest.h>

#include <micm/process/troe_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

TEST(TroeRateConstant, CalculateWithMinimalArugments)
{
  micm::State<micm::Matrix> state{ 0, 0, 1 };
  state.conditions_[0].temperature_ = 301.24;  // [K]
  state.conditions_[0].air_density_ = 42.2;    // [mol mol-1]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::TroeRateConstantParameters troe_params;
  troe_params.k0_A_ = 1.0;
  troe_params.kinf_A_ = 1.0;
  micm::TroeRateConstant troe{ troe_params };
  auto k = troe.calculate(state.conditions_[0], params);
  double k0 = 1.0;
  double kinf = 1.0;
  EXPECT_EQ(k, 42.2 * k0 / (1.0 + 42.2 * k0 / kinf) * pow(0.6, 1.0 / (1 + pow(log10(42.2 * k0 / kinf), 2))));
}

TEST(TroeRateConstant, CalculateWithAllArugments)
{
  micm::State<micm::Matrix> state{ 0, 0, 1 };
  double temperature = 301.24;
  state.conditions_[0].temperature_ = temperature;  // [K]
  state.conditions_[0].air_density_ = 42.2;         // [mol mol-1]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::TroeRateConstant troe{ micm::TroeRateConstantParameters{ .k0_A_ = 1.2,
                                                                 .k0_B_ = 2.3,
                                                                 .k0_C_ = 302.3,
                                                                 .kinf_A_ = 2.6,
                                                                 .kinf_B_ = -3.1,
                                                                 .kinf_C_ = 402.1,
                                                                 .Fc_ = 0.9,
                                                                 .N_ = 1.2 } };
  auto k = troe.calculate(state.conditions_[0], params);
  double k0 = 1.2 * exp(302.3 / temperature) * pow(temperature / 300.0, 2.3);
  double kinf = 2.6 * exp(402.1 / temperature) * pow(temperature / 300.0, -3.1);
  EXPECT_EQ(k, 42.2 * k0 / (1.0 + 42.2 * k0 / kinf) * pow(0.9, 1.0 / (1.0 + 1.0 / 1.2 * pow(log10(42.2 * k0 / kinf), 2))));
}

TEST(TroeRateConstant, AnalyticalTroeExampleAB)
{
  // based off of the troe rate constants in the analytical integration test:
  micm::State<micm::Matrix> state{ 0, 0, 1 };
  state.conditions_[0].temperature_ = 301.24;  // [K]
  state.conditions_[0].air_density_ = 42.2;    // [mol mol-1]

  auto params = state.custom_rate_parameters_[0].begin();

  micm::TroeRateConstantParameters troe_params;
  troe_params.k0_A_ = 4.0e-18;
  micm::TroeRateConstant troe{ troe_params };

  auto k = troe.calculate(state.conditions_[0], params);

  double k_0 = 4.0e-18;
  double k_inf = 1;
  double k1 =
      k_0 * 42.2 / (1.0 + k_0 * 42.2 / k_inf) * pow(0.6, 1.0 / (1.0 + 1.0 / 1.0 * pow(log10(k_0 * 42.2 / k_inf), 2)));

  EXPECT_EQ(k, k1);
}

TEST(TroeRateConstant, AnalyticalTroeExampleBC)
{
  // based off of the troe rate constants in the analytical integration test:
  micm::State<micm::Matrix> state{ 0, 0, 1 };
  state.conditions_[0].temperature_ = 301.24;  // [K]
  state.conditions_[0].air_density_ = 42.2;    // [mol mol-1]

  auto params = state.custom_rate_parameters_[0].begin();

  micm::TroeRateConstantParameters troe_params{
    .k0_A_ = 1.2e-12, .k0_B_ = 167, .k0_C_ = 3, .kinf_A_ = 136, .kinf_B_ = 5, .kinf_C_ = 24, .Fc_ = 0.9, .N_ = 0.8
  };
  micm::TroeRateConstant troe{ troe_params };

  auto k = troe.calculate(state.conditions_[0], params);

  double k_0 = 1.2e-12 * exp(3.0 / 301.24) * pow(301.24 / 300.0, 167.0);
  double k_inf = 136.0 * exp(24.0 / 301.24) * pow(301.24 / 300.0, 5.0);
  double k1 =
      k_0 * 42.2 / (1.0 + k_0 * 42.2 / k_inf) * pow(0.9, 1.0 / (1.0 + (1.0 / 0.8) * pow(log10(k_0 * 42.2 / k_inf), 2)));

  EXPECT_EQ(k, k1);
}