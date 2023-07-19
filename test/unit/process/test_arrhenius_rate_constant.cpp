#include <gtest/gtest.h>

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

TEST(ArrheniusRateConstant, CalculateWithSystem)
{
  micm::ArrheniusRateConstant zero{};
  micm::State<micm::Matrix> state{ 0, 0, 1 };
  state.conditions_[0].temperature_ = 301.24;  // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  auto k = zero.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 1, 0.01);

  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;

  micm::ArrheniusRateConstant basic(parameters);
  k = basic.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 1, 0.01);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  micm::ArrheniusRateConstant o1d(parameters);
  k = o1d.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 2.2e-10, 0.01);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  micm::ArrheniusRateConstant hox(parameters);
  k = hox.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 3e-11 * std::exp(-200 / 301.24), 0.01);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  micm::ArrheniusRateConstant clox(parameters);
  k = clox.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 2.6e-12 * std::exp(-350 / 301.24), 0.01);
}

TEST(ArrheniusRateConstant, CalculateWithPrescribedArugments)
{
  micm::State<micm::Matrix> state{ 0, 0, 1 };
  state.conditions_[0].temperature_ = 301.24;  // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 1, 0.01);

  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;

  micm::ArrheniusRateConstant basic(parameters);
  k = basic.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 1, 0.01);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  micm::ArrheniusRateConstant o1d(parameters);
  k = o1d.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 2.2e-10, 0.01);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  micm::ArrheniusRateConstant hox(parameters);
  k = hox.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 3e-11 * std::exp(200 / 301.24), 0.01);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  micm::ArrheniusRateConstant clox(parameters);
  k = clox.calculate(state.conditions_[0], params);
  EXPECT_NEAR(k, 2.6e-12 * std::exp(-350 / 301.24), 0.01);
}
