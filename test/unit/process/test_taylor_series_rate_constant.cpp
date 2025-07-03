#include <micm/process/taylor_series_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

TEST(TaylorSeriesRateConstant, CalculateWithSystem)
{
  micm::TaylorSeriesRateConstant zero{};
  micm::Conditions conditions = {
    .temperature_ = 301.24  // [K]
  };

  auto k = zero.Calculate(conditions);
  EXPECT_NEAR(k, 1, 0.01);

  micm::TaylorSeriesRateConstantParameters parameters;
  parameters.A_ = 1;

  micm::TaylorSeriesRateConstant basic(parameters);
  k = basic.Calculate(conditions);
  EXPECT_NEAR(k, 1, 0.01);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  micm::TaylorSeriesRateConstant o1d(parameters);
  k = o1d.Calculate(conditions);
  EXPECT_NEAR(k, 2.2e-10, 0.01);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  parameters.coefficients_ = { 12.5, 1.3e-2, 5.2e-4 };  // Taylor series coefficients
  micm::TaylorSeriesRateConstant hox(parameters);
  k = hox.Calculate(conditions);
  EXPECT_NEAR(k, 3e-11 * std::exp(-200 / 301.24) * (12.5 + 1.3e-2 * 301.24 + 5.2e-4 * std::pow(301.24, 2)), 0.01);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  parameters.coefficients_ = { 1.0, 4.3e-1, 7.3e-3 };  // Taylor series coefficients
  micm::TaylorSeriesRateConstant clox(parameters);
  k = clox.Calculate(conditions);
  EXPECT_NEAR(k, 2.6e-12 * std::exp(-350 / 301.24) * (1.0 + 4.3e-1 * 301.24 + 7.3e-3 * std::pow(301.24, 2)), 0.01);
}

TEST(TaylorSeriesRateConstant, CalculateWithPrescribedArugments)
{
  micm::Conditions conditions = {
    .temperature_ = 301.24  // [K]
  };

  micm::TaylorSeriesRateConstant zero{};
  auto k = zero.Calculate(conditions);
  EXPECT_NEAR(k, 1, 0.01);

  micm::TaylorSeriesRateConstantParameters parameters;
  parameters.A_ = 1;

  micm::TaylorSeriesRateConstant basic(parameters);
  k = basic.Calculate(conditions);
  EXPECT_NEAR(k, 1, 0.01);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  micm::TaylorSeriesRateConstant o1d(parameters);
  k = o1d.Calculate(conditions);
  EXPECT_NEAR(k, 2.2e-10, 0.01);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  parameters.coefficients_ = { 12.5, 1.3e-2, 5.2e-4 };  // Taylor series coefficients
  micm::TaylorSeriesRateConstant hox(parameters);
  k = hox.Calculate(conditions);
  EXPECT_NEAR(k, 3e-11 * std::exp(200 / 301.24) * (12.5 + 1.3e-2 * 301.24 + 5.2e-4 * std::pow(301.24, 2)), 0.01);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  parameters.coefficients_ = { 1.0, 4.3e-1, 7.3e-3 };  // Taylor series coefficients
  micm::TaylorSeriesRateConstant clox(parameters);
  k = clox.Calculate(conditions);
  EXPECT_NEAR(k, 2.6e-12 * std::exp(-350 / 301.24) * (1.0 + 4.3e-1 * 301.24 + 7.3e-3 * std::pow(301.24, 2)), 0.01);
}
