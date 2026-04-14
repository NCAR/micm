#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/process/rate_constant/taylor_series_rate_constant.hpp>
#include <micm/system/conditions.hpp>

#include <gtest/gtest.h>

constexpr double TOLERANCE = 1e-13;

TEST(TaylorSeriesRateConstant, CalculateWithSystem)
{
  micm::Conditions conditions = {
    .temperature_ = 301.24  // [K]
  };

  micm::TaylorSeriesRateConstantParameters params{};
  double k;

  // Default parameters: A=1, n_coefficients=1, coefficients[0]=1 → k=1
  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  double expected = 1;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  params.A_ = 1;
  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  params.A_ = 2.2e-10;
  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  expected = 2.2e-10;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // O + HO2 -> OH + O2
  params.A_               = 3e-11;
  params.C_               = -200;
  params.coefficients_[0] = 12.5;
  params.coefficients_[1] = 1.3e-2;
  params.coefficients_[2] = 5.2e-4;
  params.n_coefficients_  = 3;
  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  expected = 3e-11 * std::exp(-200 / 301.24) * (12.5 + 1.3e-2 * 301.24 + 5.2e-4 * std::pow(301.24, 2));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // OH + HCl → H2O + Cl
  params.A_               = 2.6e-12;
  params.C_               = -350;
  params.coefficients_[0] = 1.0;
  params.coefficients_[1] = 4.3e-1;
  params.coefficients_[2] = 7.3e-3;
  params.n_coefficients_  = 3;
  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  expected = 2.6e-12 * std::exp(-350 / 301.24) * (1.0 + 4.3e-1 * 301.24 + 7.3e-3 * std::pow(301.24, 2));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(TaylorSeriesRateConstant, CalculateWithPrescribedArguments)
{
  micm::Conditions conditions = {
    .temperature_ = 301.24  // [K]
  };

  micm::TaylorSeriesRateConstantParameters params{};
  double k;

  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  double expected = 1;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  params.A_ = 1;
  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  params.A_ = 2.2e-10;
  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  expected = 2.2e-10;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // O + HO2 -> OH + O2
  params.A_               = 3e-11;
  params.C_               = -200;
  params.coefficients_[0] = 12.5;
  params.coefficients_[1] = 1.3e-2;
  params.coefficients_[2] = 5.2e-4;
  params.n_coefficients_  = 3;
  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  expected = 3e-11 * std::exp(-200 / 301.24) * (12.5 + 1.3e-2 * 301.24 + 5.2e-4 * std::pow(301.24, 2));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);

  // OH + HCl → H2O + Cl
  params.A_               = 2.6e-12;
  params.C_               = -350;
  params.coefficients_[0] = 1.0;
  params.coefficients_[1] = 4.3e-1;
  params.coefficients_[2] = 7.3e-3;
  params.n_coefficients_  = 3;
  micm::CalculateTaylorSeries(&params, 1, conditions.temperature_, conditions.pressure_, &k);
  expected = 2.6e-12 * std::exp(-350 / 301.24) * (1.0 + 4.3e-1 * 301.24 + 7.3e-3 * std::pow(301.24, 2));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}
