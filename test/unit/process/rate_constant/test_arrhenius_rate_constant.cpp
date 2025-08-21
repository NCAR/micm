#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

constexpr double TOLERANCE = 1e-13;

TEST(ArrheniusRateConstant, CalculateWithSystem)
{
  micm::ArrheniusRateConstant zero{};
  micm::Conditions conditions = {
    .temperature_ = 301.24  // [K]
  };

  auto k = zero.Calculate(conditions);
  EXPECT_NEAR(k, 1, TOLERANCE);

  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;

  micm::ArrheniusRateConstant basic(parameters);
  k = basic.Calculate(conditions);
  EXPECT_NEAR(k, 1, TOLERANCE);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  micm::ArrheniusRateConstant o1d(parameters);
  k = o1d.Calculate(conditions);
  EXPECT_NEAR(k, 2.2e-10, TOLERANCE);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  micm::ArrheniusRateConstant hox(parameters);
  k = hox.Calculate(conditions);
  EXPECT_NEAR(k, 3e-11 * std::exp(-200 / 301.24), TOLERANCE);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  micm::ArrheniusRateConstant clox(parameters);
  k = clox.Calculate(conditions);
  EXPECT_NEAR(k, 2.6e-12 * std::exp(-350 / 301.24), TOLERANCE);
}

TEST(ArrheniusRateConstant, CalculateWithPrescribedArugments)
{
  micm::Conditions conditions = {
    .temperature_ = 301.24  // [K]
  };

  micm::ArrheniusRateConstant zero{};
  auto k = zero.Calculate(conditions);
  EXPECT_NEAR(k, 1, TOLERANCE);

  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;

  micm::ArrheniusRateConstant basic(parameters);
  k = basic.Calculate(conditions);
  EXPECT_NEAR(k, 1, TOLERANCE);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  micm::ArrheniusRateConstant o1d(parameters);
  k = o1d.Calculate(conditions);
  EXPECT_NEAR(k, 2.2e-10, TOLERANCE);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  micm::ArrheniusRateConstant hox(parameters);
  k = hox.Calculate(conditions);
  EXPECT_NEAR(k, 3e-11 * std::exp(-200 / 301.24), TOLERANCE);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  micm::ArrheniusRateConstant clox(parameters);
  k = clox.Calculate(conditions);
  EXPECT_NEAR(k, 2.6e-12 * std::exp(-350 / 301.24), TOLERANCE);
}
