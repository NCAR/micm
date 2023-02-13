#include <micm/process/arrhenius_rate_constant.hpp>

#include <gtest/gtest.h>

TEST(ArrheniusRateConstant, DefaultConstructor){
  micm::ArrheniusRateConstant arrhenius{};
}

TEST(ArrheniusRateConstant, CalculateWithSystem){
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(micm::System());
  EXPECT_NEAR(k, 1, 0.01);

  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;

  micm::ArrheniusRateConstant basic(parameters);
  k = basic.calculate(micm::System());
  EXPECT_NEAR(k, 1, 0.01);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  micm::ArrheniusRateConstant o1d(parameters);
  k = o1d.calculate(micm::System());
  EXPECT_NEAR(k, 2.2e-10, 0.01);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  micm::ArrheniusRateConstant hox(parameters);
  k = hox.calculate(micm::System());
  // at present the system cannot be used to determine the temperature so they should be different
  EXPECT_NE(k, 5.9e-11);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  micm::ArrheniusRateConstant clox(parameters);
  k = clox.calculate(micm::System());
  // at present the system cannot be used to determine the temperature so they should be different
  EXPECT_NE(k, 8e-13);
}

TEST(ArrheniusRateConstant, CalculateWithPrescribedArugments){
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(1.0, 1.0);
  EXPECT_NEAR(k, 1, 0.01);

  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;

  micm::ArrheniusRateConstant basic(parameters);
  k = basic.calculate(1.0, 1.0);
  EXPECT_NEAR(k, 1, 0.01);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  parameters.A_ = 2.2e-10;
  micm::ArrheniusRateConstant o1d(parameters);
  k = o1d.calculate(1.0, 1.0);
  EXPECT_NEAR(k, 2.2e-10, 0.01);

  // O + HO2 -> OH + O2
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  micm::ArrheniusRateConstant hox(parameters);
  k = hox.calculate(298, 1.0);
  // at present the system cannot be used to determine the temperature so they should be different
  EXPECT_NE(k, 5.9e-11);

  // OH + HCl → H2O + Cl
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  micm::ArrheniusRateConstant clox(parameters);
  k = clox.calculate(298, 100'000);
  // at present the system cannot be used to determine the temperature so they should be different
  EXPECT_NE(k, 8e-13);
}