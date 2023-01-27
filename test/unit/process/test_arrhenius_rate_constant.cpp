#include <micm/process/arrhenius_rate_constant.hpp>

#include <gtest/gtest.h>

TEST(ArrheniusRateConstant, DefaultConstructor){
  micm::ArrheniusRateConstant arrhenius{};
}

TEST(ArrheniusRateConstant, CalculateWithSystem){
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(micm::System());
  EXPECT_NEAR(k, 0, 0.01);

  micm::ArrheniusRateConstant basic(1, 0, 0, 0, 0);
  k = basic.calculate(micm::System());
  EXPECT_NEAR(k, 1, 0.01);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  micm::ArrheniusRateConstant o1d(2.2e-10, 0, 0, 0, 0);
  k = o1d.calculate(micm::System());
  EXPECT_NEAR(k, 2.2e-10, 0.01);

  // O + HO2 -> OH + O2
  micm::ArrheniusRateConstant hox(3e-11, 0, -200, 0, 0);
  k = hox.calculate(micm::System());
  // at present the system cannot be used to determine the temperature so they should be different
  EXPECT_NE(k, 5.9e-11);

  // OH + HCl → H2O + Cl
  micm::ArrheniusRateConstant clox(2.6e-12, 0, -350, 0, 0);
  k = clox.calculate(micm::System());
  // at present the system cannot be used to determine the temperature so they should be different
  EXPECT_NE(k, 8e-13);
}

TEST(ArrheniusRateConstant, CalculateWithPrescribedArugments){
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(1.0, 1.0);
  EXPECT_NEAR(k, 0, 0.01);

  micm::ArrheniusRateConstant basic(1, 0, 0, 0, 0);
  k = basic.calculate(1.0, 1.0);
  EXPECT_NEAR(k, 1, 0.01);

  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  // O(1D) + H2O → OH + OH
  micm::ArrheniusRateConstant o1d(2.2e-10, 0, 0, 0, 0);
  k = o1d.calculate(1.0, 1.0);
  EXPECT_NEAR(k, 2.2e-10, 0.01);

  // O + HO2 -> OH + O2
  micm::ArrheniusRateConstant hox(3e-11, 0, 200, 0, 0);
  k = hox.calculate(298, 1.0);
  EXPECT_NEAR(k, 5.9e-11, 0.01);

  // OH + HCl → H2O + Cl
  micm::ArrheniusRateConstant clox(2.6e-12, 0, -350, 0, 0);
  k = clox.calculate(298, 100'000);
  EXPECT_NEAR(k, 8e-13, 0.01);
}