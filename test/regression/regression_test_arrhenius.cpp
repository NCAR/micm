#include <micm/process/arrhenius_rate_constant.hpp>

#include <gtest/gtest.h>

extern "C" double arrhenius_rate(
  double temperature, 
  double pressure, 
  double a, 
  double b, 
  double c, 
  double d, 
  double e
  );

TEST(ArrheniusRateRegressionTest, AllZeros){
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(1.0, 1.0);
  auto k_fortran = arrhenius_rate(1, 1, 0, 0, 0, 0, 0);
  EXPECT_NEAR(k, 0, 0.01);
  EXPECT_NEAR(k, k_fortran, 0.01);

}

TEST(ArrheniusRateRegressionTest, A1RestZero){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;
  micm::ArrheniusRateConstant basic(parameters);
  auto k = basic.calculate(1.0, 1.0);
  auto k_fortran = arrhenius_rate(1, 1, 1, 0, 0, 0, 0);
  EXPECT_NEAR(k, 1, 0.01);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

// values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
TEST(ArrheniusRateRegressionTest, O1D){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 2.2e-10;
  micm::ArrheniusRateConstant o1d(parameters);
  auto k = o1d.calculate(1.0, 1.0);
  auto k_fortran = arrhenius_rate(1, 1, 2.2e-10, 0, 0, 0, 0);
  EXPECT_NEAR(k, 2.2e-10, 0.01);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

// O + HO2 -> OH + O2
TEST(ArrheniusRateRegressionTest, HOx){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  micm::ArrheniusRateConstant hox(parameters);
  auto k = hox.calculate(298, 1.0);
  auto k_fortran = arrhenius_rate(298, 1, 3e-11, 0, 200, 0, 0);
  EXPECT_NEAR(k, 5.9e-11, 0.01);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

// OH + HCl → H2O + Cl
TEST(ArrheniusRateRegressionTest, ClOx){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  micm::ArrheniusRateConstant clox(parameters);
  auto k = clox.calculate(298, 100'000);
  auto k_fortran = arrhenius_rate(298, 100'000, 2.6e-12, 0, -350, 0, 0);
  EXPECT_NEAR(k, 8e-13, 0.01);
  EXPECT_NEAR(k, k_fortran, 0.01);
}