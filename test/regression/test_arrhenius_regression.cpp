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
  micm::ArrheniusRateConstant basic(1, 0, 0, 0, 0);
  auto k = basic.calculate(1.0, 1.0);
  auto k_fortran = arrhenius_rate(1, 1, 1, 0, 0, 0, 0);
  EXPECT_NEAR(k, 1, 0.01);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

TEST(ArrheniusRateRegressionTest, O1D){
  // values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
  micm::ArrheniusRateConstant o1d(2.2e-10, 0, 0, 0, 0);
  auto k = o1d.calculate(1.0, 1.0);
  auto k_fortran = arrhenius_rate(1, 1, 2.2e-10, 0, 0, 0, 0);
  EXPECT_NEAR(k, 2.2e-10, 0.01);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

TEST(ArrheniusRateRegressionTest, HOx){
  // O + HO2 -> OH + O2
  micm::ArrheniusRateConstant hox(3e-11, 0, 200, 0, 0);
  auto k = hox.calculate(298, 1.0);
  auto k_fortran = arrhenius_rate(298, 1, 3e-11, 0, 200, 0, 0);
  EXPECT_NEAR(k, 5.9e-11, 0.01);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

TEST(ArrheniusRateRegressionTest, ClOx){
  // OH + HCl → H2O + Cl
  micm::ArrheniusRateConstant clox(2.6e-12, 0, -350, 0, 0);
  auto k = clox.calculate(298, 100'000);
  auto k_fortran = arrhenius_rate(298, 100'000, 2.6e-12, 0, -350, 0, 0);
  EXPECT_NEAR(k, 8e-13, 0.01);
  EXPECT_NEAR(k, k_fortran, 0.01);
}