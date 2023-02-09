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

double call_fortran_rate(double temperature, double pressure, micm::ArrheniusRateConstantParameters params){
  return arrhenius_rate(
    temperature, pressure,
    params.A_, params.B_, params.C_, params.D_, params.E_
  );
}

TEST(ArrheniusRateRegressionTest, AllZeros){
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(1.0, 1.0);
  auto k_fortran = call_fortran_rate(1, 1, micm::ArrheniusRateConstantParameters());
  EXPECT_NEAR(k, k_fortran, 0.01);

}

TEST(ArrheniusRateRegressionTest, A1RestZero){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;
  micm::ArrheniusRateConstant basic(parameters);
  auto k = basic.calculate(1.0, 1.0);
  auto k_fortran = call_fortran_rate(1, 1, parameters);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

// values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
TEST(ArrheniusRateRegressionTest, O1D){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 2.2e-10;
  micm::ArrheniusRateConstant o1d(parameters);
  auto k = o1d.calculate(1.0, 1.0);
  auto k_fortran = call_fortran_rate(298, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

// O + HO2 -> OH + O2
TEST(ArrheniusRateRegressionTest, HOx){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  micm::ArrheniusRateConstant hox(parameters);
  auto k = hox.calculate(298, 1.0);
  auto k_fortran = call_fortran_rate(298, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

// OH + HCl → H2O + Cl
TEST(ArrheniusRateRegressionTest, ClOx){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  micm::ArrheniusRateConstant clox(parameters);
  auto k = clox.calculate(298, 100'000);
  auto k_fortran = call_fortran_rate(298, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

// k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
TEST(ArrheniusRateRegressionTest, N2_O1D_1){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 2.15e-11;
  parameters.C_ = 110;
  micm::ArrheniusRateConstant rate(parameters);
  auto k = rate.calculate(298, 100'000);
  auto k_fortran = call_fortran_rate(298, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

// k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
TEST(ArrheniusRateRegressionTest, O1D_O2_1){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 3.3e-11;
  parameters.C_ = 55;
  micm::ArrheniusRateConstant rate(parameters);
  auto k = rate.calculate(298, 100'000);
  auto k_fortran = call_fortran_rate(298, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 0.01);
}


// k_O_O3_1: O + O3 -> 2*O2
TEST(ArrheniusRateRegressionTest, O_O3_1){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 8e-12;
  parameters.C_ = -2060;
  micm::ArrheniusRateConstant rate(parameters);
  auto k = rate.calculate(298, 100'000);
  auto k_fortran = call_fortran_rate(298, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 0.01);
}

// k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
TEST(ArrheniusRateRegressionTest, M_O_O2_1){
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 6e-34;
  parameters.B_ = -2.4;
  micm::ArrheniusRateConstant rate(parameters);
  auto k = rate.calculate(273.14, 100'000);
  auto k_fortran = call_fortran_rate(273.14, 100'000, parameters);
  std::cout << k << " " << k_fortran << std::endl;
  EXPECT_EQ(k, k_fortran);
}